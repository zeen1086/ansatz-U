#!/usr/bin/env python3
"""Fit v/B/full-complex-N labels with one explicit occupation correction."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import tomllib

import numpy as np
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE/'local'))
from features import FeatureSpec, features
from kernels import KernelBasis
from model import LocalModel
from symmetry import to_channels, to_matrix
from workflow import SINGLE_RESPONSE_SCHEMA, save_model
from training_data import initial_reference, load_dataset

SCHEMA = 'LANIO2_SINGLE_RESPONSE_FIT_V1'
DEFAULTS = dict(
    features=dict(radial_lengths_A=[2., 4.], scalar_degree=4, vector_degree=3),
    kernels=dict(radial_lengths_A=[.6, 1., 1.5, 2., 3., 4.]),
    loss=dict(scalar_scale_eV=1., magnetic_scale_eV=1., occupation_scale=.1, regularization=1e-8),
    optimizer=dict(max_iterations=200, gradient_tolerance=1e-6, loss_tolerance=1e-12,
                   checkpoint_every=10, seed=0))


def digest(path: Path) -> str:
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def read_config(path: Path) -> dict:
    with path.open('rb') as stream:
        raw = tomllib.load(stream)
    if set(raw) - {'data', *DEFAULTS}:
        raise ValueError('Unknown fit configuration section')
    data = raw.get('data', {})
    if set(data) - {'labels', 'weights'} or not isinstance(data.get('labels'), list) or not data['labels']:
        raise ValueError('[data] must provide a nonempty labels list and optional weights')
    config = {'data': dict(labels=[str((path.parent/label).resolve()) for label in data['labels']],
                           weights=data.get('weights'))}
    for section, defaults in DEFAULTS.items():
        override = raw.get(section, {})
        if set(override) - set(defaults):
            raise ValueError(f'Unknown {section} setting')
        config[section] = defaults | override
    for section in ('features', 'kernels'):
        lengths = np.asarray(config[section]['radial_lengths_A'], dtype=float)
        if lengths.ndim != 1 or not len(lengths) or not np.isfinite(lengths).all() or np.any(lengths <= 0):
            raise ValueError('Positive finite radial lengths required')
    if config['features']['scalar_degree'] not in (1, 4) or config['features']['vector_degree'] not in (1, 3, 4):
        raise ValueError('Supported feature degrees: scalar 1/4, vector 1/3/4')
    for key, value in config['loss'].items():
        if not np.isfinite(value) or value < 0 or (key != 'regularization' and value == 0):
            raise ValueError('Positive loss scales and nonnegative regularization required')
    optimizer = config['optimizer']
    for key in ('max_iterations', 'checkpoint_every', 'seed'):
        if type(optimizer[key]) is not int or optimizer[key] < (0 if key == 'seed' else 1):
            raise ValueError(f'Invalid optimizer {key}')
    for key in ('gradient_tolerance', 'loss_tolerance'):
        if not np.isfinite(optimizer[key]) or optimizer[key] <= 0:
            raise ValueError(f'Invalid optimizer {key}')
    return config


class FitProblem:
    """Exact full-grid quadratic field objective plus the current model gradient."""

    def __init__(self, dataset, config: dict):
        self.dataset, self.config = dataset, config
        self.basis = KernelBasis(dataset.samples[0].structure, dataset.grid,
                                 lengths=config['kernels']['radial_lengths_A'])
        nc = self.basis.ncoeff
        spec = FeatureSpec(6., tuple(config['features']['radial_lengths_A']),
                           scalar_degree=config['features']['scalar_degree'],
                           vector_degree=config['features']['vector_degree'])
        reference = initial_reference(dataset, nc)
        self.model = LocalModel(spec, np.ones((nc, 2)), reference,
                                p_representation=self.basis.representation,
                                seed=config['optimizer']['seed'])
        self.model.v_channels = nc
        raw = [features(sample.structure, spec) for sample in dataset.samples]
        for key in self.model.feature_scales:
            square = sum(w*np.mean(phi[key]**2, axis=(0, 1, 3)) for w, phi in zip(dataset.weights, raw))
            rms = np.sqrt(square)
            self.model.feature_scales[key] = np.where(rms > 1e-10, rms, 1.)
        self.phi = [self.model.prepare_features(phi) for phi in raw]
        nr, _ = reference.values(np.array([[0., 0., 1.]]), self.model.scale)
        eigenvalues, vectors = np.linalg.eigh(nr[0])
        # Only the initial logit needs finite endpoints; labels are never clipped.
        eigenvalues = np.clip(eigenvalues, 1e-6, 1-1e-6)
        qref = (vectors * np.log(eigenvalues/(1-eigenvalues))) @ vectors.conj().T
        channels = to_channels(qref)
        for sector in ('se', 'vo'):
            self.model.params['b_Q_' + sector] = channels[sector][:, 0 if sector == 'se' else 2].copy()
        for key in self.model.params:
            if key.startswith(('Q_', 'p_')):
                self.model.params[key] *= .05
        self.kernel_parameters = self.basis.initial_parameters()
        self.targets = [self.basis.compact(sample.field_block()) for sample in dataset.samples]
        self.parameters = {'model/' + k: v for k, v in self.model.params.items()}
        self.parameters.update({'kernel/' + k: v for k, v in self.kernel_parameters.items()})
        self.keys = sorted(self.parameters)
        self.shapes = {key: self.parameters[key].shape for key in self.keys}
        self.sizes = {key: self.parameters[key].size for key in self.keys}

    def pack(self, values: dict | None = None) -> np.ndarray:
        values = self.parameters if values is None else values
        return np.concatenate([values[key].ravel() for key in self.keys])

    def unpack(self, vector: np.ndarray) -> None:
        if vector.shape != (sum(self.sizes.values()),) or not np.isfinite(vector).all():
            raise ValueError('Invalid optimizer parameter vector')
        cursor = 0
        for key in self.keys:
            value = vector[cursor:cursor+self.sizes[key]].reshape(self.shapes[key]).copy()
            self.parameters[key] = value
            owner, name = key.split('/')
            (self.model.params if owner == 'model' else self.kernel_parameters)[name] = value
            cursor += self.sizes[key]

    def evaluate(self, vector: np.ndarray):
        self.unpack(vector)
        loss = self.config['loss']
        grads = {key: np.zeros(shape) for key, shape in self.shapes.items()}
        objective, reports = 0., []
        for weight, sample, phi, target in zip(self.dataset.weights, self.dataset.samples, self.phi, self.targets):
            pred = self.model.predict(sample.structure, prepared_features=phi, normalized=True)
            field_loss, gp, gkernel, metrics = self.basis.loss_and_grad(
                self.model.physical_p(pred), self.kernel_parameters, target,
                loss['scalar_scale_eV'], loss['magnetic_scale_eV'])
            delta = pred.N - sample.N
            n_mse = float(np.vdot(delta, delta).real / delta.size)
            gn = 2*delta/(delta.size*loss['occupation_scale']**2)
            gmodel = self.model.backward(pred, gp*self.model.scale, gn)
            objective += weight*(field_loss + n_mse/loss['occupation_scale']**2)
            for key, value in gmodel.items():
                grads['model/'+key] += weight*value
            for key, value in gkernel.items():
                grads['kernel/'+key] += weight*value
            reports.append(dict(name=sample.name, weight=float(weight),
                                v_RMS_eV=float(np.sqrt(metrics['v_mse_eV2'])),
                                B_vector_RMS_eV=float(np.sqrt(metrics['B_vector_mse_eV2'])),
                                B_component_RMS_eV=float(np.sqrt(metrics['B_vector_mse_eV2']/3)),
                                N_element_RMS=float(np.sqrt(n_mse))))
        regularization = loss['regularization']*float(vector @ vector)/len(vector)
        gradient = self.pack(grads) + 2*loss['regularization']*vector/len(vector)
        objective += regularization
        if not np.isfinite(objective) or not np.isfinite(gradient).all():
            raise ValueError('Nonfinite training objective or gradient')
        return objective, gradient, dict(samples=reports, regularization=regularization)

    def objective(self, vector: np.ndarray):
        try:
            value, gradient, _ = self.evaluate(vector)
            return value, gradient
        except ValueError as exc:
            if 'sigmoid saturated' not in str(exc) and 'sigmoid lost strict bounds' not in str(exc):
                raise
            # Reject an out-of-domain line-search trial without clipping predictions.
            return np.inf, np.zeros_like(vector)


def atomic_json(path: Path, record: dict) -> None:
    temporary = path.with_suffix(path.suffix+'.partial')
    temporary.write_text(json.dumps(record, indent=2, allow_nan=False)+'\n')
    temporary.replace(path)


def identity(config: dict, dataset) -> dict:
    code = [HERE/'fit.py', HERE/'training_data.py',
            *(HERE/'local'/name for name in ('features.py', 'model.py', 'symmetry.py', 'workflow.py', 'kernels.py'))]
    return dict(schema=SCHEMA, model_schema=SINGLE_RESPONSE_SCHEMA, config=config,
                dataset=dataset.identity, code={str(p.relative_to(HERE)): digest(p) for p in code})


def save_checkpoint(path: Path, problem: FitProblem, iteration: int, binding: dict) -> None:
    vector = problem.pack()
    metadata = dict(identity=binding, iteration=iteration, keys=problem.keys,
                    parameters_sha256=hashlib.sha256(vector.tobytes()).hexdigest())
    temporary = path.with_suffix('.partial')
    with temporary.open('wb') as stream:
        np.savez(stream, parameters=vector, metadata=json.dumps(metadata, allow_nan=False))
    temporary.replace(path)


def read_checkpoint(path: Path, problem: FitProblem, binding: dict) -> int:
    with np.load(path, allow_pickle=False) as data:
        meta = json.loads(str(data['metadata']))
        if meta.get('identity') != binding or meta.get('keys') != problem.keys:
            raise ValueError('Checkpoint model schema, labels, configuration or code changed')
        vector = data['parameters'].copy()
    if hashlib.sha256(vector.tobytes()).hexdigest() != meta.get('parameters_sha256'):
        raise ValueError('Checkpoint parameter checksum mismatch')
    if type(meta.get('iteration')) is not int or meta['iteration'] < 0:
        raise ValueError('Invalid checkpoint iteration')
    problem.unpack(vector)
    return meta['iteration']


def resume_history(path: Path, iteration: int) -> None:
    """Remove iterations whose parameters did not reach the last checkpoint."""
    if not path.exists():
        if iteration:
            raise ValueError('Checkpoint iteration has no corresponding history')
        return
    rows = [json.loads(line) for line in path.read_text().splitlines()]
    if [row.get('iteration') for row in rows] != list(range(1, len(rows)+1)) or len(rows) < iteration:
        raise ValueError('Training history and checkpoint iterations disagree')
    if len(rows) > iteration:
        temporary = path.with_suffix('.partial')
        temporary.write_text(''.join(json.dumps(row, allow_nan=False)+'\n' for row in rows[:iteration]))
        temporary.replace(path)


def train(config_path: Path, output: Path, resume: bool = False) -> dict:
    config_path, output = Path(config_path).resolve(), Path(output).absolute()
    config = read_config(config_path)
    if output.is_symlink() or (output.exists() and not resume):
        raise FileExistsError(f'Refusing to replace existing training output: {output}')
    if resume and not output.is_dir():
        raise ValueError('Resume requires an existing training directory')
    dataset = load_dataset([Path(p) for p in config['data']['labels']], config['data']['weights'])
    binding = identity(config, dataset)
    if resume and json.loads((output/'RUN.json').read_text()) != binding:
        raise ValueError('Resume requires unchanged single-response code, settings and labels')
    problem = FitProblem(dataset, config)
    iteration = read_checkpoint(output/'checkpoint.npz', problem, binding) if resume else 0
    if resume:
        resume_history(output/'history.jsonl', iteration)
    initial, _, _ = problem.evaluate(problem.pack())
    if not resume:
        output.mkdir(parents=True, exist_ok=False)
        atomic_json(output/'RUN.json', binding)
        save_checkpoint(output/'checkpoint.npz', problem, iteration, binding)
    history = output/'history.jsonl'
    options = config['optimizer']

    def accepted(vector):
        nonlocal iteration
        iteration += 1
        value, gradient, _ = problem.evaluate(vector)
        with history.open('a') as stream:
            stream.write(json.dumps(dict(iteration=iteration, loss=value,
                                        gradient_inf=float(np.max(abs(gradient)))), allow_nan=False)+'\n')
        if iteration % options['checkpoint_every'] == 0:
            save_checkpoint(output/'checkpoint.npz', problem, iteration, binding)

    result = minimize(problem.objective, problem.pack(), jac=True, method='L-BFGS-B', callback=accepted,
                      options=dict(maxiter=options['max_iterations'], gtol=options['gradient_tolerance'],
                                   ftol=options['loss_tolerance'], maxls=40))
    value, gradient, metrics = problem.evaluate(result.x)
    for sample in dataset.identity['samples']:
        directory = Path(sample['path'])
        if (digest(directory/'LABELS.json') != sample['manifest_sha256']
                or any(digest(directory/name) != expected for name, expected in sample['files'].items())):
            raise ValueError('Training labels changed during optimization; no model published')
    save_checkpoint(output/'checkpoint.npz', problem, iteration, binding)
    bundle = output/f'model_{iteration:06d}'
    parameter_hash = hashlib.sha256(problem.pack().tobytes()).hexdigest()
    report = dict(schema=SCHEMA, model_schema=SINGLE_RESPONSE_SCHEMA, parameters_sha256=parameter_hash,
                  prediction_rule='N0=sigmoid(Q0); N1=sigmoid(Q0+A(N0-Nref)); p=p0+K(N1-Nref)',
                  iterations=iteration, new_iterations=int(result.nit), initial_loss=initial, final_loss=value,
                  optimizer_success=bool(result.success), optimizer_status=int(result.status),
                  optimizer_message=str(result.message), gradient_inf=float(np.max(abs(gradient))),
                  gradient_converged=bool(np.max(abs(gradient)) <= options['gradient_tolerance']),
                  metrics=metrics, bundle=str(bundle), physical_release=False,
                  evaluation_scope='training labels only; no held-out or material-accuracy claim',
                  restart='parameter restart; L-BFGS history is not restored')
    if not bundle.exists():
        with tempfile.TemporaryDirectory(prefix='.model-', dir=output) as temporary:
            stage = Path(temporary)
            save_model(problem.model, stage/'model', training_record=report)
            problem.basis.export(stage, problem.kernel_parameters)
            sample = dataset.samples[0].structure
            np.savez(stage/'structure.npz', cell=sample.cell, positions=sample.positions,
                     species=sample.species, ni=sample.ni, spins=sample.spins)
            report['bundle_files'] = {str(p.relative_to(stage)): digest(p) for p in stage.rglob('*') if p.is_file()}
            atomic_json(stage/'FIT_REPORT.json', report)
            stage.rename(bundle)
    else:
        previous = json.loads((bundle/'FIT_REPORT.json').read_text())
        if not resume or previous.get('parameters_sha256') != parameter_hash:
            raise FileExistsError(f'Refusing to replace a different model bundle: {bundle}')
        required = {'C_v_local.npy', 'C_B_local.npy', 'local_grid.npz', 'KERNELS.json',
                    'model/model.json', 'model/parameters.npz', 'structure.npz'}
        hashes = previous.get('bundle_files', {})
        if set(hashes) != required or any(digest(bundle/name) != hashes[name] for name in required):
            raise ValueError('Existing model bundle checksum mismatch')
        report['bundle_files'] = hashes
    atomic_json(output/'RESULT.json', report)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--resume', action='store_true', help='Restart parameters with identical labels/settings/code')
    args = parser.parse_args()
    report = train(args.config, args.output, args.resume)
    print(json.dumps(report, indent=2))
    if not report['optimizer_success']:
        raise SystemExit(2)


if __name__ == '__main__':
    main()
