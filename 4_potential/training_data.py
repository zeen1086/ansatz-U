"""Validated V1 label packages and references for single-response fitting."""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import sys
from typing import Sequence

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent/'local'))
from features import Structure
from model import Reference
from symmetry import REP, require, to_channels

RY_EV = 13.605693122994
SCHEMA = 'LANIO2_V1_SCF_LABELS_V1'
FILES = ('scalar.npy', 'B.npy', 'N.npy', 'moments.npy', 'atom_ids.npy', 'grid.npy', 'structure.npz')
UNITS = dict(scalar='Ry', B='Ry', N='dimensionless', moments='mu_B', cell='angstrom', positions='angstrom')
ORBITALS = ['dz2', 'dxz', 'dyz', 'dx2-y2', 'dxy']


def _sha(path: Path) -> str:
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


@dataclass(frozen=True)
class Sample:
    name: str
    directory: Path
    structure: Structure
    N: np.ndarray
    moments: np.ndarray
    scalar_Ry: np.ndarray
    B_Ry: np.ndarray
    manifest: dict

    def field_block(self, start: int = 0, stop: int | None = None) -> np.ndarray:
        """Materialize only the requested scalar/B grid block, in eV."""
        return np.column_stack((self.scalar_Ry[start:stop], self.B_Ry[start:stop])) * RY_EV


@dataclass(frozen=True)
class Dataset:
    samples: tuple[Sample, ...]
    grid: np.ndarray
    weights: np.ndarray
    identity: dict


def _array(directory: Path, name: str, shape: tuple, dtype) -> np.ndarray:
    value = np.load(directory/name, mmap_mode='r', allow_pickle=False)
    require(value.shape == shape and value.dtype == np.dtype(dtype), f'{name}: shape or dtype mismatch')
    require(np.isfinite(value).all(), f'{name}: nonfinite labels')
    return value


def _load_sample(directory: Path) -> tuple[Sample, np.ndarray, dict]:
    manifest_path = directory/'LABELS.json'
    manifest = json.loads(manifest_path.read_text())
    require(manifest.get('schema') == SCHEMA, 'unsupported label schema')
    require(manifest.get('units') == UNITS, 'label units mismatch')
    require(manifest.get('grid_order') == 'x fastest; ix + nx*(iy + ny*iz)', 'unsupported label grid order')
    require(manifest.get('field_convention') == 'v*I + Bx*sigma_x + By*sigma_y + Bz*sigma_z',
            'unsupported field convention')
    require(manifest.get('scalar_components') == 'local ionic + Hartree + XC; full-cell mean subtracted',
            'unsupported scalar gauge/components')
    require(manifest.get('B_components') == 'XC Pauli coefficients; no explicit constraint field or Hubbard term',
            'unsupported magnetic-field components')
    hashes = manifest.get('files', {})
    require(set(hashes) == set(FILES), 'incomplete label file manifest')
    for name in FILES:
        require(isinstance(hashes[name], str) and len(hashes[name]) == 64, 'invalid label hash')
        require(_sha(directory/name) == hashes[name], f'label checksum mismatch: {name}')

    grid = np.load(directory/'grid.npy', allow_pickle=False)
    require(grid.shape == (3,) and np.issubdtype(grid.dtype, np.integer) and np.all(grid > 0),
            'invalid label grid')
    require(np.array_equal(np.asarray(manifest.get('grid')), grid), 'manifest grid differs from grid.npy')
    ngrid = int(np.prod(grid, dtype=np.int64))
    require(ngrid > 0, 'invalid label grid product')
    scalar = _array(directory, 'scalar.npy', (ngrid,), np.float64)
    magnetic = _array(directory, 'B.npy', (ngrid, 3), np.float64)
    require(abs(float(scalar.mean())) < 1e-8, 'scalar labels must have zero full-cell mean')
    n = _array(directory, 'N.npy', (4, 10, 10), np.complex128)
    require(np.max(abs(n - n.conj().swapaxes(-1, -2))) < 1e-8, 'N labels must be Hermitian')
    eigenvalues = np.linalg.eigvalsh(n)
    require(eigenvalues.min() >= -1e-7 and eigenvalues.max() <= 1 + 1e-7,
            'N labels outside occupation bounds')
    moments = _array(directory, 'moments.npy', (4, 3), np.float64)
    amplitudes = np.linalg.norm(moments, axis=1)
    require(np.all(amplitudes >= 1e-8), 'zero Ni moment has no spin direction')
    atom_ids = np.load(directory/'atom_ids.npy', allow_pickle=False)
    require(atom_ids.shape == (4,) and np.issubdtype(atom_ids.dtype, np.integer), 'invalid atom_ids')
    require(np.array_equal(np.asarray(manifest.get('atom_ids')), atom_ids), 'manifest atom order mismatch')
    with np.load(directory/'structure.npz', allow_pickle=False) as values:
        require(set(values.files) == {'cell', 'positions', 'species', 'ni', 'spins'}, 'structure keys mismatch')
        structure = Structure(**{key: values[key].copy() for key in values.files})
    structure.validate()
    require(len(structure.positions) == 16 and len(structure.ni) == 4, 'require four-Ni La4Ni4O8 labels')
    require([np.count_nonzero(structure.species == z) for z in (57, 28, 8)] == [4, 4, 8],
            'label stoichiometry mismatch')
    require(np.array_equal(atom_ids, structure.ni + 1), 'N/moment atom order differs from the structure')
    require(np.allclose(structure.spins, moments/amplitudes[:, None], rtol=0, atol=1e-10),
            'structure spins differ from actual moment directions')
    cell = structure.cell
    require(np.allclose(cell, np.diag(np.diag(cell)), atol=1e-8, rtol=0)
            and np.all(np.diag(cell) > 0) and abs(cell[0, 0]-cell[1, 1]) < 1e-8,
            'require the axis-aligned tetragonal V1 cell')
    registered = np.array([[.25, .25, .5], [.75, .25, .5], [.25, .75, .5], [.75, .75, .5]])
    fractional = structure.positions[structure.ni]/np.diag(cell)
    matches = np.all(np.isclose(fractional[:, None], registered[None], atol=1e-8, rtol=0), axis=2)
    require(np.all(matches.sum(axis=0) == 1) and np.all(matches.sum(axis=1) == 1),
            'Ni positions differ from the registered V1 layout')
    occupation = manifest.get('occupation', {})
    require(occupation.get('functional') == 'PBE' and occupation.get('U_eV') == 1.
            and occupation.get('projector') == 'ortho-atomic', 'unsupported occupation provenance')
    require(occupation.get('orbital_order') == ORBITALS
            and occupation.get('spin_order') == 'up five orbitals, down five orbitals',
            'occupation orbital/spin order mismatch')
    identity = dict(path=str(directory), manifest_sha256=_sha(manifest_path), files=hashes)
    return Sample(directory.name, directory, structure, n, moments, scalar, magnetic, manifest), grid, identity


def load_dataset(paths: Sequence[Path], weights: Sequence[float] | None = None) -> Dataset:
    """Load one immutable, common-geometry set; weights come from the fit config."""
    directories = tuple(Path(path).resolve() for path in paths)
    require(len(directories) > 0 and len(set(directories)) == len(directories),
            'require a nonempty list of distinct label directories')
    require(len({directory.name for directory in directories}) == len(directories),
            'label directory names must be unique')
    weight = np.ones(len(directories)) if weights is None else np.asarray(weights, dtype=float)
    require(weight.shape == (len(directories),) and np.isfinite(weight).all() and np.all(weight > 0),
            'require one finite positive weight per label package')
    weight = weight/weight.max()
    weight /= weight.sum()
    samples, identities, common_grid = [], [], None
    for directory in directories:
        sample, grid, identity = _load_sample(directory)
        if samples:
            reference = samples[0].structure
            require(np.array_equal(grid, common_grid), 'training grids differ between configurations')
            for name in ('cell', 'positions'):
                require(np.allclose(getattr(sample.structure, name), getattr(reference, name),
                                    rtol=0, atol=1e-9), f'training {name} differ between configurations')
            for name in ('species', 'ni'):
                require(np.array_equal(getattr(sample.structure, name), getattr(reference, name)),
                        f'training {name} differ between configurations')
        else:
            common_grid = grid.copy()
        samples.append(sample)
        identities.append(identity)
    weight.flags.writeable = False
    common_grid.flags.writeable = False
    return Dataset(tuple(samples), common_grid, weight,
                   dict(schema=SCHEMA, samples=identities, weights=weight.tolist()))


def initial_reference(dataset: Dataset, ncoeff: int) -> Reference:
    """Project weighted measured N onto the D4h scalar and spin-aligned reference."""
    require(isinstance(ncoeff, int) and not isinstance(ncoeff, bool) and ncoeff > 0,
            'positive integer field coefficient count required')
    scalar, spin = np.zeros(15), np.zeros(15)
    for weight, sample in zip(dataset.weights, dataset.samples):
        channels = to_channels(sample.N)
        scalar += weight*channels['se'][:, :, 0].mean(axis=0)
        spin += weight*np.einsum('nac,nc->a', channels['vo'], sample.structure.spins)/len(sample.structure.ni)
    scalar = np.einsum('gab,b->a', REP['se'], scalar)/len(REP['se'])
    spin = np.einsum('gab,b->a', REP['vo'], spin)/len(REP['vo'])
    identity = hashlib.sha256(json.dumps(dataset.identity, sort_keys=True).encode()).hexdigest()
    return Reference(scalar, spin, np.zeros(ncoeff), np.zeros(ncoeff), 'label-set-sha256:' + identity)
