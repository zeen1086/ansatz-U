"""Trainable compact spatial kernels for the single-response V1 model."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from features import Structure
from symmetry import SPATIAL, finite, require


class KernelBasis:
    """D4h A1g/B1g nodal basis and exact full-grid quadratic field losses.

    Each radial scale supplies [1, z², x²+y², x²-y²]. Coordinates are
    normalized by the cutoff. Node values exclude the native cosine window.
    """

    def __init__(self, structure: Structure, grid: np.ndarray,
                 lengths=(.6, 1., 1.5, 2., 3., 4.), cutoff: float = 6.):
        structure.validate()
        self.grid = np.asarray(grid, dtype=int)
        require(self.grid.shape == (3,) and np.all(self.grid > 0)
                and np.array_equal(self.grid, grid), 'positive integer field grid required')
        self.cell = np.asarray(structure.cell, dtype=float).copy()
        require(np.allclose(self.cell, np.diag(np.diag(self.cell)), rtol=0, atol=1e-10)
                and np.all(np.diag(self.cell) > 0), 'axis-aligned positive kernel cell required')
        require(np.isclose(self.cell[0, 0], self.cell[1, 1], rtol=0, atol=1e-10)
                and self.grid[0] == self.grid[1], 'D4h requires equal x/y grid and lattice')
        require(cutoff == 6., 'native V1 export requires a 6 angstrom cutoff')
        self.cutoff = float(cutoff)
        self.lengths = np.asarray(lengths, dtype=float)
        require(self.lengths.ndim == 1 and len(self.lengths) > 0
                and np.all(np.isfinite(self.lengths)) and np.all(self.lengths > 0)
                and len(set(self.lengths)) == len(self.lengths), 'unique positive radial lengths required')
        self.ncoeff = self.nbasis = 4 * len(self.lengths)
        self.ni = len(structure.ni)
        self.spacing = np.diag(self.cell) / self.grid
        self.half = np.ceil(self.cutoff / self.spacing).astype(int)
        axes = [np.arange(-h, h+1, dtype=np.int32) for h in self.half]
        cube = np.stack(np.meshgrid(*axes, indexing='ij'), axis=-1).reshape(-1, 3)
        # Keep every vertex needed by trilinear interpolation inside the cutoff.
        radius = np.linalg.norm(cube * self.spacing, axis=1)
        self.nodes = np.ascontiguousarray(cube[radius <= self.cutoff + np.linalg.norm(self.spacing)])
        xyz = self.nodes * self.spacing
        r2 = np.sum(xyz * xyz, axis=1)
        x, y, z = (xyz / self.cutoff).T
        polynomial = np.column_stack((np.ones(len(x)), z*z, x*x+y*y, x*x-y*y))
        radial = np.exp(-r2[:, None] / self.lengths**2)
        self.values = np.ascontiguousarray((radial[:, :, None] * polynomial[:, None, :]).reshape(len(x), -1))
        self.characters = np.tile([1., 1., 1., -1.], len(self.lengths))
        self.mask = self.characters[:, None] == self.characters[None, :]
        self.representation = np.array([np.diag(np.where(self.characters > 0, 1.,
                                             1. if action[0, 0] != 0 else -1.))
                                        for action in SPATIAL])
        centers = structure.positions[structure.ni] / self.spacing
        require(np.allclose(centers, np.rint(centers), rtol=0, atol=1e-8),
                'Ni centers must coincide with field-grid nodes')
        centers = np.rint(centers).astype(int)
        radius = np.sqrt(r2)
        inside = radius < self.cutoff
        nodal = self.values[inside] * np.cos(np.pi*radius[inside, None]/(2*self.cutoff))**2
        self.design = np.zeros((int(np.prod(self.grid)), self.ni*self.nbasis))
        for site, center in enumerate(centers):
            indices = (self.nodes[inside] + center) % self.grid
            linear = indices[:, 0] + self.grid[0] * (indices[:, 1] + self.grid[1]*indices[:, 2])
            np.add.at(self.design[:, site*self.nbasis:(site+1)*self.nbasis], linear, nodal)
        self.mean = self.design.mean(axis=0)
        self.gram_B = self.design.T @ self.design / len(self.design)
        self.gram_v = self.gram_B - np.outer(self.mean, self.mean)

    def initial_parameters(self) -> dict[str, np.ndarray]:
        return dict(W_v=np.eye(self.ncoeff), W_B=np.eye(self.ncoeff), bias=np.zeros(self.nbasis))

    def project(self, parameters: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
        require(set(parameters) == {'W_v', 'W_B', 'bias'}, 'kernel parameter keys mismatch')
        result = {}
        for name in ('W_v', 'W_B'):
            value = np.asarray(parameters[name], dtype=float)
            require(value.shape == (self.nbasis, self.ncoeff), 'kernel matrix shape mismatch')
            finite(value, name)
            result[name] = value * self.mask
        bias = np.asarray(parameters['bias'], dtype=float)
        require(bias.shape == (self.nbasis,), 'kernel bias shape mismatch')
        finite(bias, 'kernel bias')
        result['bias'] = bias * (self.characters > 0)
        return result

    def compact(self, fields_eV: np.ndarray) -> dict[str, np.ndarray | float]:
        fields = np.asarray(fields_eV, dtype=float)
        require(fields.shape == (len(self.design), 4), 'field label shape mismatch')
        finite(fields, 'field labels')
        require(abs(fields[:, 0].mean()) < 1e-8, 'scalar labels must have full-cell mean zero')
        target = self.design.T @ fields / len(fields)
        # Centering D also removes the small floating-point mean in the label.
        target[:, 0] -= self.mean * fields[:, 0].mean()
        return dict(v=target[:, 0], B=target[:, 1:],
                    norm_v=float(np.mean(fields[:, 0]**2)),
                    norm_B=float(np.mean(np.sum(fields[:, 1:]**2, axis=1))))

    def _coefficients(self, p: np.ndarray, parameters: dict) -> tuple[np.ndarray, np.ndarray, dict]:
        require(p.shape == (self.ni, self.ncoeff, 4), 'physical model coefficient shape mismatch')
        finite(p, 'physical model coefficients')
        params = self.project(parameters)
        v = p[:, :, 0] @ params['W_v'].T + params['bias']
        b = np.einsum('nac,ba->nbc', p[:, :, 1:], params['W_B'])
        return v.reshape(-1), b.reshape(-1, 3), params

    def fields(self, p: np.ndarray, parameters: dict) -> np.ndarray:
        v, b, _ = self._coefficients(p, parameters)
        scalar = self.design @ v
        scalar -= scalar.mean()
        return np.column_stack((scalar, self.design @ b))

    def loss_and_grad(self, p: np.ndarray, parameters: dict, target: dict,
                      scalar_scale: float, magnetic_scale: float):
        require(np.isfinite(scalar_scale) and scalar_scale > 0
                and np.isfinite(magnetic_scale) and magnetic_scale > 0,
                'positive finite field loss scales required')
        v, b, params = self._coefficients(p, parameters)
        rv = self.gram_v @ v - target['v']
        rb = self.gram_B @ b - target['B']
        mse_v = float(v @ (rv-target['v']) + target['norm_v'])
        mse_B = float(np.sum(b*(rb-target['B'])) + target['norm_B'])
        gv = (2*rv/scalar_scale**2).reshape(self.ni, self.nbasis)
        gb = (2*rb/magnetic_scale**2).reshape(self.ni, self.nbasis, 3)
        gp = np.empty_like(p)
        gp[:, :, 0] = gv @ params['W_v']
        gp[:, :, 1:] = np.einsum('nbc,ba->nac', gb, params['W_B'])
        gradients = self.project(dict(W_v=gv.T @ p[:, :, 0],
                                     W_B=np.einsum('nbc,nac->ba', gb, p[:, :, 1:]),
                                     bias=gv.sum(axis=0)))
        loss = mse_v/scalar_scale**2 + mse_B/magnetic_scale**2
        return loss, gp, gradients, dict(v_mse_eV2=max(0., mse_v), B_vector_mse_eV2=max(0., mse_B))

    def export(self, path: Path, parameters: dict) -> None:
        params = self.project(parameters)
        path = Path(path)
        names = ('C_v_local.npy', 'C_B_local.npy', 'local_grid.npz', 'KERNELS.json')
        require(not any((path/name).exists() for name in names), 'refusing to replace spatial kernels')
        path.mkdir(parents=True, exist_ok=True)
        cv = np.column_stack((self.values @ params['bias'], self.values @ params['W_v']))
        cb = self.values @ params['W_B']
        np.save(path/'C_v_local.npy', np.ascontiguousarray(cv), allow_pickle=False)
        np.save(path/'C_B_local.npy', np.ascontiguousarray(cb), allow_pickle=False)
        np.savez(path/'local_grid.npz', nodes=self.nodes, spacing_A=self.spacing, half=self.half)
        record = dict(schema='V1_D4H_NODAL_KERNELS_V1', cutoff_A=self.cutoff,
                      scalar_bias_coefficient_per_Ni=1., radial_lengths_A=self.lengths.tolist(),
                      polynomial=['1', '(z/cutoff)^2', '(x^2+y^2)/cutoff^2', '(x^2-y^2)/cutoff^2'],
                      representation='A1g,A1g,A1g,B1g per radial length',
                      node_values='pre-window; native evaluator applies cos(pi*r/(2*cutoff))^2',
                      parameterization='same-representation linear kernel mixtures; new single-response fit')
        (path/'KERNELS.json').write_text(json.dumps(record, indent=2) + '\n')
