"""Full-grid losses, spatial covariance and native export for new V1 kernels."""
from itertools import product
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'4_potential/local'))
sys.path.insert(0, str(ROOT))
from features import FeatureSpec, Structure
from kernels import KernelBasis
from model import LocalModel, Reference
from symmetry import SPATIAL
from workflow import save_model
from local_vfunc import export


class KernelTests(unittest.TestCase):
    def setUp(self):
        self.structure = Structure(np.diag([8., 8., 4.]),
                                   np.array([[2., 2., 2.], [6., 2., 2.], [2., 6., 2.], [6., 6., 2.]]),
                                   np.full(4, 28), np.arange(4), np.tile([0., 0., 1.], (4, 1)))
        self.basis = KernelBasis(self.structure, np.array([8, 8, 4]), lengths=(1.2, 2.4))
        rng = np.random.default_rng(173)
        self.p = rng.normal(scale=.4, size=(4, self.basis.ncoeff, 4))
        self.parameters = self.basis.project(dict(W_v=rng.normal(size=(8, 8)),
                                                  W_B=rng.normal(size=(8, 8)), bias=rng.normal(size=8)))

    def test_compact_loss_and_both_gradients_match_grid_objective(self):
        b = self.basis
        rng = np.random.default_rng(15)
        fields = rng.normal(scale=.2, size=(len(b.design), 4))
        fields[:, 0] -= fields[:, 0].mean()
        target = b.compact(fields)
        loss, gp, grads, metrics = b.loss_and_grad(self.p, self.parameters, target, .3, .4)
        residual = b.fields(self.p, self.parameters) - fields
        expected_v = np.mean(residual[:, 0]**2)
        expected_B = np.mean(np.sum(residual[:, 1:]**2, axis=1))
        self.assertAlmostEqual(loss, expected_v/.3**2 + expected_B/.4**2, places=11)
        self.assertAlmostEqual(metrics['v_mse_eV2'], expected_v, places=13)
        self.assertAlmostEqual(metrics['B_vector_mse_eV2'], expected_B, places=13)
        h = 1e-5
        dp = rng.normal(size=self.p.shape)
        forward = b.loss_and_grad(self.p + h*dp, self.parameters, target, .3, .4)[0]
        backward = b.loss_and_grad(self.p - h*dp, self.parameters, target, .3, .4)[0]
        np.testing.assert_allclose(np.vdot(gp, dp), (forward-backward)/(2*h), rtol=1e-8, atol=1e-8)
        for key in self.parameters:
            direction = rng.normal(size=self.parameters[key].shape)
            plus = {k: v.copy() for k, v in self.parameters.items()}
            minus = {k: v.copy() for k, v in self.parameters.items()}
            plus[key] += h*direction; minus[key] -= h*direction
            forward = b.loss_and_grad(self.p, plus, target, .3, .4)[0]
            backward = b.loss_and_grad(self.p, minus, target, .3, .4)[0]
            np.testing.assert_allclose(np.vdot(grads[key], direction), (forward-backward)/(2*h),
                                       rtol=1e-8, atol=1e-8)
        self.assertTrue(np.all(grads['W_v'][~b.mask] == 0))
        self.assertTrue(np.all(grads['bias'][b.characters < 0] == 0))

    def test_periodic_sum_and_d4h_nodal_covariance(self):
        b = self.basis
        lookup = {tuple(node): i for i, node in enumerate(b.nodes)}
        for action, representation in zip(SPATIAL, b.representation):
            rows = [lookup[tuple(node)] for node in b.nodes @ action.T]
            np.testing.assert_allclose(b.values[rows], b.values @ representation.T, atol=2e-16)
        flat = np.arange(len(b.design))
        points = np.column_stack((flat % 8, flat//8 % 8, flat//64)) * b.spacing
        oracle = np.zeros_like(b.design)
        for site, center in enumerate(self.structure.positions):
            for image in product(range(-2, 3), range(-2, 3), range(-3, 4)):
                d = points - center + np.array(image) @ self.structure.cell
                radius = np.linalg.norm(d, axis=1)
                inside = radius < 6.
                if not np.any(inside):
                    continue
                xyz = d[inside]/6.
                poly = np.column_stack((np.ones(inside.sum()), xyz[:, 2]**2,
                                        xyz[:, 0]**2+xyz[:, 1]**2, xyz[:, 0]**2-xyz[:, 1]**2))
                radial = np.exp(-radius[inside, None]**2 / b.lengths**2)
                values = (radial[:, :, None]*poly[:, None, :]).reshape(inside.sum(), -1)
                values *= np.cos(np.pi*radius[inside, None]/12.)**2
                oracle[inside, site*8:(site+1)*8] += values
        np.testing.assert_allclose(b.design, oracle, atol=1e-15)
        fields = b.fields(self.p, self.parameters)
        self.assertLess(abs(fields[:, 0].mean()), 1e-14)
        self.assertGreater(np.linalg.norm(fields[:, 1:].mean(axis=0)), 1e-3)

    def test_export_matches_native_interpolation_and_z_images(self):
        compiler = shutil.which('gfortran')
        if compiler is None:
            self.skipTest('gfortran unavailable for native kernel interpolation check')
        b = self.basis
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); source = root/'source'
            b.export(source, self.parameters)
            model = LocalModel(FeatureSpec(6., (2.,)), np.ones((8, 2)),
                               Reference(np.zeros(15), np.zeros(15), np.zeros(8), np.zeros(8), 'synthetic'),
                               p_representation=b.representation)
            save_model(model, source/'model')
            structure = root/'structure.npz'
            np.savez(structure, **vars(self.structure))
            export(source, ROOT/'4_potential/local', structure, root/'native', 0.)
            command = [compiler, '-O0', '-J'+str(root), '-I'+str(root),
                       str(ROOT/'common/numerics.f90'), str(ROOT/'common/npy_io.f90'),
                       str(ROOT/'6_RVR/local_vfunc.f90'), str(ROOT/'tests/sample_local_vfunc.f90'),
                       '-o', str(root/'sample')]
            command += ['-framework', 'Accelerate'] if sys.platform == 'darwin' else ['-lblas']
            subprocess.run(command, cwd=root, check=True, capture_output=True, text=True)
            points = np.array([[0., 0., 0.], [.2, -.7, 1.3], [5.9, .1, 0.], [6.1, 0., 0.], [1., 2., -4.6]])
            np.savetxt(root/'points', points, header=str(len(points)), comments='')
            subprocess.run([str(root/'sample'), str(root/'native'), str(root/'points'), str(root/'values')],
                           check=True, capture_output=True, text=True)
            actual = np.fromfile(root/'values', '<f8').reshape(len(points), 2, 17)
            values = np.column_stack((np.load(source/'C_v_local.npy'), np.load(source/'C_B_local.npy')))
            lookup = {tuple(node): i for i, node in enumerate(b.nodes)}

            def local(point):
                radius = np.linalg.norm(point)
                if radius >= 6.:
                    return np.zeros(17)
                pos = point/b.spacing; low = np.floor(pos).astype(int); fraction = pos-low
                value = np.zeros(17)
                for shift in product((0, 1), repeat=3):
                    row = lookup.get(tuple(low+shift))
                    if row is not None:
                        weight = np.prod(np.where(shift, fraction, 1-fraction))
                        value += weight*values[row]
                return value*np.cos(np.pi*radius/12.)**2

            for i, point in enumerate(points):
                np.testing.assert_allclose(actual[i, 0], local(point), atol=1e-14)
                plane = sum((local(point+[0, 0, 4*z]) for z in range(-4, 5)), np.zeros(17))
                np.testing.assert_allclose(actual[i, 1], plane, atol=1e-14)
            with self.assertRaisesRegex(ValueError, 'refusing to replace'):
                b.export(source, self.parameters)


if __name__ == '__main__':
    unittest.main()
