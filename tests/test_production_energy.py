"""Small numerical regressions for the imported production kernels and file boundary."""
from pathlib import Path
import json
import os
import struct
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'7_spin'))
import complete_readout
import density_map
import kz_production
import native_input
import summarize_energy


class ProductionEnergyTest(unittest.TestCase):
    def test_external_config_and_reused_case_paths(self):
        repo = Path(__file__).resolve().parents[1]
        config = json.loads((repo/'examples/7_spin/production.json').read_text())
        with tempfile.TemporaryDirectory() as tmp:
            config['input_root'] = str(Path(tmp)/'inputs')
            config['complete_energy']['work'] = str(Path(tmp)/'output')
            config['complete_energy']['reuse_cases'] = {'4': str(Path(tmp)/'existing')}
            path = Path(tmp)/'run.json'
            path.write_text(json.dumps(config))
            command = ('import complete_energy as c; '
                       'print(c.ROOT); print(c.oldcase(4, 1)); print(c.oldcase(6, 1))')
            env = dict(os.environ, ANSATZ_RUN_CONFIG=str(path),
                       PYTHONPATH=str(repo/'7_spin'))
            result = subprocess.run([sys.executable, '-B', '-c', command], cwd=tmp,
                                    env=env, capture_output=True, text=True, check=True)
            self.assertEqual(result.stdout.splitlines(), [str(Path(tmp)/'inputs'),
                str(Path(tmp)/'existing/case_01'), str(Path(tmp)/'output/coefficients/L06/case_01')])

    def test_shared_and_spin_layout_and_legacy_paths(self):
        repo = Path(__file__).resolve().parents[1]
        config = json.loads((repo/'examples/7_spin/production.json').read_text())
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            shared, spins = root/'shared', root/'spins'
            config.update(work=str(shared), spin_work=str(spins), input_root=str(spins))
            config['complete_energy'].update(work=str(spins/'complete'), shared_work=str(shared/'complete'))
            inputs = spins/'inputs/L04'
            inputs.mkdir(parents=True)
            (inputs/'CASE.json').write_text(json.dumps({'angle_configurations': [{}, {}]}))
            path = root/'config.json'
            command = ('import common as c, run_material as r, complete_energy as e; '
                       'print(c.SOURCE); print(c.WORK); print(c.SPINS); '
                       'print(r.case_root(4, 1, 2, .001)[2]); print(e.SHARED); print(e.folder(4, 1))')
            env = dict(os.environ, ANSATZ_RUN_CONFIG=str(path), PYTHONPATH=str(repo/'7_spin'))
            for legacy in (False, True):
                if legacy:
                    del config['spin_work']
                    del config['complete_energy']['shared_work']
                path.write_text(json.dumps(config))
                result = subprocess.run([sys.executable, '-B', '-c', command], cwd=tmp,
                    env=env, capture_output=True, text=True, check=True)
                output = shared if legacy else spins
                energy_shared = spins/'complete' if legacy else shared/'complete'
                self.assertEqual(result.stdout.splitlines(), list(map(str, [repo, shared, output,
                    output/'J/L04/mesh2_h0.001/case_01', energy_shared, spins/'complete/L04/case_01'])))

    def test_complete_functional_and_fixed_electron_response(self):
        complete_readout.energy_selfcheck()
        complete_readout.selfcheck()
        complete_readout.full_functional_selfcheck()
        summarize_energy.selfcheck()

    def test_gaussian_free_energy_fixed_number_gradient(self):
        import common
        self.assertEqual(common.SIGMA, complete_readout.SIGMA)
        self.assertEqual(common.RY, complete_readout.RY)
        spectra = np.array([[-.3, -.02, .1, .4], [-.2, .03, .2, .5]])
        slopes = np.array([[.2, -.4, .1, .3], [-.1, .5, -.3, .2]])
        weights = np.array([.3, .7])
        def read(t):
            e = spectra+t*slopes
            mu = common.chemical_potential(e, weights, 2)
            rows = [complete_readout.fermi(row, mu) for row in e]
            return weights@np.array([r[2] for r in rows]), np.array([r[0] for r in rows])
        h = 1e-5
        derivative = (read(h)[0]-read(-h)[0])/(2*h)
        self.assertAlmostEqual(derivative, weights@np.sum(read(0)[1]*slopes, axis=1), places=8)

    def test_native_energy_columns_units_and_double_counting(self):
        import complete_energy
        # Exercise the production ENERGY.json writer with synthetic QE output.
        rho = np.array([[2., .1, -.2, .3]])
        v = np.array([[.4, -.1, .2, .5]])
        n = np.eye(10, dtype=complex)[None]*.3
        n[0, 0, 1] = .02j
        n[0, 1, 0] = -.02j
        vu = np.eye(10)[None]*.17
        bare = -4.
        band = bare+complete_readout.integral(rho, v, 1.)+complete_readout.trace(n, vu)
        free = band-.07
        qe = np.array([.2, -.3, -.01, 41., 99., 88., -.4, .13])
        ry = complete_readout.RY
        expected, _ = complete_readout.energy(free, 0., rho, v, n, vu,
            qe[0]*ry, qe[1]*ry, qe[2]*ry, qe[6]*ry, 1.)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            native = root/'L04/native'
            native.mkdir(parents=True)
            for i in range(11):
                case = root/f'L04/case_{i:02}'
                case.mkdir()
                np.savez(case/'DENSITY.npz', N=n)
                (case/'ELECTRONIC.json').write_text(json.dumps(dict(F=free, band=band, bare=bare)))
                np.savetxt(native/f'case_{i:02}_ENERGY.dat', qe[None])
            def fake_native(*args, **kwargs):
                kwargs['stdout'].write('GGA_BRANCH_REEXPORT_COMPLETE')
                kwargs['stdout'].flush()
            with patch.object(complete_energy, 'R', root), \
                 patch.object(complete_energy, 'input_file'), \
                 patch.object(complete_energy.subprocess, 'run', side_effect=fake_native):
                complete_energy.native(4)
            result = json.loads((root/'L04/case_00/ENERGY.json').read_text())
            self.assertAlmostEqual(result['energy'], expected, places=12)
            self.assertAlmostEqual(result['terms']['minus_active_potential'], bare-band, places=12)
            self.assertTrue(result['constant_reference_excluded'])

    def test_density_projection_and_fft_gauge(self):
        density_map.selfcheck()
        density_map.fft_selfcheck()

    def test_reference_projection_against_direct_grid(self):
        rng = np.random.default_rng(1620)
        grid = [4, 4, 2]
        p = density_map.Projection.__new__(density_map.Projection)
        p.ix = density_map.indices(grid)
        p.bits = density_map.BITS[[0, 2, 1, 3]]
        p.block = 3
        p.w = rng.normal(size=(8, 12))+1j*rng.normal(size=(8, 12))
        full = np.empty((32, 12), complex)
        for cell, ix in zip(density_map.BITS, p.ix):
            full[ix] = p.w*np.repeat((-1.)**(p.bits@cell), 3)
        fields = rng.normal(size=(4, 32))
        direct = np.array([full.conj().T@(f[:, None]*full) for f in fields])
        np.testing.assert_allclose(p.general(fields), direct, atol=1e-12)

    def test_density_record_roundtrip_and_truncation(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'density.bin'
            values = np.arange(density_map.NG*4, dtype=float).reshape(-1, 4)
            native_input.write_density(path, values, 1)
            np.testing.assert_array_equal(native_input.read_density(path, 1), values)
            path.write_bytes(path.read_bytes()[:-1])
            with self.assertRaises(AssertionError):
                native_input.read_density(path, 1)

    def test_lossless_kz_slice_and_rejection(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp)/'full'
            target = Path(tmp)/'slice'
            prefix = (b'LANIO2_LRVR_V2  '+struct.pack('<8i', 1, 1, 1, 30, 2, 1, 0, 0)
                      +b'a'*64+b'b'*64+bytes(80+16+16))
            size = 16*(1+3)*30*30
            payload = np.arange(size//8, dtype='<f8').tobytes()
            source.write_bytes(prefix+struct.pack('<2d', 0., .125)+bytes(size)+payload)
            with patch.object(kz_production.shutil, 'disk_usage') as disk:
                disk.return_value.free = 100*1024**3
                kz_production.stage_slice(source, target, 1)
            _, sectors, offset, actual_size = kz_production.layout(target)
            self.assertEqual(sectors, (.125,))
            self.assertEqual(actual_size, size)
            self.assertEqual(target.read_bytes()[offset:], payload)
            with self.assertRaises(ValueError):
                kz_production.stage_slice(source, Path(tmp)/'bad', -1)
            target.write_bytes(target.read_bytes()[:-1])
            with self.assertRaises(ValueError):
                kz_production.layout(target)


if __name__ == '__main__':
    unittest.main()
