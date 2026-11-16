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

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'7_spin/production'))
import complete_readout
import density_map
import energy_torque
import kz_production
import native_input
import reference_projection
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
                       PYTHONPATH=str(repo/'7_spin/production'))
            result = subprocess.run([sys.executable, '-B', '-c', command], cwd=tmp,
                                    env=env, capture_output=True, text=True, check=True)
            self.assertEqual(result.stdout.splitlines(), [str(Path(tmp)/'inputs'),
                str(Path(tmp)/'existing/case_01'), str(Path(tmp)/'output/coefficients/L06/case_01')])

    def test_complete_functional_and_fixed_electron_response(self):
        energy_torque.selfcheck()
        complete_readout.selfcheck()
        complete_readout.full_functional_selfcheck()
        summarize_energy.selfcheck()

    def test_density_projection_and_fft_gauge(self):
        density_map.selfcheck()
        density_map.fft_selfcheck()

    def test_reference_projection_against_direct_grid(self):
        rng = np.random.default_rng(1620)
        grid = [4, 4, 2]
        p = reference_projection.Projection.__new__(reference_projection.Projection)
        p.ix = reference_projection.indices(grid)
        p.bits = reference_projection.BITS[[0, 2, 1, 3]]
        p.block = 3
        p.w = rng.normal(size=(8, 12))+1j*rng.normal(size=(8, 12))
        full = np.empty((32, 12), complex)
        for cell, ix in zip(reference_projection.BITS, p.ix):
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
