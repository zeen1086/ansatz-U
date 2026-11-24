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


class ProductionEnergyTest(unittest.TestCase):
    def test_merged_commands_import_without_calculation_and_dispatch_selfchecks(self):
        repo = Path(__file__).resolve().parents[1]
        config = json.loads((repo/'examples/7_spin/production.json').read_text())
        config.pop('complete_energy')
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'config.json'
            path.write_text(json.dumps(config))
            env = dict(os.environ, ANSATZ_RUN_CONFIG=str(path),
                       PYTHONPATH=os.pathsep.join((str(repo), str(repo/'7_spin'))))
            script = '''from pathlib import Path
from unittest.mock import patch
with patch.object(Path, 'mkdir', side_effect=AssertionError('created a directory')), \\
     patch('subprocess.run', side_effect=AssertionError('launched a calculation')):
    import native_input, complete_readout
'''
            subprocess.run([sys.executable, '-B', '-c', script], cwd=tmp, env=env,
                           check=True, capture_output=True, text=True)
            for name, command in (('native_input.py', 'check'), ('complete_readout.py', 'summary')):
                result = subprocess.run([sys.executable, '-B', str(repo/'7_spin'/name), '--help'],
                                        cwd=tmp, env=env, check=True, capture_output=True, text=True)
                self.assertIn(command, result.stdout)
            subprocess.run([sys.executable, '-B', str(repo/'7_spin/complete_readout.py'), 'selfcheck'],
                           cwd=tmp, env=env, check=True, capture_output=True, text=True)
            self.assertEqual(list(Path(tmp).iterdir()), [path])

    def test_merged_sources_bind_cache_and_reject_changed_identity(self):
        import complete_energy as driver
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            reference, wannier = root/'reference', root/'wannier'
            for folder in (reference/'results', reference/'inputs/reference_actual', wannier):
                folder.mkdir(parents=True)
            for name in ('basis/bloch.npy', 'basis/pure_H0.au', 'RVR/PROJECTION.json',
                         'RVR/B_general', 'model/V_func_tmp/V_FUNC.json'):
                path = root/name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(b'identity fixture')
            with patch.multiple(driver, WORK=root, INTERFACE=root/'model', REF=wannier,
                                MIN=reference, SHARED=root):
                identity = driver.cache_identity()
                self.assertIn('code/complete_readout.py', identity['files'])
                self.assertIn('code/native_input.py', identity['files'])
                source_names = {name[5:] for name in identity['files'] if name.startswith('code/')}
                spin_sources = {path.name for path in (Path(driver.__file__).parent).glob('*.py')}
                self.assertTrue(source_names-{'spin_coordinates.py', 'frozen_io.f90', 'cal_Hsk.f90'} <= spin_sources)
                driver.bind_cache(identity)
                driver.bind_cache(identity)
                changed = json.loads(json.dumps(identity))
                changed['files']['code/complete_readout.py'] = 'changed source identity'
                with self.assertRaisesRegex(ValueError, 'Shared inputs changed'):
                    driver.bind_cache(changed)

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
                       'print(c.ROOT); print(c.coefficient_case(4, 1)); print(c.coefficient_case(6, 1))')
            env = dict(os.environ, ANSATZ_RUN_CONFIG=str(path),
                       PYTHONPATH=str(repo/'7_spin'))
            result = subprocess.run([sys.executable, '-B', '-c', command], cwd=tmp,
                                    env=env, capture_output=True, text=True, check=True)
            self.assertEqual(result.stdout.splitlines(), [str(Path(tmp)/'inputs'),
                str(Path(tmp)/'existing/case_01'), str(Path(tmp)/'output/coefficients/L06/case_01')])

    def test_shared_and_spin_layout_requires_explicit_paths(self):
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
            path.write_text(json.dumps(config))
            result = subprocess.run([sys.executable, '-B', '-c', command], cwd=tmp,
                env=env, capture_output=True, text=True, check=True)
            self.assertEqual(result.stdout.splitlines(), list(map(str, [repo, shared, spins,
                spins/'J/L04/mesh2_h0.001/case_01', shared/'complete', spins/'complete/L04/case_01'])))
            for section, key in ((config, 'spin_work'), (config['complete_energy'], 'shared_work')):
                value = section.pop(key)
                path.write_text(json.dumps(config))
                result = subprocess.run([sys.executable, '-B', '-c', command], cwd=tmp,
                    env=env, capture_output=True, text=True)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(key, result.stderr)
                section[key] = value

    def test_complete_functional_and_fixed_electron_response(self):
        complete_readout.energy_selfcheck()
        complete_readout.selfcheck()
        complete_readout.full_functional_selfcheck()
        complete_readout.curvature_selfcheck()
        n = np.eye(10, dtype=complex)*.4
        n[1, 5], n[5, 1] = .04+.07j, .04-.07j
        direction = np.zeros_like(n)
        direction[1, 5], direction[5, 1] = .2-.3j, .2+.3j
        ueff, step = 2.75, 1e-5
        _, potential = complete_readout.hubbard(n, ueff)
        np.testing.assert_allclose(potential, ueff*(.5*np.eye(10)-n))
        eplus, eminus = (complete_readout.hubbard(n+s*step*direction, ueff)[0] for s in (1, -1))
        self.assertAlmostEqual((eplus-eminus)/(2*step), np.trace(potential@direction).real, places=9)

    def test_complete_spectrum_is_invariant_under_a_scalar_gauge_shift(self):
        import complete_energy
        from common import chemical_potential

        length, k = 2, np.array([.5, 0., .125])
        ueff = 2.75
        cells = length**2
        dimension = 60*cells
        eye = np.eye(dimension)
        wave = np.diag(np.linspace(.97, 1.03, 30)).astype(complex)
        wave[0, 1] = .03j
        psi = np.array([(1.+.01*i)*wave for i in range(cells)])
        overlap, transform, gram_max = complete_readout.metric(psi)
        base = np.diag(np.linspace(-1., 1., dimension))

        def assemble(coeff, _, directories, *args, **kwargs):
            # cal_Hsk adds the coefficient offset to the native diagonal.
            offset = struct.unpack_from('<d', coeff.read_bytes(), 120)[0]
            directories[0].mkdir()
            np.asarray(base+offset*eye, dtype='<c16').ravel(order='F').tofile(
                directories[0]/'H.k000001.h')

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            inputs = root/'inputs'/f'L{length:02}'
            inputs.mkdir(parents=True)
            case = {'readout_sites_native_1based': [1]}
            (inputs/'CASE.json').write_text(json.dumps({'angle_configurations': [case, case]}))
            np.save(root/'transform.npy', transform)
            np.save(root/'overlap.npy', overlap)
            (root/'metric.json').write_text(json.dumps({'gram_max': gram_max}))
            observations = []
            with patch.multiple(complete_energy, ROOT=root, R=root/'output', ENERGY={},
                                points=lambda *_: [(1, k)], native_batch=assemble,
                                barrier=lambda *_: None):
                for case_index, offset in enumerate((-.7, .3)):
                    coeff = complete_energy.coefficient_case(length, case_index)/'coeff'/'value.bin'
                    coeff.parent.mkdir(parents=True)
                    coeff.write_bytes(struct.pack('<16s8i64s2d', b'LANIO2_COEFF_V2 ',
                        1, 1, length, cells, 0, 0, 0, 0, complete_energy.CONFIG['model_id'].encode(),
                        ueff, offset) + bytes(cells*(4*8 + 100*16)))
                    complete_energy.spectrum(length, case_index, 1, {1: root/'wave.npy'}, {1: root}, root/'rvr', ueff)
                    with np.load(complete_energy.folder(length, case_index)/'k_0001.npz') as data:
                        e, Q, H = data['e'], data['Q'], data['H']
                    mu = chemical_potential(e[None], np.ones(1), 21*cells)
                    f = complete_readout.fermi(e, mu)[0]
                    rho = density_map.fold_density(psi, (Q*f)@Q.conj().T, length)
                    observations.append((H, e, mu, f, rho))
                # An existing spectrum cannot be reused with a mismatched coefficient header.
                with self.assertRaisesRegex(ValueError, 'Coefficient Ueff'):
                    complete_energy.spectrum(length, 1, 1, {1: root/'wave.npy'}, {1: root}, root/'rvr', 1.)
        before, after = observations
        for name, a, b, shift in zip(('H', 'e', 'mu', 'occupations', 'density'),
                                     before, after, (overlap, 1., 1., 0., 0.)):
            with self.subTest(quantity=name):
                np.testing.assert_allclose(b-a, shift, atol=1e-12, rtol=0)

    def test_complete_energy_rejects_unbound_result_with_unknown_gauge(self):
        import complete_energy
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            case = root/'L06/case_00'
            case.mkdir(parents=True)
            np.savez(case/'k_0000.npz', H=np.eye(1))
            with patch.multiple(complete_energy, R=root, SIZES=[6], compute_node=lambda: None):
                with self.assertRaisesRegex(ValueError, 'unknown energy gauge'):
                    complete_energy.main(0)

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
        n[0, 0, 6] = .02j
        n[0, 6, 0] = -.02j
        vu = np.eye(10)[None]*.17
        bare = -4.
        band = bare+complete_readout.integral(rho, v, 1.)+complete_readout.trace(n, vu)
        free = band-.07
        qe = np.array([.2, -.3, -.01, 41., 99., 88., -.4, .13])
        ry = complete_readout.RY
        ueff = 2.75
        expected, _ = complete_readout.energy(free, 0., rho, v, n, vu,
            qe[0]*ry, qe[1]*ry, qe[2]*ry, qe[6]*ry, 1., ueff)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            native = root/'L04/native'
            native.mkdir(parents=True)
            binary=root/'qe';binary.write_bytes(b'fixture')
            for i in range(11):
                case = root/f'L04/case_{i:02}'
                case.mkdir()
                np.savez(case/'DENSITY.npz', N=n)
                (case/'DENSITY.bin').write_bytes(b'density')
                (case/'ELECTRONIC.json').write_text(json.dumps(dict(F=free, band=band, bare=bare)))
                np.savetxt(native/f'case_{i:02}_ENERGY.dat', qe[None])
            def fake_native(*args, **kwargs):
                for i in range(11):
                    np.savetxt(kwargs['cwd']/f'case_{i:02}_ENERGY.dat',qe[None])
                    (kwargs['cwd']/f'case_{i:02}_TRUE.bin').write_bytes(b'field')
                kwargs['stdout'].write('GGA_BRANCH_REEXPORT_COMPLETE')
                kwargs['stdout'].flush()
            with patch.object(complete_energy, 'R', root), \
                 patch.dict(complete_energy.ENERGY, {'native_binary':str(binary)}), \
                 patch.object(complete_energy, 'read_density', return_value=np.zeros((1,4))), \
                 patch.object(complete_energy, 'input_file'), \
                 patch.object(complete_energy.subprocess, 'run', side_effect=fake_native):
                complete_energy.native(4, ueff)
            result = json.loads((root/'L04/case_00/ENERGY.json').read_text())
            self.assertAlmostEqual(result['energy'], expected, places=12)
            self.assertEqual(result['Ueff_eV'], ueff)
            self.assertAlmostEqual(result['terms']['EU'], .5*ueff*(10*(.3-.3**2)-2*.02**2), places=12)
            self.assertAlmostEqual(result['terms']['minus_active_potential'], bare-band, places=12)
            self.assertTrue(result['constant_reference_excluded'])

    def test_global_spin_frame_in_case_and_density_assembly(self):
        import complete_energy as driver
        import run_material

        # A 120-degree rotation about (1,1,1) cycles x -> y -> z -> x.
        local = dict(perturb_site_native_1based=1, angle_deg=.25)
        rotated = dict(local, global_spin_rotation_deg=[120/np.sqrt(3)]*3)
        spins = np.array([[0.,0.,1.],[0.,0.,-1.],[0.,0.,-1.],[0.,0.,1.]])
        moments = np.array([[.01,-.02,.03],[-.02,.03,.01],[.04,.01,-.03],[.03,-.01,.02]])
        orbital = np.eye(5, dtype=complex)
        orbital[1,2], orbital[2,1] = .2j, -.2j
        def shells(vectors):
            return np.array([np.kron([[z, x-1j*y],[x+1j*y, -z]], orbital) for x,y,z in vectors])

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            inputs = root/'inputs/L02'
            inputs.mkdir(parents=True)
            (inputs/'CASE.json').write_text(json.dumps({'angle_configurations':[local, rotated]}))
            np.savez(inputs/'structure.npz', spins=spins)
            ref = root/'reference'
            ref.mkdir()
            density = np.column_stack((np.full(4,41.),moments))[:,None,:]
            np.savez(ref/'background.npz', density=density, N=shells(moments))
            with patch.object(run_material,'ROOT',root), patch.multiple(
                    driver, ROOT=root, R=root, SHARED=root, NG=1, cells_to_grid=lambda a,L:a[:,0]):
                original_spins=run_material.spin_structure(2,local)['spins']
                self.assertGreater(original_spins[0,0], 0.)
                np.testing.assert_allclose(run_material.spin_structure(2,rotated)['spins'],
                                           original_spins[:,[2,0,1]], atol=1e-15)
                for i, expected in enumerate((moments,moments[:,[2,0,1]])):
                    case=driver.folder(2,i)
                    case.mkdir(parents=True)
                    active=.25*np.broadcast_to(np.eye(10),(4,10,10))
                    for z in range(8):
                        np.savez(case/f'density_z{z}.npz', density=np.zeros_like(density),
                                 N=active/8, bare=0., F=0., band=0., count=0.)
                    driver.assemble_density(2,i,1.)
                    with np.load(case/'DENSITY.npz') as output:
                        np.testing.assert_allclose(output['N'], active+shells(expected), atol=1e-15)
                    field=native_input.read_density(case/'DENSITY.bin',2)
                    np.testing.assert_allclose(field[:,1:], expected, atol=1e-15)
                    np.testing.assert_array_equal(field[:,0], density[:,0,0])
                for invalid in ([0,180], [0,float('nan'),0]):
                    with self.assertRaisesRegex(ValueError,'rotation vector'):
                        run_material.spin_structure(2,dict(local,global_spin_rotation_deg=invalid))

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
