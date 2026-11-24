"""Synthetic SCF label-contract checks; no QE executable or physical data needed."""
import json
from pathlib import Path
import struct
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

import qe_labels
import extract_labels


def record(stream, values):
    payload = values.tobytes(order='F')
    marker = struct.pack('<i', len(payload))
    stream.write(marker + payload + marker)


class LabelTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.save = self.root/'fixture.save'
        self.save.mkdir()
        self.cell = np.diag([8., 8., 4.])
        self.grid = np.array([4, 4, 2], dtype=np.int32)
        self.names = ['La']*4 + ['Ni1', 'Ni2', 'Ni3', 'Ni4'] + ['O']*8
        self.expected_N = np.array([np.eye(10, dtype=complex)*(.3 + .1*i) for i in range(4)])
        for i, n in enumerate(self.expected_N):
            n[0, 3] = .02 + .01j*(i+1)
            n[3, 0] = n[0, 3].conjugate()
            n[2, 8] = -.03 + .02j
            n[8, 2] = n[2, 8].conjugate()
        self.write_save()
        self.log = self.root/'scf.out'
        self.log.write_text('Program PWSCF v.7.6 starts\n'
                            'number of atoms/cell = 16\n'
                            'Fixed quantization axis for GGA: 0.000000 0.000000 1.000000\n'
                            '! total energy = -20.0000000000 Ry\n'
                            'convergence has been achieved in 9 iterations\n'
                            'JOB DONE.\n')
        self.fields = np.column_stack([np.arange(32.) + 7.,
                                      np.arange(32.)*.1 + 2.,
                                      np.arange(32.)*-.2 - 3.,
                                      np.arange(32.)*.3 + 4.])
        self.moments = np.zeros((16, 5))
        self.moments[:, 0] = np.arange(1, 17)
        self.moments[:, 4] = 2.
        self.moments[4:8, 1:4] = [[.1, .2, .7], [-.2, .1, .6], [.3, -.1, .5], [.2, .3, -.6]]

    def write_save(self):
        species = ''.join(f'<species name="{name}"><starting_magnetization>{.05 if name.startswith("Ni") else 0}</starting_magnetization>'
                          '<spin_teta>0</spin_teta><spin_phi>0</spin_phi></species>'
                          for name in dict.fromkeys(self.names))
        ni_positions = {'Ni1': '2 2 2', 'Ni2': '6 2 2', 'Ni3': '2 6 2', 'Ni4': '6 6 2'}
        atoms = ''.join(f'<atom index="{i+1}" name="{name}">{ni_positions.get(name, f"{i%4} {(i//4)%4} {i%2}")}</atom>'
                        for i, name in enumerate(self.names))
        hubbard = ''.join(f'<Hubbard_U specie="{name}" label="3d">{1/qe_labels.HARTREE_EV:.17g}</Hubbard_U>'
                          for name in self.names if name.startswith('Ni'))
        for n, atom in zip(self.expected_N, [i for i, name in enumerate(self.names) if name.startswith('Ni')]):
            magnitudes = ' '.join(f'{v:.17g}' for v in np.abs(n.T).ravel(order='F'))
            hubbard += (f'<Hubbard_ns_nc index="{atom+1}" specie="{self.names[atom]}" label="3d" '
                        f'rank="2" dims="10 10" order="F" spin="1">{magnitudes}</Hubbard_ns_nc>')
        xml = f'''<espresso Units="Hartree atomic units">
<input><control_variables><calculation>scf</calculation><prefix>fixture</prefix></control_variables>
<spin><noncolin>true</noncolin><spinorbit>false</spinorbit><lsda>false</lsda></spin>
<dft><functional>PBE</functional></dft>
<atomic_species>{species}</atomic_species><atomic_structure><atomic_positions>{atoms}</atomic_positions></atomic_structure></input>
<output><convergence_info><scf_conv><convergence_achieved>true</convergence_achieved><n_scf_steps>9</n_scf_steps><scf_error>1e-12</scf_error></scf_conv></convergence_info>
<atomic_structure nat="16"><cell><a1>8 0 0</a1><a2>0 8 0</a2><a3>0 0 4</a3></cell><atomic_positions>{atoms}</atomic_positions></atomic_structure>
<dft><functional>PBE</functional><dftU><lda_plus_u_kind>0</lda_plus_u_kind><U_projection_type>ortho-atomic</U_projection_type>{hubbard}</dftU></dft>
<basis_set><fft_grid nr1="4" nr2="4" nr3="2"/><fft_smooth nr1="2" nr2="2" nr3="2"/></basis_set>
<total_energy><etot>-10</etot></total_energy></output></espresso>'''
        (self.save/'data-file-schema.xml').write_text(xml)
        raw = np.zeros((5, 5, 4, 16), dtype=complex)
        for n, atom in zip(self.expected_N, [i for i, name in enumerate(self.names) if name.startswith('Ni')]):
            for a in range(10):
                for b in range(10):
                    raw[a%5, b%5, 2*(a//5)+b//5, atom] = n[b, a]
        (self.save/'occup.txt').write_text('\n'.join(f'({z.real:.17g},{z.imag:.17g})' for z in raw.ravel(order='F')) + '\n')
        (self.save/'charge-density.dat').write_bytes(b'synthetic density identity')

    def write_fields(self, directory):
        directory.mkdir(exist_ok=True)
        with (directory/'V_LABELS.bin').open('wb') as stream:
            record(stream, self.grid.astype('<i4'))
            record(stream, self.fields.astype('<f8'))
        (directory/'LABELS.meta').write_text('V1_LABEL_FIELDS_RY\n4 4 2\n1 0 0 1\n8 0 0 0 8 0 0 0 4\n')
        np.savetxt(directory/'moments.dat', self.moments, fmt='%.17g')

    def test_saved_complex_occupations_and_arbitrary_atom_ids(self):
        self.names = ['Ni1', 'La', 'O', 'Ni2', 'O', 'La', 'Ni3', 'O', 'La', 'O', 'Ni4', 'La', 'O', 'O', 'O', 'O']
        self.write_save()
        saved = qe_labels.read_save(self.save)
        np.testing.assert_array_equal(saved['ni'], [0, 3, 6, 10])
        np.testing.assert_allclose(saved['N'], self.expected_N, atol=1e-16, rtol=0)
        self.assertEqual(saved['N'].dtype, np.complex128)
        np.testing.assert_allclose(saved['cell'], self.cell*qe_labels.BOHR_ANGSTROM)
        np.testing.assert_array_equal(saved['grid_dense'], self.grid)
        self.assertEqual(saved['metadata']['orbital_order'], ['dz2', 'dxz', 'dyz', 'dx2-y2', 'dxy'])

    def test_reject_incomplete_unconverged_or_invalid_save(self):
        path = self.save/'data-file-schema.xml'
        original = path.read_text()
        bad = [original.replace('<convergence_achieved>true', '<convergence_achieved>false'),
               original.replace('index="16"', 'index="15"'),
               original.replace('<spinorbit>false', '<spinorbit>true'),
               original.replace('<a1>8 0 0', '<a1>NaN 0 0'),
               original.replace('name="Ni1">2 2 2', 'name="Ni1">3 2 2')]
        for text in bad:
            with self.subTest(text=text[-80:]):
                path.write_text(text)
                with self.assertRaises(ValueError):
                    qe_labels.read_save(self.save)
        path.write_text(original)
        (self.save/'occup.txt').write_text('(0,0)\n')
        with self.assertRaises(ValueError):
            qe_labels.read_save(self.save)
        self.expected_N[0, 0, 0] = 1.2
        self.write_save()
        with self.assertRaises(ValueError):
            qe_labels.read_save(self.save)

    def test_reject_non_pbe_functional(self):
        path = self.save/'data-file-schema.xml'
        path.write_text(path.read_text().replace('<functional>PBE', '<functional>PZ'))
        with self.assertRaisesRegex(ValueError, 'PBE functional'):
            qe_labels.read_save(self.save)

    def test_reject_occupations_from_a_different_saved_state(self):
        qe_labels.read_save(self.save)
        path = self.save/'data-file-schema.xml'
        original = path.read_text()
        self.expected_N *= .9
        self.write_save()
        path.write_text(original)
        with self.assertRaisesRegex(ValueError, 'XML Hubbard matrix magnitudes'):
            qe_labels.read_save(self.save)

    def test_scf_log_matches_saved_state_and_gga_axis(self):
        saved = qe_labels.read_save(self.save)
        branch = extract_labels.read_scf_log(self.log, saved)
        self.assertTrue(branch['lsign'])
        np.testing.assert_allclose(branch['ux'], [0, 0, 1])
        original = self.log.read_text()
        for text in (original.replace('JOB DONE.', ''),
                     original.replace('convergence has been achieved', 'convergence NOT achieved'),
                     original.replace('-20.0000000000', '-21.0000000000'),
                     original.replace('atoms/cell = 16', 'atoms/cell = 4'),
                     original.replace('0.000000 0.000000 1.000000', '1.000000 0.000000 0.000000')):
            self.log.write_text(text)
            with self.assertRaises(ValueError):
                extract_labels.read_scf_log(self.log, saved)
        self.log.write_text(original.replace('Fixed quantization axis for GGA: 0.000000 0.000000 1.000000\n', ''))
        branch = extract_labels.read_scf_log(self.log, saved)
        self.assertFalse(branch['lsign'])

    def test_field_grid_channel_order_gauge_and_actual_moments(self):
        saved = qe_labels.read_save(self.save)
        branch = extract_labels.read_scf_log(self.log, saved)
        raw = self.root/'raw'
        self.write_fields(raw)
        fields, moments = extract_labels.read_fields(raw, saved, branch)
        np.testing.assert_array_equal(fields, self.fields)
        stage = self.root/'package'
        stage.mkdir()
        extract_labels.write_package(stage, saved, branch, fields, moments, {'source': 'synthetic'})
        np.testing.assert_array_equal(np.load(stage/'scalar.npy'), self.fields[:, 0]-self.fields[:, 0].mean())
        np.testing.assert_array_equal(np.load(stage/'B.npy'), self.fields[:, 1:])
        np.testing.assert_allclose(np.load(stage/'N.npy'), self.expected_N)
        np.testing.assert_array_equal(np.load(stage/'moments.npy'), self.moments[4:8, 1:4])
        np.testing.assert_array_equal(np.load(stage/'atom_ids.npy'), [5, 6, 7, 8])
        np.testing.assert_array_equal(np.load(stage/'grid.npy'), self.grid)
        with np.load(stage/'structure.npz') as geometry:
            np.testing.assert_allclose(geometry['spins'], self.moments[4:8, 1:4]/np.linalg.norm(self.moments[4:8, 1:4], axis=1)[:, None])
        self.assertTrue(json.loads((stage/'LABELS.json').read_text()))

    def test_reject_mixed_gga_grid_geometry_atoms_and_nonfinite_fields(self):
        saved = qe_labels.read_save(self.save)
        branch = extract_labels.read_scf_log(self.log, saved)
        raw = self.root/'raw'
        self.write_fields(raw)
        meta = (raw/'LABELS.meta').read_text()
        for text in (meta.replace('1 0 0 1', '0 0 0 0'),
                     meta.replace('4 4 2', '4 4 3'),
                     meta.replace('8 0 0 0 8', '9 0 0 0 8')):
            (raw/'LABELS.meta').write_text(text)
            with self.assertRaises(ValueError):
                extract_labels.read_fields(raw, saved, branch)
        self.write_fields(raw)
        self.moments[15, 0] = 15
        np.savetxt(raw/'moments.dat', self.moments)
        with self.assertRaises(ValueError):
            extract_labels.read_fields(raw, saved, branch)
        self.moments[15, 0] = 16
        self.fields[0, 2] = np.nan
        self.write_fields(raw)
        with self.assertRaises(ValueError):
            extract_labels.read_fields(raw, saved, branch)

    def adapter(self):
        executable = self.root/'export_labels.x'
        executable.write_text('#!/bin/sh\nexit 0\n')
        executable.chmod(0o700)
        return executable

    def test_extract_publishes_complete_package_and_preserves_sources(self):
        executable = self.adapter()
        original = {p: p.read_bytes() for p in self.save.iterdir()}
        def run_adapter(command, **kwargs):
            self.assertEqual(command, [str(executable.resolve()), '-in', 'extract.in'])
            request = (kwargs['cwd']/'extract.in').read_text()
            self.assertIn('source_lsign=.true.', request)
            self.assertIn('source_ux=0, 0, 1', request)
            self.write_fields(kwargs['cwd'])
            kwargs['stdout'].write('V1_LABEL_EXPORT_DONE\n')
        output = self.root/'published'
        with patch.object(extract_labels.subprocess, 'run', side_effect=run_adapter) as run:
            result = extract_labels.extract(self.save, self.log, executable, output)
            self.assertEqual(run.call_count, 1)
        self.assertEqual(result['units']['scalar'], 'Ry')
        np.testing.assert_array_equal(np.load(output/'B.npy'), self.fields[:, 1:])
        self.assertTrue((output/'LABELS.json').is_file())
        self.assertFalse((output/'V_LABELS.bin').exists())
        self.assertEqual(original, {p: p.read_bytes() for p in self.save.iterdir()})
        with patch.object(extract_labels.subprocess, 'run') as run:
            with self.assertRaises(FileExistsError):
                extract_labels.extract(self.save, self.log, executable, output)
            run.assert_not_called()
        self.assertEqual(result, json.loads((output/'LABELS.json').read_text()))

    def test_adapter_failure_or_changed_source_never_publishes(self):
        executable = self.adapter()
        output = self.root/'failed'
        for failure in (subprocess.CalledProcessError(1, [str(executable)]),
                        subprocess.TimeoutExpired([str(executable)], 900)):
            def failed_adapter(command, **kwargs):
                kwargs['stdout'].write(''.join(f'adapter line {i}\n' for i in range(25)))
                kwargs['stdout'].write('synthetic failure marker\n')
                raise failure
            with patch.object(extract_labels.subprocess, 'run', side_effect=failed_adapter):
                with self.assertRaisesRegex(RuntimeError, 'synthetic failure marker') as caught:
                    extract_labels.extract(self.save, self.log, executable, output)
            self.assertNotIn('adapter line 0\n', str(caught.exception))
            self.assertFalse(output.exists())
            self.assertFalse(list(self.root.glob('.failed-*')))
        def changed_source(command, **kwargs):
            self.write_fields(kwargs['cwd'])
            kwargs['stdout'].write('V1_LABEL_EXPORT_DONE\n')
            (self.save/'charge-density.dat').write_bytes(b'changed while extracting')
        with patch.object(extract_labels.subprocess, 'run', side_effect=changed_source):
            with self.assertRaisesRegex(ValueError, 'changed during extraction'):
                extract_labels.extract(self.save, self.log, executable, output)
        self.assertFalse(output.exists())
        self.assertFalse(list(self.root.glob('.failed-*')))


if __name__ == '__main__':
    unittest.main()
