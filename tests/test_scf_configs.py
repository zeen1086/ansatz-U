"""Input-generation checks; no QE executable or scientific payload is required."""
import importlib.util
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('generate_configs', ROOT/'1_scf/generate_configs.py')
configs = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(configs)

BASE = """&CONTROL
 calculation='scf', prefix='base', pseudo_dir='/old/path/with/slash',
 outdir='./previous/', wfcdir='/shared/wavefunctions', title='a quoted ! comment', max_seconds=42,
/
&SYSTEM
 ibrav=0, nat=4, ntyp=3, ecutwfc=70, ecutrho=280, nbnd=120,
 noncolin=.true., lspinorb=.false., occupations='smearing',
 smearing='gaussian', degauss=0.01, starting_magnetization(2)=0.05,
 nr1=40, nr2=40, nr3=36,
/
&ELECTRONS
 conv_thr=1d-10, mixing_beta=0.03, diagonalization='cg',
 startingpot='file', startingwfc='file',
/
ATOMIC_SPECIES
La 138.90547 La.upf
Ni 58.6934 Ni.upf
O 15.9994 O.upf
CELL_PARAMETERS angstrom
3.9592 0 0
0 3.9592 0
0 0 3.3745
ATOMIC_POSITIONS crystal
La 0 0 0
Ni .5 .5 .5
O 0 .5 .5
O .5 0 .5
K_POINTS automatic
8 4 8 0 0 0
HUBBARD (ortho-atomic)
U Ni-3d 1.0
"""


class ScfConfigurationTests(unittest.TestCase):
    def test_bundled_pseudopotentials_match_frozen_training_records(self):
        expected = {
            'La': '7a5b05be056c96292827043b6c564de17805566e789d043250c57bba1cd10c57',
            'Ni': '07ac84022756676cfc9f0a328f68bd6b19c25876998ca0289ffb682092306946',
            'O': '16c5996fd611f61de83ea416b5248967460ef668837ee05f5b8f3d3dfac79da2',
        }
        for name, digest in expected.items():
            self.assertEqual(hashlib.sha256((ROOT/f'1_scf/pseudo/{name}.upf').read_bytes()).hexdigest(), digest)

    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.pseudo = self.root/'pseudo'
        self.pseudo.mkdir()
        for name in ('La', 'Ni', 'O'):
            (self.pseudo/f'{name}.upf').write_text('test fixture\n')
        self.scf = self.root/'base.in'
        self.scf.write_text(BASE)

    def test_primitive_expansion_preserves_settings_and_site_order(self):
        output = self.root/'cases'
        manifest = configs.generate(self.scf, output, self.pseudo, ['Neel_anchor', 'fm'])
        prepared = configs.read_scf((output/'Neel_anchor/scf.in').read_text())
        system = prepared['groups']['SYSTEM']
        self.assertEqual([system[k] for k in ('nat', 'ntyp', 'nbnd', 'nr1', 'nr2', 'nr3')],
                         ['16', '6', '480', '80', '80', '36'])
        self.assertEqual(system['ecutwfc'], '70')
        self.assertEqual(system['ecutrho'], '280')
        self.assertEqual(prepared['groups']['ELECTRONS']['conv_thr'], '1d-10')
        self.assertEqual(prepared['groups']['ELECTRONS']['mixing_beta'], '0.03')
        self.assertEqual(prepared['groups']['ELECTRONS']['startingpot'], "'atomic'")
        self.assertEqual(prepared['groups']['CONTROL']['title'], "'a quoted ! comment'")
        self.assertEqual(prepared['groups']['CONTROL']['outdir'], "'./tmp'")
        self.assertNotIn('wfcdir', prepared['groups']['CONTROL'])
        self.assertNotIn('/shared/wavefunctions', (output/'Neel_anchor/scf.in').read_text())
        self.assertIn('wfcdir removed', manifest['scratch'])
        self.assertEqual(prepared['cards']['K_POINTS'], ('automatic', ['4 2 8 0 0 0']))
        cell = np.array([[float(x) for x in line.split()] for line in prepared['cards']['CELL_PARAMETERS'][1]])
        np.testing.assert_allclose(cell, np.diag([7.9184, 7.9184, 3.3745]), atol=1e-12)
        ni = [line.split() for line in prepared['cards']['ATOMIC_POSITIONS'][1] if line.startswith('Ni')]
        self.assertEqual([row[0] for row in ni], ['Ni1', 'Ni2', 'Ni3', 'Ni4'])
        np.testing.assert_array_equal(np.array([row[1:] for row in ni], float), configs.NI_POSITIONS)
        self.assertEqual(prepared['cards']['HUBBARD'][1], [f'U Ni{i}-3d 1.0' for i in range(1, 5)])
        np.testing.assert_allclose(manifest['configurations']['Neel_anchor']['target_directions'],
                                   [[0, 0, 1], [0, 0, -1], [0, 0, -1], [0, 0, 1]], atol=2e-16)
        self.assertNotIn('constrained_magnetization', configs.read_scf((output/'fm/scf.in').read_text())['groups']['SYSTEM'])
        self.assertEqual(json.loads((output/'CONFIGURATIONS.json').read_text()), manifest)

    def test_four_ni_input_and_explicit_kpoints_are_not_expanded_again(self):
        base = configs.prepare(configs.read_scf(BASE), self.pseudo)
        text, _ = configs.render(base, [[0, 0]] * 4)
        text = text.replace('K_POINTS automatic\n4 2 8 0 0 0', 'K_POINTS crystal\n2\n0 0 0 0.5\n0 0 .5 .5')
        self.scf.write_text(text)
        configs.generate(self.scf, self.root/'cases', self.pseudo, ['Neel_anchor'])
        result = configs.read_scf((self.root/'cases/Neel_anchor/scf.in').read_text())
        self.assertEqual(result['groups']['SYSTEM']['nbnd'], '480')
        self.assertEqual(result['cards']['CELL_PARAMETERS'], configs.read_scf(text)['cards']['CELL_PARAMETERS'])
        self.assertEqual(result['cards']['K_POINTS'], ('crystal', ['2', '0 0 0 0.5', '0 0 .5 .5']))

    def test_tetragonal_cell_and_default_v1_hubbard(self):
        text = BASE.replace('ibrav=0', 'ibrav=6, celldm(1)=7.48180367, celldm(3)=0.85231865')
        text = text.replace('CELL_PARAMETERS angstrom\n3.9592 0 0\n0 3.9592 0\n0 0 3.3745\n', '')
        text = text.replace('HUBBARD (ortho-atomic)\nU Ni-3d 1.0\n', '')
        result = configs.prepare(configs.read_scf(text), self.pseudo)
        self.assertAlmostEqual(result['cell'][0, 0], 2*7.48180367*configs.BOHR_ANGSTROM)
        self.assertEqual(result['hubbard'], ['1.0'] * 4)
        self.assertIn('added', result['hubbard_source'])

    def test_invalid_inputs_fail_before_creating_output(self):
        for text in (BASE.replace('8 4 8 0 0 0', '7 4 8 0 0 0'),
                     BASE.replace('8 4 8 0 0 0', '8 4 8 1 0 0'),
                     BASE.replace('U Ni-3d 1.0', 'U Ni-3d 2.0'),
                     BASE.replace('HUBBARD (ortho-atomic)', 'HUBBARD (atomic)'),
                     BASE.replace('lspinorb=.false.', 'lspinorb=.true.'),
                     BASE.replace('nat=4,', 'nat=4, tot_magnetization=4,'),
                     BASE.replace("pseudo_dir='/old/path/with/slash'", "pseudo_dir='bad/'/../PP/'"),
                     BASE + '\nATOMIC_FORCES\nLa 0 0 0\n',
                     BASE.replace('nat=4,', 'nat=4, nat=4,'),
                     BASE.replace('nat=4,', 'nat=4, starting_charge(2)=0.3,')):
            with self.subTest(text=text[:50]):
                self.scf.write_text(text)
                with self.assertRaises(ValueError):
                    configs.generate(self.scf, self.root/'cases', self.pseudo, ['fm'])
                self.assertFalse((self.root/'cases').exists())
        self.scf.write_text(BASE)
        (self.pseudo/'Ni.upf').unlink()
        with self.assertRaises(FileNotFoundError):
            configs.generate(self.scf, self.root/'cases', self.pseudo, ['fm'])

    def test_existing_output_and_invalid_case_selection_are_rejected(self):
        output = self.root/'cases'
        output.mkdir()
        (output/'keep').write_text('user data')
        with self.assertRaises(FileExistsError):
            configs.generate(self.scf, output, self.pseudo, ['fm'])
        self.assertEqual((output/'keep').read_text(), 'user data')
        for cases in ([], ['unknown'], ['cri3_full_01'], ['ni56_pp'], ['fm', 'fm']):
            with self.assertRaises(ValueError):
                configs.generate(self.scf, self.root/'new', self.pseudo, cases)
        self.assertFalse((self.root/'new').exists())


if __name__ == '__main__':
    unittest.main(verbosity=2)
