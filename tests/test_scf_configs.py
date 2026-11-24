"""Input-generation checks; no QE executable or scientific payload is required."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('generate_configs', ROOT/'4_potential/generate_configs.py')
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

    def test_sampling_rules_seed_and_fixed_patterns(self):
        first = configs.configurations(seed=1729)
        self.assertEqual(first, configs.configurations(seed=1729))
        second = configs.configurations(seed=1730)
        expected = (['fm', 'ni56_pp', 'ni56_pm']
                    + [f'cri3_full_{i:02}' for i in range(1, 10)]
                    + [f'seki_upper_{i:02}' for i in range(1, 10)]
                    + ['rel78_p2', 'rel78_m2', 'rel78_p6', 'rel78_m6', 'Neel_anchor'])
        self.assertEqual(list(first), expected)
        for name, spec in first.items():
            angles = np.asarray(spec['angles'])
            self.assertEqual(angles.shape, (4, 2))
            self.assertTrue(np.isfinite(angles).all())
            if name.startswith(('cri3_full_', 'seki_upper_')):
                self.assertNotEqual(spec, second[name])
                self.assertEqual(len(set(map(tuple, angles))), 4)
                self.assertTrue(np.all(angles >= 0))
                self.assertTrue(np.all(angles[:, 1] < 360))
                if name.startswith('cri3_full_'):
                    self.assertTrue(np.all(angles[:, 0] < 180))
                    self.assertTrue(np.any(angles != np.floor(angles)))
                else:
                    self.assertTrue(np.all(angles[:, 0] < 90))
                    np.testing.assert_array_equal(angles, np.floor(angles))
            else:
                self.assertEqual(spec, second[name])
            self.assertEqual(spec['constrained'], name != 'fm')
        np.testing.assert_array_equal(first['ni56_pp']['angles'], [[5, 0], [5, 0], [0, 0], [0, 0]])
        np.testing.assert_array_equal(first['ni56_pm']['angles'], [[5, 0], [5, 180], [0, 0], [0, 0]])
        np.testing.assert_array_equal(first['rel78_m6']['angles'], [[0, 0], [0, 0], [6, 180], [6, 0]])

    def test_relative_pseudopotentials_and_seeded_subset_match_full_generation(self):
        self.scf.write_text(BASE.replace("pseudo_dir='/old/path/with/slash'", "pseudo_dir='pseudo'"))
        manifest = configs.generate(self.scf, self.root/'all', seed=1729)
        self.assertEqual(manifest['seed'], 1729)
        self.assertIn('sampling', manifest)
        self.assertEqual(len(list((self.root/'all').glob('*/scf.in'))), 26)
        names = ['ni56_pm', 'cri3_full_04', 'seki_upper_09']
        subset = configs.generate(self.scf, self.root/'subset', cases=names, seed=1729)
        for name in names:
            self.assertEqual((self.root/f'all/{name}/scf.in').read_text(),
                             (self.root/f'subset/{name}/scf.in').read_text())
            self.assertEqual(manifest['configurations'][name], subset['configurations'][name])
            parsed = configs.read_scf((self.root/f'subset/{name}/scf.in').read_text())
            self.assertEqual(configs.unquote(parsed['groups']['CONTROL']['pseudo_dir']),
                             str(self.pseudo.resolve()))

    def test_invalid_seed_is_rejected_before_writing(self):
        for seed in (-1, 1.5, '2', True):
            with self.subTest(seed=seed):
                with self.assertRaises((TypeError, ValueError)):
                    configs.generate(self.scf, self.root/'cases', self.pseudo, seed=seed)
                self.assertFalse((self.root/'cases').exists())

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
        for cases in ([], ['unknown'], ['fm', 'fm']):
            with self.assertRaises(ValueError):
                configs.generate(self.scf, self.root/'new', self.pseudo, cases)
        self.assertFalse((self.root/'new').exists())


if __name__ == '__main__':
    unittest.main(verbosity=2)
