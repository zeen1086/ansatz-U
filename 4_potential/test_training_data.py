"""Synthetic label-loading checks for the single-response training route."""
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from extract_labels import sha256, write_package
from training_data import RY_EV, initial_reference, load_dataset
from model import Reference


def label_fixture(path: Path, spins: np.ndarray | None = None) -> dict:
    path.mkdir()
    if spins is None:
        spins = np.array([[0., 0., 1.], [1., 0., 0.], [0., 1., 0.], [0., 0., -1.]])
    cell = np.diag([8., 8., 4.])
    positions = np.zeros((16, 3))
    positions[4:8] = [[2., 2., 2.], [6., 2., 2.], [2., 6., 2.], [6., 6., 2.]]
    species = np.array([57]*4 + [28]*4 + [8]*8)
    ns, nv = np.zeros(15), np.zeros(15)
    ns[:5] = [.7, .8, .8, .75, .65]
    nv[:5] = [.1, .05, .05, .07, .02]
    reference = Reference(ns, nv, np.zeros(2), np.zeros(2), 'synthetic')
    n = reference.values(spins, np.ones((2, 4)))[0]
    saved = dict(cell=cell, positions=positions, species=species, ni=np.arange(4, 8),
                 N=n, grid_dense=np.array([4, 4, 2]),
                 metadata=dict(functional='PBE', U_eV=1., projector='ortho-atomic',
                               orbital_order=['dz2', 'dxz', 'dyz', 'dx2-y2', 'dxy'],
                               spin_order='up five orbitals, down five orbitals'))
    fields = np.arange(128, dtype=float).reshape(32, 4)/100
    moments = spins*np.array([.7, .8, .9, 1.])[:, None]
    manifest = write_package(path, saved, dict(lsign=False, ux=[0., 0., 0.]), fields, moments, {})
    return dict(saved=saved, fields=fields, manifest=manifest, reference=reference)


def rewrite_array(path: Path, name: str, value: np.ndarray) -> None:
    np.save(path/name, value, allow_pickle=False)
    manifest = json.loads((path/'LABELS.json').read_text())
    manifest['files'][name] = sha256(path/name)
    (path/'LABELS.json').write_text(json.dumps(manifest))


class TrainingDataTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.a = self.root/'first'
        self.b = self.root/'second'
        self.fixture = label_fixture(self.a)
        label_fixture(self.b, np.array([[0., 0., -1.], [0., 1., 0.], [1., 0., 0.], [0., 0., 1.]]))

    def test_packages_are_bound_and_fields_convert_only_on_request(self):
        data = load_dataset([self.a, self.b], [1., 3.])
        np.testing.assert_array_equal(data.weights, [.25, .75])
        np.testing.assert_array_equal(data.grid, [4, 4, 2])
        self.assertIsInstance(data.samples[0].scalar_Ry, np.memmap)
        self.assertIsInstance(data.samples[0].B_Ry, np.memmap)
        expected = self.fixture['fields'].copy()
        expected[:, 0] -= expected[:, 0].mean()
        np.testing.assert_allclose(data.samples[0].field_block(4, 9), expected[4:9]*RY_EV)
        self.assertEqual(data.identity['samples'][0]['manifest_sha256'], sha256(self.a/'LABELS.json'))
        self.assertEqual(set(data.identity['samples'][0]['files']),
                         {'scalar.npy', 'B.npy', 'N.npy', 'moments.npy', 'atom_ids.npy', 'grid.npy', 'structure.npz'})

    def test_initial_reference_recovers_spin_aligned_d4h_reference(self):
        data = load_dataset([self.a, self.b], [2., 1.])
        actual = initial_reference(data, 2)
        expected = self.fixture['reference']
        np.testing.assert_allclose(actual.n_scalar, expected.n_scalar, atol=2e-16)
        np.testing.assert_allclose(actual.n_spin, expected.n_spin, atol=2e-16)
        np.testing.assert_array_equal(actual.p_scalar, [0., 0.])
        np.testing.assert_array_equal(actual.p_spin, [0., 0.])
        self.assertTrue(actual.source.startswith('label-set-sha256:'))
        n = actual.values(data.samples[0].structure.spins, np.ones((2, 4)))[0]
        np.testing.assert_allclose(n, data.samples[0].N, atol=2e-16)

    def test_changed_payload_and_invalid_metadata_are_rejected(self):
        path = self.a/'B.npy'
        values = np.load(path)
        values[0, 1] += .1
        np.save(path, values)
        with self.assertRaisesRegex(ValueError, 'checksum'):
            load_dataset([self.a])
        rewrite_array(self.a, 'B.npy', values)
        manifest_path = self.a/'LABELS.json'
        original = json.loads(manifest_path.read_text())
        for key, value in [('schema', 'feedback_labels'), ('units', {'scalar': 'eV'}),
                           ('grid', [8, 4, 2]), ('atom_ids', [1, 2, 3, 4])]:
            with self.subTest(key=key):
                manifest = dict(original)
                manifest[key] = value
                manifest_path.write_text(json.dumps(manifest))
                with self.assertRaises(ValueError):
                    load_dataset([self.a])
        manifest_path.write_text(json.dumps(original))
        load_dataset([self.a])

    def test_physical_invalidity_is_rejected_even_with_correct_hash(self):
        originals = {name: np.load(self.a/name).copy() for name in ('scalar.npy', 'B.npy', 'N.npy', 'moments.npy')}
        cases = []
        bad = originals['scalar.npy'] + 1.
        cases.append(('scalar.npy', bad))
        bad = originals['B.npy'].copy(); bad[0, 0] = np.nan
        cases.append(('B.npy', bad))
        bad = originals['N.npy'].copy(); bad[0, 0, 1] = .1j
        cases.append(('N.npy', bad))
        bad = originals['N.npy'].copy(); bad[0, 0, 0] = 1.2
        cases.append(('N.npy', bad))
        bad = originals['moments.npy'].copy(); bad[0] = 0.
        cases.append(('moments.npy', bad))
        bad = originals['moments.npy'].copy(); bad[0] *= -1
        cases.append(('moments.npy', bad))
        for name, value in cases:
            with self.subTest(file=name):
                rewrite_array(self.a, name, value)
                with self.assertRaises(ValueError):
                    load_dataset([self.a])
                rewrite_array(self.a, name, originals[name])

    def test_common_geometry_and_weights_are_required(self):
        path = self.b/'structure.npz'
        with np.load(path) as values:
            structure = {key: values[key].copy() for key in values.files}
        structure['positions'][0, 0] += .1
        np.savez(path, **structure)
        manifest = json.loads((self.b/'LABELS.json').read_text())
        manifest['files']['structure.npz'] = sha256(path)
        (self.b/'LABELS.json').write_text(json.dumps(manifest))
        with self.assertRaisesRegex(ValueError, 'positions differ'):
            load_dataset([self.a, self.b])
        for weights in ([0.], [-1.], [np.nan], [1., 2.]):
            with self.subTest(weights=weights):
                with self.assertRaises(ValueError):
                    load_dataset([self.a], weights)
        with self.assertRaises(ValueError):
            load_dataset([self.a, self.a])
        with self.assertRaises(ValueError):
            load_dataset([])


if __name__ == '__main__':
    unittest.main()
