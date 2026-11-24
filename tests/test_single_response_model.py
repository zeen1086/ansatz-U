"""One effective occupation response, its gradients and legacy weight import."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'4_potential/local'))
from features import FeatureSpec, Structure, features
from model import LocalModel, Reference
from symmetry import REP, SECTORS, matrix_sigmoid, to_channels, to_matrix
from workflow import SINGLE_RESPONSE_SCHEMA, load_model, save_model


class SingleResponseModelTests(unittest.TestCase):
    def setUp(self):
        spins = np.array([[0., 0., 1.], [.6, .8, 0.], [.2, -.3, .9], [-.3, .8, .5]])
        spins /= np.linalg.norm(spins, axis=1)[:, None]
        self.s = Structure(np.eye(3)*6.,
                           np.array([[0., 0., 0.], [1.1, .2, 0.], [.2, 1.2, .4], [.3, .2, 1.3]]),
                           np.full(4, 28), np.arange(4), spins)
        n_scalar, n_spin = np.zeros(15), np.zeros(15)
        n_scalar[0], n_spin[0] = .1, .05
        self.ref = Reference(n_scalar, n_spin, np.array([.2]), np.array([-.1]), 'synthetic')
        self.m = LocalModel(FeatureSpec(2.4, (1.2,)), np.array([[1.5, .8]]), self.ref,
                            p_representation=np.ones((16, 1, 1)), ueff=2.25, seed=7)
        self.m.params['K_se'][0, 0] = .7
        self.m.params['K_vo'][0, 0] = -.4
        rng = np.random.default_rng(29)
        for sector, diagonal in (('se', .3), ('vo', -.4)):
            self.m.params['A_' + sector] = diagonal*np.eye(15) + rng.normal(scale=.1, size=(15, 15))

    def test_two_sigmoids_and_one_effective_response(self):
        m, s = self.m, self.s
        phi = m.prepare_features(features(s, m.spec))
        q0 = to_matrix(m._head('Q', phi, s.spins))
        p0 = m.pack_p(m._head('p', phi, s.spins))
        n0 = matrix_sigmoid(q0)[0]
        nref, _ = self.ref.values(s.spins, m.scale)
        delta0 = to_channels(n0 - nref)
        response = {k: np.zeros((len(s.ni), a, c)) for k, (a, c) in SECTORS.items()}
        for sector in ('se', 'vo'):
            rep = REP[sector]
            amap = sum(g @ m.params['A_' + sector] @ g.T for g in rep) / len(rep)
            response[sector] = np.einsum('ab,nbc->nac', amap, delta0[sector])
        n1 = matrix_sigmoid(q0 + to_matrix(response))[0]
        delta1 = to_channels(n1 - nref)
        correction = np.concatenate((.7*delta1['se'][:, :1], -.4*delta1['vo'][:, :1]), axis=-1)
        with patch('model.matrix_sigmoid', wraps=matrix_sigmoid) as sigmoid:
            pred = m.predict(s, max_iterations=0)
        self.assertEqual(sigmoid.call_count, 2)
        np.testing.assert_allclose(pred.N, n1, atol=2e-15)
        np.testing.assert_allclose(pred.p, p0 + correction, atol=2e-15)
        np.testing.assert_allclose(pred.VU, m.ueff*(.5*np.eye(10)-n1))
        self.assertGreater(np.max(abs(n1-n0)), 1e-3)
        self.assertEqual(pred.iterations, 0)
        self.assertTrue(pred.converged)
        m.params['b_p_se'] += 3.
        m.params['K_se'] *= 2.
        changed = m.predict(s)
        np.testing.assert_array_equal(changed.N, pred.N)
        self.assertGreater(np.max(abs(changed.p-pred.p)), 1.)

    def test_full_complex_occupation_bounds_and_spin_symmetries(self):
        pred = self.m.predict(self.s)
        np.testing.assert_allclose(pred.N, pred.N.conj().swapaxes(-1, -2), atol=2e-15)
        eigenvalues = np.linalg.eigvalsh(pred.N)
        self.assertTrue(np.all((eigenvalues > 0) & (eigenvalues < 1)))
        channels = to_channels(pred.N)
        for sector in ('so', 've'):
            self.assertGreater(np.linalg.norm(channels[sector]), 1e-7)
        rotated = self.s.spins[:, [1, 0, 2]] * [-1., 1., 1.]
        u = np.kron(np.diag(np.exp([-1j*np.pi/4, 1j*np.pi/4])), np.eye(5))
        time_reversal = np.kron(np.array([[0., 1.], [-1., 0.]]), np.eye(5))
        for spins, expected in ((rotated, u @ pred.N @ u.conj().T),
                                (-self.s.spins, time_reversal @ pred.N.conj() @ time_reversal.T)):
            s = Structure(self.s.cell, self.s.positions, self.s.species, self.s.ni, spins)
            np.testing.assert_allclose(self.m.predict(s).N, expected, atol=3e-14)

    def test_field_and_occupation_gradients_through_both_sigmoids(self):
        m, s = self.m, self.s
        raw = features(s, m.spec)
        pred = m.predict(s, prepared_features=raw)
        rng = np.random.default_rng(41)
        gp = rng.normal(size=pred.p.shape)
        gn = rng.normal(size=pred.N.shape) + 1j*rng.normal(size=pred.N.shape)
        gn = (gn + gn.conj().swapaxes(-1, -2))/2
        gradients = m.backward(pred, gp, gn)

        def loss():
            value = m.predict(s, prepared_features=raw)
            return np.vdot(gp, value.p).real + np.vdot(gn, value.N).real

        for key in ('A_se', 'A_vo', 'Q_se', 'Q_so', 'Q_vo', 'Q_ve', 'K_se', 'K_vo',
                    'b_Q_se', 'b_Q_vo', 'p_se', 'p_vo', 'b_p_se', 'b_p_vo'):
            with self.subTest(parameter=key):
                original = m.params[key].copy()
                direction = rng.normal(size=original.shape)
                direction /= np.linalg.norm(direction)
                analytical = np.vdot(gradients[key], direction).real
                h = 1e-5
                try:
                    m.params[key] = original + h*direction
                    plus = loss()
                    m.params[key] = original - h*direction
                    minus = loss()
                finally:
                    m.params[key] = original
                self.assertGreater(abs(analytical), 1e-9)
                np.testing.assert_allclose(analytical, (plus-minus)/(2*h), rtol=2e-5, atol=2e-9)

    def legacy_archive(self, directory, schema):
        save_model(self.m, directory)
        meta = json.loads((directory/'model.json').read_text())
        meta.update(schema=schema, rho=.25, feedback_coordinates={'se': [[2.]], 'vo': [[.5]]})
        (directory/'model.json').write_text(json.dumps(meta))
        weights = {k: v.copy() for k, v in self.m.params.items() if not k.startswith('A_')}
        for sector, kvalue, lvalue in (('se', 8., 5.), ('vo', 3., -2.)):
            weights['K_' + sector] = np.zeros((1, 15))
            weights['L_' + sector] = np.zeros((15, 1))
            weights['K_' + sector][0, 0] = kvalue
            weights['L_' + sector][0, 0] = lvalue
            # The dz2-dxz orbital channel averages to zero under local D4h.
            weights['K_' + sector][0, 5] = 7.
            weights['L_' + sector][5, 0] = -4.
        np.savez(directory/'parameters.npz', **weights)
        norm_k = np.sqrt((2*8.)**2 + 3*(.5*3.)**2)
        norm_l = np.sqrt((5./2)**2 + 3*(-2./.5)**2)
        return weights, np.sqrt(4*.25/(norm_k*norm_l))

    def test_coupled_import_retains_scaled_LK_and_roundtrips(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp)/'old'
            weights, scale = self.legacy_archive(source, 'LANIO2_ACTIVE_LOCAL_GRID_FEEDBACK_V2')
            imported = load_model(source)
            self.assertFalse(any(k.startswith('L_') for k in imported.params))
            for sector in ('se', 'vo'):
                np.testing.assert_allclose(imported.params['K_' + sector], scale*weights['K_' + sector])
                expected_a = np.zeros((15, 15))
                expected_a[0, 0] = scale**2*weights['L_' + sector][0, 0]*weights['K_' + sector][0, 0]
                np.testing.assert_allclose(imported.params['A_' + sector], expected_a, atol=1e-15)
                np.testing.assert_array_equal(imported.params['Q_' + sector], weights['Q_' + sector])
            target = Path(tmp)/'new'
            save_model(imported, target)
            self.assertEqual(json.loads((target/'model.json').read_text())['schema'], SINGLE_RESPONSE_SCHEMA)
            saved = load_model(target)
            for key in imported.params:
                np.testing.assert_array_equal(saved.params[key], imported.params[key])
            before, after = imported.predict(self.s), saved.predict(self.s)
            np.testing.assert_array_equal(after.N, before.N)
            np.testing.assert_array_equal(after.p, before.p)

    def test_warm_import_preserves_compensated_head(self):
        for schema in ('LANIO2_GEOMETRIC_FEEDBACK_V1', 'LANIO2_SEPARATE_LOCAL_GRID_FEEDBACK_V1'):
            with self.subTest(schema=schema), tempfile.TemporaryDirectory() as tmp:
                source = Path(tmp)/'warm'
                weights, scale = self.legacy_archive(source, schema)
                with self.assertRaises(ValueError):
                    load_model(source)
                imported = load_model(source, warm=True)
                for sector, column, reference in (('se', 0, self.ref.p_scalar), ('vo', 1, self.ref.p_spin)):
                    lmap = np.zeros_like(weights['L_' + sector])
                    lmap[0, 0] = scale*weights['L_' + sector][0, 0]
                    np.testing.assert_allclose(imported.params['Q_' + sector],
                                               weights['Q_' + sector] + lmap @ weights['p_' + sector])
                    np.testing.assert_allclose(imported.params['b_Q_' + sector],
                                               weights['b_Q_' + sector] + lmap @ (
                                                   weights['b_p_' + sector] - reference/self.m.scale[:, column]))
                    expected_a = np.zeros((15, 15))
                    expected_a[0, 0] = scale**2*weights['L_' + sector][0, 0]*weights['K_' + sector][0, 0]
                    np.testing.assert_allclose(imported.params['A_' + sector], expected_a, atol=1e-15)
                meta = json.loads((source/'model.json').read_text())
                meta['p_feature_maps'] = {}
                for sector in ('se', 'vo'):
                    size = weights['p_' + sector].shape[1]
                    meta['p_feature_maps'][sector] = {
                        'plus': np.eye(size).tolist(), 'minus': np.zeros((size, 0)).tolist()}
                (source/'model.json').write_text(json.dumps(meta))
                with self.assertRaisesRegex(ValueError, 'p_feature_maps cannot be compensated exactly'):
                    load_model(source, warm=True)

    def test_saved_one_way_requires_original_coupled_archive(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp)/'one_way'
            save_model(self.m, source)
            meta = json.loads((source/'model.json').read_text())
            meta['schema'] = 'LANIO2_ACTIVE_LOCAL_GRID_ONE_WAY_V1'
            (source/'model.json').write_text(json.dumps(meta))
            np.savez(source/'parameters.npz', **{k: v for k, v in self.m.params.items() if not k.startswith('A_')})
            with self.assertRaisesRegex(ValueError, '(?i)(coupled|feedback)'):
                load_model(source)


if __name__ == '__main__':
    unittest.main()
