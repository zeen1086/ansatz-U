"""Direct complex occupation, its one-way field response and frozen weight import."""
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
from symmetry import matrix_sigmoid, to_channels, to_matrix
from workflow import load_model, save_model


class OneWayModelTests(unittest.TestCase):
    def setUp(self):
        self.s = Structure(np.eye(3)*6., np.array([[0.,0.,0.], [1.1,.2,0.]]),
                           np.array([28,28]), np.array([0,1]), np.array([[0.,0.,1.], [.6,.8,0.]]))
        self.ref = Reference(np.zeros(15), np.zeros(15), np.zeros(1), np.zeros(1), 'synthetic')
        self.m = LocalModel(FeatureSpec(2.4,(1.2,)), np.array([[1.5,.8]]), self.ref,
                            p_representation=np.ones((16,1,1)), ueff=2.25, seed=7)
        self.m.params['K_se'][0,0] = .7
        self.m.params['K_vo'][0,0] = -.4

    def test_one_sigmoid_and_one_field_correction(self):
        m, s = self.m, self.s
        phi = m.prepare_features(features(s, m.spec))
        q0 = to_matrix(m._head('Q', phi, s.spins))
        p0 = m.pack_p(m._head('p', phi, s.spins))
        n = matrix_sigmoid(q0)[0]
        nref, _ = self.ref.values(s.spins, m.scale)
        delta = to_channels(n - nref)
        correction = np.concatenate((.7*delta['se'][:,:1], -.4*delta['vo'][:,:1]), axis=-1)
        with patch('model.matrix_sigmoid', wraps=matrix_sigmoid) as sigmoid:
            pred = m.predict(s, max_iterations=0)
        sigmoid.assert_called_once()
        np.testing.assert_allclose(pred.N, n)
        np.testing.assert_allclose(pred.p, p0 + correction)
        np.testing.assert_allclose(pred.VU, m.ueff*(.5*np.eye(10)-n))
        self.assertGreater(np.max(abs(pred.p-p0)), 1e-3)
        self.assertEqual(pred.iterations, 0)
        m.params['b_p_se'] += 3.
        m.params['K_se'] *= 2.
        changed = m.predict(s)
        np.testing.assert_array_equal(changed.N, pred.N)
        self.assertGreater(np.max(abs(changed.p-pred.p)), 1.)

    def test_field_loss_reaches_occupation_head_through_K(self):
        m, s = self.m, self.s
        pred = m.predict(s)
        gp = np.ones_like(pred.p)
        gradients = m.backward(pred, gp, np.zeros_like(pred.N))
        for key, index in [('b_p_se', (0,)), ('b_Q_se', (0,)), ('K_se', (0,0))]:
            with self.subTest(key=key):
                h = 1e-5
                original = m.params[key][index]
                m.params[key][index] = original + h
                plus = m.predict(s).p.sum()
                m.params[key][index] = original - h
                minus = m.predict(s).p.sum()
                m.params[key][index] = original
                self.assertGreater(abs(gradients[key][index]), 1e-4)
                self.assertAlmostEqual(gradients[key][index], (plus-minus)/(2*h), places=8)

    def test_coupled_archive_keeps_its_scaled_K_and_drops_L(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp)/'old'
            save_model(self.m, source)
            meta = json.loads((source/'model.json').read_text())
            meta.update(schema='LANIO2_ACTIVE_LOCAL_GRID_FEEDBACK_V2', rho=.25)
            (source/'model.json').write_text(json.dumps(meta))
            weights = self.m.params.copy()
            weights.update(K_se=np.zeros((1,15)), K_vo=np.zeros((1,15)),
                           L_se=np.zeros((15,1)), L_vo=np.zeros((15,1)))
            weights['K_se'][0,0], weights['L_se'][0,0] = 8., 5.
            np.savez(source/'parameters.npz', **weights)
            imported = load_model(source)
            self.assertFalse(any(k.startswith('L_') for k in imported.params))
            self.assertAlmostEqual(imported.params['K_se'][0,0], 8.*np.sqrt(4*.25/(8.*5.)))
            np.testing.assert_array_equal(imported.predict(self.s).N, self.m.predict(self.s).N)
            target = Path(tmp)/'new'
            save_model(imported, target)
            np.testing.assert_array_equal(load_model(target).predict(self.s).p, imported.predict(self.s).p)


if __name__ == '__main__':
    unittest.main()
