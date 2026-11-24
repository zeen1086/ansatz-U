"""Small end-to-end fits of synthetic labels using only the explicit response."""
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

import fit
from test_training_data import label_fixture, rewrite_array
from training_data import load_dataset
from workflow import SINGLE_RESPONSE_SCHEMA
from predict import Bundle
from symmetry import matrix_sigmoid

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import local_vfunc


class FitTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.a, self.b = self.root/'first', self.root/'second'
        label_fixture(self.a)
        label_fixture(self.b, np.array([[0., 0., -1.], [0., 1., 0.], [1., 0., 0.], [0., 0., 1.]]))
        self.config = self.root/'fit.toml'
        self.config.write_text('''[data]
labels = ["first", "second"]
weights = [1.0, 2.0]
[features]
radial_lengths_A = [1.5]
scalar_degree = 1
vector_degree = 1
[kernels]
radial_lengths_A = [1.5]
[loss]
scalar_scale_eV = 1.0
magnetic_scale_eV = 1.0
occupation_scale = 0.1
regularization = 1e-8
[optimizer]
max_iterations = 4
gradient_tolerance = 1e-8
loss_tolerance = 1e-14
checkpoint_every = 1
seed = 4
''')

    def problem(self):
        config = fit.read_config(self.config)
        dataset = load_dataset([Path(p) for p in config['data']['labels']], config['data']['weights'])
        return fit.FitProblem(dataset, config)

    def test_joint_gradient_and_two_sigmoid_prediction(self):
        problem = self.problem()
        rng = np.random.default_rng(234)
        vector = problem.pack()
        vector += rng.normal(scale=.005, size=vector.shape)
        loss, gradient, _ = problem.evaluate(vector)
        self.assertTrue(np.isfinite(loss))
        groups = ['model/A_se', 'model/A_vo', 'model/K_se', 'model/K_vo',
                  'model/Q_se', 'model/Q_so', 'model/Q_vo', 'model/Q_ve',
                  'model/p_se', 'model/p_vo', 'kernel/W_v', 'kernel/W_B', 'kernel/bias']
        for key in groups:
            with self.subTest(parameter=key):
                directions = {k: np.zeros(shape) for k, shape in problem.shapes.items()}
                directions[key] = rng.normal(size=problem.shapes[key])
                direction = problem.pack(directions)
                direction /= np.linalg.norm(direction)
                step = 2e-5
                plus = problem.evaluate(vector + step*direction)[0]
                minus = problem.evaluate(vector - step*direction)[0]
                np.testing.assert_allclose(gradient @ direction, (plus-minus)/(2*step),
                                           rtol=3e-4, atol=2e-7)
        problem.unpack(vector)
        with patch('model.matrix_sigmoid', wraps=matrix_sigmoid) as sigmoid:
            problem.model.predict(problem.dataset.samples[0].structure,
                                  prepared_features=problem.phi[0], normalized=True)
        self.assertEqual(sigmoid.call_count, 2)
        self.assertFalse(any(key.startswith('L_') for key in problem.model.params))

    def test_short_training_export_and_parameter_resume(self):
        output = self.root/'training'
        report = fit.train(self.config, output)
        self.assertLess(report['final_loss'], report['initial_loss'])
        self.assertEqual(report['model_schema'], SINGLE_RESPONSE_SCHEMA)
        self.assertFalse(report['physical_release'])
        bundle = Path(report['bundle'])
        for name in local_vfunc.DATA:
            self.assertTrue((bundle/name).is_file(), name)
        self.assertEqual(set(report['bundle_files']), {*local_vfunc.DATA, 'structure.npz'})
        self.assertEqual(report['bundle_files'],
                         {name: fit.digest(bundle/name) for name in report['bundle_files']})
        meta = json.loads((bundle/'model/model.json').read_text())
        self.assertEqual(meta['schema'], SINGLE_RESPONSE_SCHEMA)
        with np.load(bundle/'model/parameters.npz') as parameters:
            self.assertFalse(any(key.startswith('L_') for key in parameters.files))
            self.assertIn('A_se', parameters.files)
            self.assertIn('A_vo', parameters.files)
        native = self.root/'exported'
        local_vfunc.export(bundle, ROOT/'4_potential/local', bundle/'structure.npz', native, 0.)
        before, after = Bundle(bundle), Bundle(native)
        structure = load_dataset([self.a]).samples[0].structure
        p0, p1 = before.model.predict(structure), after.model.predict(structure)
        np.testing.assert_array_equal(p0.p, p1.p)
        np.testing.assert_array_equal(p0.N, p1.N)
        np.testing.assert_array_equal(before.cv, after.cv)
        np.testing.assert_array_equal(before.cb, after.cb)
        with self.assertRaises(FileExistsError):
            fit.train(self.config, output)
        resumed = fit.train(self.config, output, resume=True)
        np.testing.assert_allclose(resumed['initial_loss'], report['final_loss'], rtol=0, atol=1e-12)
        self.assertGreaterEqual(resumed['iterations'], report['iterations'])
        self.assertLessEqual(resumed['final_loss'], resumed['initial_loss'] + 1e-12)

    def test_resume_rejects_old_schema_or_changed_inputs(self):
        output = self.root/'training'
        fit.train(self.config, output)
        checkpoint = output/'checkpoint.npz'
        with np.load(checkpoint, allow_pickle=False) as data:
            vector = data['parameters'].copy()
            metadata = json.loads(str(data['metadata']))
        changed = vector.copy()
        changed[0] += .01
        np.savez(checkpoint, parameters=changed, metadata=json.dumps(metadata))
        with self.assertRaisesRegex(ValueError, 'parameter checksum'):
            fit.train(self.config, output, resume=True)
        np.savez(checkpoint, parameters=vector, metadata=json.dumps(metadata))
        old = json.loads(json.dumps(metadata))
        old['identity']['model_schema'] = 'LANIO2_ACTIVE_LOCAL_GRID_FEEDBACK_V2'
        np.savez(checkpoint, parameters=vector, metadata=json.dumps(old))
        with self.assertRaisesRegex(ValueError, 'Checkpoint model schema'):
            fit.train(self.config, output, resume=True)
        np.savez(checkpoint, parameters=vector, metadata=json.dumps(metadata))
        config = self.config.read_text()
        self.config.write_text(config.replace('seed = 4', 'seed = 5'))
        with self.assertRaisesRegex(ValueError, 'unchanged'):
            fit.train(self.config, output, resume=True)
        self.config.write_text(config)
        b = np.load(self.a/'B.npy').copy()
        b[0, 0] += .01
        rewrite_array(self.a, 'B.npy', b)
        with self.assertRaisesRegex(ValueError, 'unchanged'):
            fit.train(self.config, output, resume=True)

    def test_resume_discards_iterations_after_last_checkpoint(self):
        self.config.write_text(self.config.read_text().replace('checkpoint_every = 1', 'checkpoint_every = 3'))
        output = self.root/'interrupted'
        initial = None

        def interrupted(function, vector, *, callback, **kwargs):
            nonlocal initial
            initial = vector.copy()
            callback(vector + .0001)
            callback(vector + .0002)
            raise KeyboardInterrupt

        with patch('fit.minimize', side_effect=interrupted):
            with self.assertRaises(KeyboardInterrupt):
                fit.train(self.config, output)
        history = output/'history.jsonl'
        self.assertEqual([json.loads(row)['iteration'] for row in history.read_text().splitlines()], [1, 2])
        with np.load(output/'checkpoint.npz', allow_pickle=False) as checkpoint:
            self.assertEqual(json.loads(str(checkpoint['metadata']))['iteration'], 0)
            np.testing.assert_array_equal(checkpoint['parameters'], initial)

        def continued(function, vector, *, callback, **kwargs):
            np.testing.assert_array_equal(vector, initial)
            callback(vector)
            return SimpleNamespace(x=vector, nit=1, success=True, status=0, message='synthetic single accepted step')

        with patch('fit.minimize', side_effect=continued):
            report = fit.train(self.config, output, resume=True)
        self.assertEqual(report['iterations'], 1)
        self.assertEqual([json.loads(row)['iteration'] for row in history.read_text().splitlines()], [1])


if __name__ == '__main__':
    unittest.main()
