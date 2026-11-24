"""Generated model and projection records contain identities and sizes only."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
from types import SimpleNamespace
import tempfile
import unittest
import numpy as np
import local_vfunc as api


class ExportRecordTests(unittest.TestCase):
    def test_export_and_projection_records_without_content_digests(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'source';code=root/'code'
            (source/'model').mkdir(parents=True);code.mkdir()
            for name in api.RUNTIME:(code/name).write_text('# fixture\n')
            (source/'model/model.json').write_text(json.dumps(dict(schema='LANIO2_ACTIVE_LOCAL_GRID_FEEDBACK_V2',p_scales=[1.],ueff=1.,scalar_gauge='cell_mean_zero')))
            np.savez(source/'model/parameters.npz',p=np.zeros(1))
            np.save(source/'C_v_local.npy',np.zeros((1,2)))
            np.save(source/'C_B_local.npy',np.zeros((1,1)))
            np.savez(source/'local_grid.npz',nodes=np.zeros((1,3),int),spacing_A=np.ones(3),half=np.zeros(3,int))
            (source/'KERNELS.json').write_text(json.dumps(dict(cutoff_A=6.,scalar_bias_coefficient_per_Ni=0.)))
            structure=root/'structure.npz'
            np.savez(structure,cell=np.diag([2.,2.,1.]),positions=np.array([[0.,0.,0.],[0.,1.,0.],[1.,0.,0.],[1.,1.,0.]]),species=np.full(4,28),ni=np.arange(4))
            model=root/'model';result=api.export(source,code,structure,model,0.)
            meta=api.verify(model)
            self.assertEqual(len(result['model_id']),64)
            self.assertEqual(set(meta),{'schema','prediction_rule','files','primitive_cell_A','Ni_fractional','primitive_motif','primitive_grid','scalar_offset_eV','scalar_gauge','Ueff_eV','ncoeff','spin_order','Hamiltonian_order','no_refit','model_id','source_bundle','source_code','native_geometry_bytes'})
            self.assertEqual(api.runtime_path(model,meta), (model/'runtime').resolve())
            self.assertEqual(set(p.name for p in (model/'runtime').iterdir()), set(api.RUNTIME))
            self.assertNotIn('runtime/predict.py', meta['files'])
            self.assertEqual((model/'runtime/model.py').read_bytes(),(api.LOCAL_CODE/'model.py').read_bytes())
            self.assertEqual(meta['prediction_rule'],'N_single_response_to_v_B')
            for rule in (None, 'unsupported'):
                invalid = meta.copy()
                if rule is None:invalid.pop('prediction_rule')
                else:invalid['prediction_rule']=rule
                (model/'V_FUNC.json').write_text(json.dumps(invalid))
                with self.assertRaisesRegex(ValueError,'Unsupported model prediction rule'):
                    api.metadata(model)
            (model/'V_FUNC.json').write_text(json.dumps(meta))
            with self.assertRaisesRegex(ValueError,'Incomplete model runtime'):
                api.runtime_path(root/'missing-runtime',meta)
            self.assertTrue(all(set(v)=={'bytes'} for v in meta['files'].values()))
            bg=root/'background';bg.write_bytes(b'fixture')
            Path(str(bg)+'.json').write_text(json.dumps(dict(kind='pure_H0_no_local_fields_no_HU',background_id='fixture')))
            wan=root/'wan';wan.mkdir();np.ones(1,dtype='<c16').tofile(wan/'wan_1')
            translations=root/'translations';translations.write_text('1\n0 0 0 1\n')
            args=SimpleNamespace(output=root/'projection',model=model,background=bg,translations=translations,kz=[0.],grid_weight=1.,nwan=1,wan_dir=wan)
            api.prepare_projection(args)
            record=json.loads((args.output/'PROJECTION.json').read_text())
            self.assertEqual(set(record),{'model_id','Wannier_files','translations_bytes','kz','background_id','grid_weight','nwan','source_background','basis_id','status','native_projection','projected_center_translations','projected_payload_bytes'})
            self.assertEqual(len(record['basis_id']),64)
            source_meta=json.loads((source/'model/model.json').read_text())
            source_meta['schema']='LANIO2_ACTIVE_LOCAL_GRID_SINGLE_RESPONSE_V1'
            (source/'model/model.json').write_text(json.dumps(source_meta))
            api.export(source,code,structure,root/'single-response',0.)
            source_meta['schema']='unsupported'
            (source/'model/model.json').write_text(json.dumps(source_meta))
            with self.assertRaisesRegex(ValueError,'Unsupported source model schema'):
                api.export(source,code,structure,root/'unsupported',0.)
            self.assertFalse((root/'unsupported').exists())
            (model/'C_v_local.npy').write_bytes(b'broken')
            with self.assertRaisesRegex(ValueError,'truncated'):api.verify(model)


class RuntimeLayoutTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        sys.path.insert(0, str(api.LOCAL_CODE))
        from features import FeatureSpec
        from model import LocalModel, Reference
        from workflow import save_model

        source = self.root/'source'
        ref = Reference(np.zeros(15), np.zeros(15), np.array([.2]),
                        np.array([-.1]), 'synthetic runtime-layout regression')
        model = LocalModel(FeatureSpec(2.4, (1.2,)), np.array([[1.5, .8]]), ref,
                           p_representation=np.ones((16, 1, 1)), seed=7)
        model.params['K_se'][0, 0] = .7
        model.params['A_se'] = .3*np.eye(15)
        save_model(model, source/'model')
        np.save(source/'C_v_local.npy', np.zeros((1, 2)))
        np.save(source/'C_B_local.npy', np.zeros((1, 1)))
        np.savez(source/'local_grid.npz', nodes=np.zeros((1, 3), int),
                 spacing_A=np.ones(3), half=np.zeros(3, int))
        (source/'KERNELS.json').write_text(json.dumps(dict(
            cutoff_A=6., scalar_bias_coefficient_per_Ni=0.)))
        self.structure = self.root/'structure.npz'
        np.savez(self.structure, cell=np.diag([4., 4., 4.]),
                 positions=np.array([[0., 0., 0.], [0., 2., 0.],
                                     [2., 0., 0.], [2., 2., 0.]]),
                 species=np.full(4, 28), ni=np.arange(4),
                 spins=np.array([[0., 0., 1.], [.6, .8, 0.],
                                 [0., -.6, .8], [-.8, 0., .6]]))
        self.current, self.legacy = self.root/'current', self.root/'legacy'
        api.export(source, api.LOCAL_CODE, self.structure, self.current, 0.)
        shutil.copytree(self.current, self.legacy)
        workflow = self.legacy/'runtime/workflow.py'
        workflow.write_text(workflow.read_text().split('\nclass Bundle:', 1)[0])
        predictor = self.legacy/'runtime/predict.py'
        predictor.write_text('''from pathlib import Path
import json
import numpy as np
from workflow import load_model
from symmetry import require

class Bundle:
    def __init__(self, path):
        path = Path(path)
        self.model = load_model(path/'model')
        self.cv = np.load(path/'C_v_local.npy', mmap_mode='r', allow_pickle=False)
        self.cb = np.load(path/'C_B_local.npy', mmap_mode='r', allow_pickle=False)
        with np.load(path/'local_grid.npz', allow_pickle=False) as data:
            self.geo = {key: data[key].copy() for key in data.files}
        self.alpha = json.loads((path/'KERNELS.json').read_text())['scalar_bias_coefficient_per_Ni']
        require(self.cv.shape == (len(self.geo['nodes']), self.model.np+1), 'scalar kernel shape')
        require(self.cb.shape == (len(self.geo['nodes']), self.model.np), 'spin kernel shape')
''')
        meta = api.metadata(self.legacy)
        for path in (workflow, predictor):
            meta['files']['runtime/'+path.name] = {'bytes': path.stat().st_size}
        (self.legacy/'V_FUNC.json').write_text(json.dumps(meta))

    def test_old_and_merged_runtime_predict_identically_and_reject_mixing(self):
        script = '''from pathlib import Path
import sys
import numpy as np
import local_vfunc as api
package, other, structure, output = map(Path, sys.argv[1:])
loader, features = api.runtime(package)
assert loader.__name__ == ('predict' if package.name == 'legacy' else 'workflow')
if package.name == 'legacy':
    import workflow
    assert not hasattr(workflow, 'Bundle')
bundle = loader.Bundle(package)
with np.load(structure) as data:
    s = features.Structure(**{key: data[key] for key in data.files})
prediction = bundle.model.predict(s)
np.savez(output, p=bundle.model.physical_p(prediction), N=prediction.N)
try:
    api.runtime(other)
except RuntimeError as error:
    assert 'Conflicting model module' in str(error)
else:
    raise AssertionError('Mixed runtime modules were accepted')
'''
        for package, other in ((self.current, self.legacy), (self.legacy, self.current)):
            api.verify(package)
            subprocess.run([sys.executable, '-B', '-c', script, str(package), str(other),
                            str(self.structure), str(self.root/(package.name+'.npz'))],
                           cwd=Path(api.__file__).resolve().parent, check=True,
                           capture_output=True, text=True)
        with np.load(self.root/'current.npz') as current, np.load(self.root/'legacy.npz') as legacy:
            for key in ('p', 'N'):
                self.assertGreater(np.linalg.norm(current[key]), 0.)
                np.testing.assert_array_equal(current[key], legacy[key])

    def test_declared_legacy_loader_is_required(self):
        current_meta, legacy_meta = api.metadata(self.current), api.metadata(self.legacy)
        self.assertEqual(api.runtime_names(current_meta), api.RUNTIME)
        self.assertEqual(api.runtime_names(legacy_meta), api.RUNTIME+('predict.py',))
        (self.legacy/'runtime/predict.py').unlink()
        with self.assertRaisesRegex(ValueError, 'Incomplete model runtime'):
            api.runtime_path(self.legacy, legacy_meta)


if __name__=='__main__':unittest.main()
