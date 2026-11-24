"""Generated model and projection records contain identities and sizes only."""
import json
from pathlib import Path
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


if __name__=='__main__':unittest.main()
