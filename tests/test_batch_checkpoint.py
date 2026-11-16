"""Check that a restart reuses completed points and never publishes a partial batch."""
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import unittest
from unittest.mock import patch
import numpy as np

import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'7_spin'))
import run_material


class BatchCheckpointTest(unittest.TestCase):
    def setUp(self):
        tmp=TemporaryDirectory();self.addCleanup(tmp.cleanup)
        self.root=Path(tmp.name);(self.root/'RVR').mkdir()
        (self.root/'RVR/PROJECTION.json').write_text(json.dumps({'basis_id':'fixture'}))
        working=patch.object(run_material,'WORK',self.root);working.start();self.addCleanup(working.stop)
        (self.root/'coeff').write_bytes(b'coefficients')
        (self.root/'RVR/B_general').write_bytes(b'rvr')
        (self.root/'basis').mkdir();(self.root/'basis/pure_H0.au').write_bytes(b'background')
        (self.root/'local').mkdir();(self.root/'local/cal_Hsk').write_bytes(b'binary')
        build=patch.object(run_material,'BUILD',self.root);build.start();self.addCleanup(build.stop)
        self.k=np.array([[.1,.2,.125,.25],[.3,.4,.125,.75]])
        self.out=[self.root/f'k_{i:04}' for i in range(2)]

    def staged(self,coeff,k,out,components,truncate=False):
        lines=['LANIO2_LOCAL_H_V2',f'1 2 {len(k)} F',run_material.CONFIG['model_id'],'fixture']
        for slot,point in enumerate(k,1):
            lines.append(' '.join(map(str,[slot,*point[:3],1/len(k),0.,0.,1.])))
            h=np.eye(2,dtype='<c16')*point[0]
            h.ravel(order='F')[:1 if truncate and slot==len(k) else 4].tofile(out/f'H.k{slot:06}.h')
        (out/'H.index').write_text('\n'.join(lines)+'\n')
        return 2.

    def test_resume_reuses_completed_point_and_publishes_only_missing_point(self):
        with patch.object(run_material,'run_native',side_effect=self.staged) as native:
            run_material.native_batch(self.root/'coeff',self.k[:1],self.out[:1],self.root,2)
            before=(self.out[0]/'COMPLETE.json').stat().st_mtime_ns
            run_material.native_batch(self.root/'coeff',self.k,self.out,self.root,2)
            np.testing.assert_array_equal(native.call_args.args[1],self.k[1:])
            self.assertEqual((self.out[0]/'COMPLETE.json').stat().st_mtime_ns,before)
            run_material.native_batch(self.root/'coeff',self.k,self.out,self.root,2)
            self.assertEqual(native.call_count,2)
        for index,out in enumerate(self.out):
            np.testing.assert_array_equal(run_material.load_h(out/'H.k000001.h',2),np.eye(2)*self.k[index,0])

    def test_missing_batch_keeps_background_cache_point_order(self):
        cache=[self.root/'cache0.bin',self.root/'cache1.bin']
        for p in cache:p.write_bytes(b'cache')
        def native(coeff,k,out,components,**kwargs):
            self.assertEqual(kwargs['background_files'],[cache[1]])
            return self.staged(coeff,k,out,components)
        with patch.object(run_material,'run_native',side_effect=lambda *a,**kw:self.staged(*a)):
            run_material.native_batch(self.root/'coeff',self.k[:1],self.out[:1],self.root,2,background_files=cache[:1])
        with patch.object(run_material,'run_native',side_effect=native) as called:
            run_material.native_batch(self.root/'coeff',self.k,self.out,self.root,2,background_files=cache)
            self.assertEqual(called.call_count,1)

    def test_rebuilds_incomplete_directory_and_rejects_truncated_new_batch(self):
        self.out[1].mkdir()
        (self.out[1]/'H.index').write_text('interrupted')
        with patch.object(run_material,'run_native',side_effect=lambda *a:self.staged(*a,truncate=True)):
            with self.assertRaisesRegex(ValueError,'Truncated'):
                run_material.native_batch(self.root/'coeff',self.k,self.out,self.root,2)
        self.assertFalse(any(p.exists() for p in self.out))
        self.assertEqual(len(list(self.root.glob('k_0001.incomplete-*'))),1)
        with patch.object(run_material,'run_native',side_effect=self.staged):
            run_material.native_batch(self.root/'coeff',self.k,self.out,self.root,2)
        self.assertTrue(all((p/'COMPLETE.json').exists() for p in self.out))

    def test_changed_coefficients_rejected_and_damaged_payload_rebuilt(self):
        with patch.object(run_material,'run_native',side_effect=self.staged) as native:
            run_material.native_batch(self.root/'coeff',self.k,self.out,self.root,2)
            (self.out[1]/'H.k000001.h').write_bytes(b'broken')
            run_material.native_batch(self.root/'coeff',self.k,self.out,self.root,2)
            self.assertEqual(len(native.call_args.args[1]),1)
            (self.root/'coeff').write_bytes(b'new spin')
            with self.assertRaisesRegex(ValueError,'inputs changed'):
                run_material.native_batch(self.root/'coeff',self.k,self.out,self.root,2)
