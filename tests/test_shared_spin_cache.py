"""Shared cache publication, input identity and reuse on a small real grid."""
from pathlib import Path
import sys
import json
import tempfile
import unittest
from unittest.mock import patch
from concurrent.futures import ThreadPoolExecutor
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'7_spin'))
import common
import complete_energy as driver


class SharedCacheTests(unittest.TestCase):
    def test_cache_reuse_checks_size_without_hashing(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'cache'
            common.shared_bundle(path,{},lambda out:np.save(out/'a.npy',np.ones(3)))
            record=json.loads((path/'CACHE.json').read_text())
            self.assertEqual(set(record['files']['a.npy']),{'bytes'})
            np.save(path/'a.npy',np.zeros(3))
            self.assertEqual(common.shared_bundle(path,{},lambda out:self.fail('rebuilt')),path)
            (path/'a.npy').write_bytes(b'truncated')
            with self.assertRaisesRegex(ValueError,'corrupt payload'):
                common.shared_bundle(path,{},lambda out:self.fail('rebuilt'))

    def test_atomic_publish_concurrent_reuse_and_rejection(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);path=root/'entry';calls=[]
            def build(out):
                calls.append(1);np.save(out/'a.npy',np.eye(3))
            with ThreadPoolExecutor(2) as pool:
                results=list(pool.map(lambda _:common.shared_bundle(path,{'source':'a'},build),range(2)))
            self.assertEqual(results,[path,path]);self.assertEqual(len(calls),1)
            stamp=(path/'a.npy').stat().st_mtime_ns
            common.shared_bundle(path,{'source':'a'},build)
            self.assertEqual((path/'a.npy').stat().st_mtime_ns,stamp)
            with self.assertRaisesRegex(ValueError,'identity mismatch'):
                common.shared_bundle(path,{'source':'b'},build)
            data=(path/'a.npy').read_bytes();(path/'a.npy').write_bytes(data[:-1])
            with self.assertRaisesRegex(ValueError,'corrupt payload'):
                common.shared_bundle(path,{'source':'a'},build)
            (root/'incomplete').mkdir()
            with self.assertRaisesRegex(ValueError,'Invalid shared cache'):
                common.shared_bundle(root/'incomplete',{},build)
            def fail(out):
                np.save(out/'a.npy',np.eye(3));raise RuntimeError('interrupted')
            with self.assertRaises(RuntimeError):common.shared_bundle(root/'failed',{},fail)
            self.assertFalse((root/'failed').exists())
            common.shared_bundle(root/'failed',{},build)
            def nonfinite(out):np.save(out/'a.npy',np.array([np.nan]))
            with self.assertRaisesRegex(ValueError,'Nonfinite'):
                common.shared_bundle(root/'nan',{},nonfinite)
            self.assertFalse((root/'nan').exists())

    def test_real_wave_metric_operators_reused(self):
        rng=np.random.default_rng(16)
        wave=np.linalg.qr(rng.normal(size=(32,30))+1j*rng.normal(size=(32,30)))[0]
        slab=np.broadcast_to(wave,(8,8,32,30)).copy()
        st=(np.zeros((1,3)),np.eye(30)[None],np.ones((1,5,30),complex)*.01)
        k=np.array([.13,.21,0.])
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'RVR').mkdir()
            (root/'RVR/PROJECTION.json').write_text('{"basis_id":"'+('b'*64)+'"}')
            with patch.object(driver,'SHARED',root),patch.object(driver,'WORK',root), \
                 patch.object(driver,'NG',32),patch.object(driver,'points',return_value=[(0,k)]):
                paths,ops=driver.wave_cache(1,0,{'basis':'fixture'},st,slab)
                np.testing.assert_allclose(np.load(paths[0]),wave[None],atol=1e-14)
                np.testing.assert_allclose(np.load(paths[0].parent/'transform.npy'),np.eye(60),atol=1e-13)
                h,a=driver.load_operators(ops[0])
                want_h,want_a=driver.operators(1,k,*st)
                np.testing.assert_array_equal(h,want_h);np.testing.assert_array_equal(a,want_a)
                before=paths[0].stat().st_mtime_ns
                with patch.object(driver,'operators',side_effect=AssertionError('recomputed')), \
                     patch.object(driver,'metric',side_effect=AssertionError('recomputed')):
                    again,_=driver.wave_cache(1,0,{'basis':'fixture'},st,slab)
                self.assertEqual(again,paths);self.assertEqual(paths[0].stat().st_mtime_ns,before)
                with self.assertRaisesRegex(ValueError,'identity mismatch'):
                    driver.wave_cache(1,0,{'basis':'changed'},st,slab)

    def test_identity_tracks_input_metadata(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            for name in ('basis/bloch.npy','basis/pure_H0.au','RVR/PROJECTION.json',
                         'RVR/B_general','V_func_tmp/V_FUNC.json','refs/data.npz',
                         'reference/results/FIELDS.npz','reference/inputs/reference_actual/DENSITY.bin'):
                p=root/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(b'input')
            with patch.object(driver,'WORK',root),patch.object(driver,'INTERFACE',root), \
                 patch.object(driver,'REF',root/'refs'),patch.object(driver,'MIN',root/'reference'):
                identity=driver.cache_identity()
                for name in ('basis/bloch.npy','basis/pure_H0.au','RVR/PROJECTION.json',
                             'V_func_tmp/V_FUNC.json','reference/results/FIELDS.npz'):
                    p=root/name;p.write_bytes(b'other')
                    self.assertNotEqual(driver.cache_identity(),identity)
                    p.write_bytes(b'input')
                current=driver.cache_identity()
                self.assertEqual(driver.cache_identity(),current)

    def test_binding_rejects_changed_inputs_and_registers_existing_data(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'input';source.write_bytes(b'first')
            identity={'input':common.input_digest(source)}
            with patch.object(driver,'SHARED',root):
                driver.bind_cache(identity);driver.bind_cache(identity)
                import os
                stamp=source.stat()
                source.write_bytes(b'other')
                os.utime(source,ns=(stamp.st_atime_ns,stamp.st_mtime_ns))
                with self.assertRaisesRegex(ValueError,'inputs changed'):
                    driver.bind_cache({'input':common.input_digest(source)})
                (root/'INPUTS.json').unlink();(root/'reference').mkdir()
                driver.bind_cache(identity)
                self.assertTrue((root/'INPUTS.json').is_file())


if __name__=='__main__':unittest.main()
