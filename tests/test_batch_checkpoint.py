"""Check that a restart reuses completed points and never publishes a partial batch."""
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import math
import unittest
from unittest.mock import patch
import numpy as np

import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'7_spin'))
import run_material
import common


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


class BandNumericsTest(unittest.TestCase):
    def test_eigensystem_reuses_valid_cache_and_rebuilds_invalid_cache(self):
        with TemporaryDirectory() as tmp:
            folder=Path(tmp);(folder/'value').mkdir()
            h=np.array([[.2,.1+.3j],[.1-.3j,-.4]],dtype='<c16')
            h.ravel(order='F').tofile(folder/'value/H.k000001.h')
            common.ensure_eigensystem(folder,2)
            cache=folder/'eigensystem.npz'
            original=cache.read_bytes()
            with patch.object(common,'load_h',side_effect=AssertionError('valid cache reread Hamiltonian')):
                common.ensure_eigensystem(folder,2)
            self.assertEqual(cache.read_bytes(),original)
            for damage in ('truncated','wrong_shape'):
                with self.subTest(damage=damage):
                    if damage=='truncated':cache.write_bytes(b'partial')
                    else:np.savez(cache,e=np.zeros(1),v=np.ones((1,1)))
                    with patch.object(common.np.linalg,'eigh',wraps=np.linalg.eigh) as solve:
                        common.ensure_eigensystem(folder,2)
                        self.assertEqual(solve.call_count,1)
                    with np.load(cache) as data:
                        np.testing.assert_allclose(data['e'],np.linalg.eigvalsh(h),atol=1e-15)
                        np.testing.assert_allclose((data['v']*data['e'])@data['v'].conj().T,h,atol=1e-15)

    def test_band_density_matches_independent_state_sum(self):
        rng=np.random.default_rng(35)
        v,_=np.linalg.qr(rng.normal(size=(4,4))+1j*rng.normal(size=(4,4)))
        mu=.23;e=mu+common.SIGMA*np.array([-3.,-.2,.6,2.])
        occupations=np.array([.5*math.erfc((ev-mu)/common.SIGMA) for ev in e])
        expected=sum(f*np.outer(v[:,i],v[:,i].conj()) for i,f in enumerate(occupations))
        free=sum(ev*f-common.SIGMA*math.exp(-((ev-mu)/common.SIGMA)**2)/(2*math.sqrt(math.pi))
                 for ev,f in zip(e,occupations))
        density,stats=common.band_density(e,v,mu)
        np.testing.assert_allclose(density,expected,rtol=0,atol=1e-15)
        self.assertAlmostEqual(stats['count'],sum(occupations),places=14)
        self.assertAlmostEqual(stats['F'],free,places=14)
