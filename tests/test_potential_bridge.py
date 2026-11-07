"""Small explicit-channel split -> fit -> periodic-grid conversion fixtures."""
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
import numpy as np
from fortran_io import read_array,read_record,write_record
from parameters import ROOT,fortran_module,load


@unittest.skipUnless(shutil.which('gfortran'),'gfortran is unavailable')
class PotentialBridgeTests(unittest.TestCase):
    def test_split_fit_and_rvr_grid(self):
        with tempfile.TemporaryDirectory() as directory:
            w=Path(directory);p=load()
            (w/'parameter_mod.f90').write_text(fortran_module(p))
            for name,sources in [('split',['4_potential/V_split.f90']),
                                 ('fit',['4_potential/fit/fit_mask.f90','4_potential/fit/fit_new.f90']),
                                 ('grid',['4_potential/make_R_vec.f90'])]:
                subprocess.run(['gfortran','-fcheck=all','-ffree-line-length-none',str(w/'parameter_mod.f90')]+
                               [str(ROOT/s) for s in sources]+['-o',str(w/name)],cwd=w,check=True,capture_output=True)
            field=np.arange(4*4*2*4,dtype='<f8').reshape((32,4),order='F')/100
            with (w/'V_HYC_dp').open('wb') as f:
                write_record(f,np.array([4,4,2],dtype='<i4'));write_record(f,field)
            (w/'map').write_text('4 1 4\n3\n1\n4\n2\n')
            (w/'polar2').write_text('0.37 40 20\n0.4 80 100\n0.6 60 130\n0.8 120 270\n')
            subprocess.run([str(w/'split'),str(w/'map')],cwd=w,check=True,capture_output=True)
            angles=np.loadtxt(w/'polar2');theta,phi=np.deg2rad(angles[:,1:]).T
            expected=angles[:,0,None]*np.array([np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta)]).T
            np.testing.assert_allclose(np.loadtxt(w/'spin_info'),expected[[2,0,3,1]],atol=1e-16,rtol=0)
            for case,(x,y) in enumerate(((0,0),(0,1),(1,0),(1,1)),1):
                with (w/f'spin{case}/V_HYC_dp').open('rb') as f:
                    np.testing.assert_array_equal(read_array(f,'<i4',(3,)),[3,3,3])
                    got=read_array(f,'<f8',(27,4))
                indices=[(x*2+ix)%4+4*((y*2+iy)%4+4*(iz%2)) for iz in range(3) for iy in range(3) for ix in range(3)]
                np.testing.assert_array_equal(got,field[indices])
            matrices={}
            for stage in ('diagonal','full'):
                dest=w/(stage+'.model')
                (w/'fit.nml').write_text(f"&fit_control n_spinor=1,n_kind=4,diagonal_iters={int(stage=='diagonal')},full_iters={int(stage=='full')},lambda_DM=0.01,fields_dir='{w}',spin_file='{w/'spin_info'}',output_file='{dest}' /\n")
                subprocess.run([str(w/'fit'),str(w/'fit.nml')],cwd=w,check=True,capture_output=True)
                with dest.open('rb') as f:
                    for _ in range(3):read_record(f)
                    matrices[stage]=read_array(f,'<f8',(27,3,1,3))
            self.assertEqual(np.max(abs(matrices['diagonal'][:,0,0,1])),0)
            self.assertGreater(np.max(abs(matrices['full'][:,0,0,1])),1e-6)
            with (w/'base').open('wb') as f:write_record(f,np.array([2,2,2],dtype='<i4'))
            subprocess.run([str(w/'grid'),str(w/'base'),str(w/'full.model'),str(w/'cropped'),str(w/'R_vec')],cwd=w,check=True,capture_output=True)
            with (w/'cropped').open('rb') as f:
                records=[read_record(f) for _ in range(5)]
                self.assertEqual(f.read(),b'')
            self.assertEqual([len(x) for x in records],[8*4*8,8*3*8,8*8*8,8*9*8,8*8])
            got=np.frombuffer(records[3],dtype='<f8').reshape((8,3,1,3),order='F')
            take=[x+3*(y+3*z) for z in range(2) for y in range(2) for x in range(2)]
            np.testing.assert_array_equal(got,matrices['full'][take])
            # Extra tensor data must fail rather than be silently truncated by sequential IO.
            wrong=w/'wrong.model';wrong.write_bytes((w/'full.model').read_bytes()+b'bad')
            failed=subprocess.run([str(w/'grid'),str(w/'base'),str(wrong),str(w/'bad-crop'),str(w/'bad-vec')],cwd=w,capture_output=True,text=True)
            self.assertNotEqual(failed.returncode,0)
            self.assertIn('format mismatch',failed.stderr)
            self.assertFalse((w/'bad-crop').exists())
            (w/'bad-restart.nml').write_text(f"&fit_control n_spinor=1,n_kind=4,diagonal_iters=0,full_iters=1,lambda_DM=0.01,fields_dir='{w}',spin_file='{w/'spin_info'}',restart_file='{wrong}',output_file='{w/'bad-fit'}' /\n")
            failed=subprocess.run([str(w/'fit'),str(w/'bad-restart.nml')],cwd=w,capture_output=True,text=True)
            self.assertNotEqual(failed.returncode,0)
            self.assertIn('format mismatch',failed.stderr)
            self.assertFalse((w/'bad-fit').exists())

    def test_missing_map_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            w=Path(directory);(w/'parameter_mod.f90').write_text(fortran_module(load()))
            subprocess.run(['gfortran',str(w/'parameter_mod.f90'),str(ROOT/'4_potential/V_split.f90'),'-o',str(w/'split')],cwd=w,check=True)
            r=subprocess.run([str(w/'split')],cwd=w,capture_output=True,text=True)
            self.assertNotEqual(r.returncode,0)
            self.assertIn('explicit-channel-map',r.stderr)


if __name__=='__main__':unittest.main(verbosity=2)
