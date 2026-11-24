"""Small explicit-channel potential splitting fixtures."""
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
import numpy as np
from tools.fortran_io import read_array,write_record
from tools.parameters import ROOT,fortran_module,load


@unittest.skipUnless(shutil.which('gfortran'),'gfortran is unavailable')
class PotentialBridgeTests(unittest.TestCase):
    def test_split_mapping_and_periodic_grid(self):
        with tempfile.TemporaryDirectory() as directory:
            w=Path(directory);p=load()
            (w/'parameter_mod.f90').write_text(fortran_module(p))
            subprocess.run(['gfortran','-fcheck=all','-ffree-line-length-none',str(w/'parameter_mod.f90'),
                            str(ROOT/'4_potential/V_split.f90'),'-o',str(w/'split')],
                           cwd=w,check=True,capture_output=True)
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

    def test_missing_map_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            w=Path(directory);(w/'parameter_mod.f90').write_text(fortran_module(load()))
            subprocess.run(['gfortran',str(w/'parameter_mod.f90'),str(ROOT/'4_potential/V_split.f90'),'-o',str(w/'split')],cwd=w,check=True)
            r=subprocess.run([str(w/'split')],cwd=w,capture_output=True,text=True)
            self.assertNotEqual(r.returncode,0)
            self.assertIn('explicit-channel-map',r.stderr)


if __name__=='__main__':unittest.main(verbosity=2)
