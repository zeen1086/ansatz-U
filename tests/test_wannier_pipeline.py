from pathlib import Path
import os
import shutil
import shlex
import subprocess
import tempfile
import unittest
import numpy as np

from tools.build import DEFAULT_LIBS
from tools.parameters import ROOT, load, fortran_module
from tools.fortran_io import write_record, read_array


class WannierPipelineTests(unittest.TestCase):
    @unittest.skipUnless((shutil.which("mpifort") or shutil.which("mpiifort")) and shutil.which("mpirun"),
                         "MPI compiler and mpirun are unavailable")
    def test_mpi_wannier_normalization(self):
        compiler=shutil.which("mpiifort") or shutil.which("mpifort")
        flags=(["-O0","-check","bounds","-mkl"] if "mpiifort" in compiler else
               ["-O0","-fcheck=all","-ffree-line-length-none","-fallow-argument-mismatch"])
        r=ROOT
        with tempfile.TemporaryDirectory() as directory:
            w=Path(directory)
            p=load();p.update(nwan=1,wan_min=1,wan_max=1,nspin=1,nk1=1,nk2=1,nk3=1,
                             convert_size_x=0,convert_size_y=0,convert_size_z=0,first_band=1,last_band=2,
                             num_bands=2,nscf_dir=str(w/'wannier'))
            wan=w/'wannier';wan.mkdir();(wan/'WFN').mkdir();(wan/'disentangled').mkdir();(wan/'make_wan').mkdir()
            (w/'parameter_mod.f90').write_text(fortran_module(p))
            for name, matrix in [('u_mat_HYC',np.ones((1,1,1),dtype='<c16')),
                                 ('u_dis_mat_HYC',np.array([1,1j],dtype='<c16').reshape((2,1,1))/np.sqrt(2))]:
             with (wan/name).open('wb') as stream:
              write_record(stream,np.array(matrix.shape,dtype='<i4'));write_record(stream,matrix)
            with (wan/'window_TF_HYC').open('wb') as stream:write_record(stream,np.ones((2,1),dtype='<i4'))
            (wan/'k_list').write_text(f'{0:12.8f}{0:12.8f}{0:12.8f}\n')
            with (wan/'V_HYC_dp').open('wb') as stream:write_record(stream,np.array([2,1,1],dtype='<i4'))
            for band,wave in enumerate(([1,1],[1,-1]),start=1):
             with (wan/'WFN'/f'wfn_k1_b{band}_vec_dp').open('wb') as stream:write_record(stream,np.array(wave,dtype='<c16'))
            env=dict(os.environ,MKL_NUM_THREADS='1',OMP_NUM_THREADS='1',I_MPI_PIN='0')
            for name,sources in [('dis',['5_wannier/disentangled/make_dis_mpi.f90']),('wan',['5_wannier/make_wan/make_wan_fun_mpi.f90'])]:
             cmd=[compiler]+flags+[str(w/'parameter_mod.f90')]+[str(r/s) for s in sources]+['-o',str(w/name)]+([] if 'mpiifort' in compiler else shlex.split(os.environ.get('LIBS', DEFAULT_LIBS)))
             subprocess.run(cmd,cwd=w,check=True)
            subprocess.run(['mpirun','-np','2',str(w/'dis')],cwd=wan/'disentangled',env=env,check=True,timeout=30)
            subprocess.run(['mpirun','-np','2',str(w/'wan'),'1'],cwd=wan/'make_wan',env=env,check=True,timeout=30)
            with (wan/'make_wan/wan_1').open('rb') as stream:a=read_array(stream,'<c16',(1,2,1))
            expected=np.array([1+1j,1-1j]).reshape((1,2,1))/2
            np.testing.assert_allclose(a,expected,atol=1e-14,rtol=0)
            self.assertAlmostEqual(float(np.sum(abs(a)**2)),1.0,places=13)


if __name__=='__main__':unittest.main(verbosity=2)
