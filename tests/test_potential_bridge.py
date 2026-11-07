from pathlib import Path
import importlib.util
import shutil
import struct
import subprocess
import tempfile
import unittest
import numpy as np
from fortran_io import read_array,read_record,write_record
from parameters import ROOT,fortran_module,load


class PotentialBridgeTests(unittest.TestCase):
    def test_export_fitted_tensor(self):
        spec=importlib.util.spec_from_file_location('export_potential',ROOT/'4_potential/fit/export_potential.py')
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as directory:
            work=Path(directory);source=work/'fit';output=work/'tensor'
            ni=8
            linear=np.arange(ni*6,dtype='<f8').reshape((ni,6),order='F')
            modulation=linear+100
            with source.open('wb') as stream:
                for array in (b'CRI3_COM_FIT_CHANNEL_V1'.ljust(32),np.array([1,1,1],dtype='<i4'),
                              np.array([ni,6,12,2,4],dtype='<i4'),np.array([10,0],dtype='<i4'),
                              np.ones(5,dtype='<f8'),linear,modulation):write_record(stream,array)
            module.export(source,output,36)
            with output.open('rb') as stream:
                for _ in range(3):self.assertEqual(read_record(stream),b'')
                got=read_array(stream,'<f8',(ni,3,12,3))
            expected=np.zeros_like(got)
            for component in range(3):expected[:,component,:,component]=np.concatenate((linear,modulation),axis=1)
            np.testing.assert_array_equal(got,expected)
            with self.assertRaises(FileExistsError):module.export(source,output,36)

    @unittest.skipUnless(shutil.which('gfortran'),'gfortran is unavailable')
    def test_two_by_two_split_labels_and_record_count(self):
        with tempfile.TemporaryDirectory() as directory:
            work=Path(directory);p=load()
            p.update(split_nspinor=6,split_super_x=2,split_super_y=2,split_super_z=1)
            (work/'parameter_mod.f90').write_text(fortran_module(p))
            subprocess.run(['gfortran','-fcheck=all','-ffree-line-length-none',str(work/'parameter_mod.f90'),
                            str(ROOT/'4_potential/V_split.f90'),'-o','split'],cwd=work,check=True,capture_output=True)
            (work/'super_atom').write_text(''.join(f'{"Cr":14s}{0.:20.10f}{0.:20.10f}{0.:20.10f}\n' for _ in range(8)))
            (work/'polar2').write_text(''.join(f'{i} 90 0\n' for i in range(1,9)))
            field=np.arange(4*4*2*4,dtype='<f8').reshape((32,4),order='F')
            with (work/'V_HYC_dp').open('wb') as stream:
                write_record(stream,np.array([4,4,2],dtype='<i4'));write_record(stream,field)
            subprocess.run([str(work/'split')],cwd=work,check=True,capture_output=True)
            spins=np.loadtxt(work/'spin_info').reshape((4,6,3))
            np.testing.assert_array_equal(spins[:,:,0],[[2,1,3,5,6,4],[4,3,1,7,8,2],[6,5,7,1,2,8],[8,7,5,3,4,6]])
            for i in range(1,5):
                with (work/f'V{i}').open('rb') as stream:
                    self.assertEqual(int(read_array(stream,'<i4',(1,))[0]),27)
                    self.assertEqual(read_array(stream,'<f8',(27,4)).shape,(27,4))

    @unittest.skipUnless(shutil.which('gfortran'),'gfortran is unavailable')
    def test_fit_with_fewer_than_ten_input_files(self):
        with tempfile.TemporaryDirectory() as directory:
            work=Path(directory);p=load();p.update(scf_dir=str(work),n_file=2,n_kind=4,fit_max_iter=20)
            (work/'parameter_mod.f90').write_text(fortran_module(p))
            modules=['wannier_weight_mod','cri3_block_mask_mod','cri3_c3_basis_mod','fit_com_ortho_cold']
            subprocess.run(['gfortran','-O0','-fcheck=all','-ffree-line-length-none',str(work/'parameter_mod.f90')]+
                           [str(ROOT/'4_potential/fit'/(name+'.f90')) for name in modules]+['-o','fit'],cwd=work,check=True,capture_output=True)
            with (work/'V_QE').open('wb') as stream:write_record(stream,np.array([3,3,2],dtype='<i4'))
            with (work/'w_density').open('wb') as stream:
                write_record(stream,b'CRI3_WANNIER_WEIGHT_V1'.ljust(32))
                write_record(stream,np.array([3,3,2,48],dtype='<i4'))
                write_record(stream,struct.pack('<4id',1,44,1,1,.1));write_record(stream,np.ones(48,dtype='<f8'))
            rng=np.random.default_rng(44)
            for case in (1,2):
                folder=work/f'angle_{case}';folder.mkdir()
                spin=rng.normal(size=(4,6,3));np.savetxt(folder/'spin_info',spin.reshape((4,18)))
                for kind in range(4):
                    field=np.zeros((48,4),dtype='<f8');field[:,1:]=spin[kind].sum(axis=0)
                    with (folder/f'V{kind+1}').open('wb') as stream:
                        write_record(stream,np.array([48],dtype='<i4'));write_record(stream,field)
            result=subprocess.run([str(work/'fit')],cwd=work,check=True,capture_output=True,text=True)
            self.assertIn('fit_com_ortho done',result.stdout)
            self.assertTrue((work/'V_fit_com_ortho').is_file())


if __name__=='__main__':unittest.main(verbosity=2)
