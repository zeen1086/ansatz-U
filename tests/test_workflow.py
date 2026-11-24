from pathlib import Path
import contextlib
import importlib.util
import io
import json
import os
import shutil
import struct
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from tools.build import build, verify
from tools.configure import configure
from tools.fortran_io import read_array, read_record, write_record
from tools.parameters import ROOT, load


def module(name, relative):
    spec = importlib.util.spec_from_file_location(name, ROOT/relative)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result

prepare = module("prepare_wannier", "5_wannier/prepare.py")


class WorkflowTests(unittest.TestCase):
    def test_qe_installer_rejects_wrong_file_size(self):
        installer=module("qe_install", "0_QE/install.py")
        with tempfile.TemporaryDirectory() as directory:
            base=Path(directory); overlay=base/"overlay"; target=base/"qe"
            (overlay/"common").mkdir(parents=True)
            (overlay/"magnetic").mkdir()
            target.mkdir()
            (overlay/"common/source.f90").write_text("modified source\n")
            (overlay/"upstream.json").write_text(json.dumps({"files":{"source.f90":5}}))
            (target/"source.f90").write_text("user changes\n")
            with patch.object(installer,"HERE",overlay):
                with self.assertRaises(ValueError): installer.install(target,"magnetic")
                self.assertEqual((target/"source.f90").read_text(),"user changes\n")
                (target/"source.f90").write_text("base\n")
                installer.install(target,"magnetic")
                self.assertEqual((target/"source.f90").read_text(),"modified source\n")

    def test_record_validation_and_split_records(self):
        values=np.arange(12,dtype="<f8").reshape((3,4),order="F")
        stream=io.BytesIO()
        write_record(stream,values)
        stream.seek(0)
        np.testing.assert_array_equal(read_array(stream,"<f8",(3,4)),values)
        split=struct.pack("<i",-3)+b"abc"+struct.pack("<i",3)+struct.pack("<i",2)+b"de"+struct.pack("<i",-2)
        self.assertEqual(read_record(io.BytesIO(split)),b"abcde")
        for bad in (split[:-1], struct.pack("<i",3)+b"abc"+struct.pack("<i",2)):
            with self.assertRaises(ValueError): read_record(io.BytesIO(bad))

    def test_wannier_text_conversion(self):
        with tempfile.TemporaryDirectory() as directory:
            base=Path(directory)
            p=load(); p.update(nwan=2,num_bands=3,nk1=2,nk2=1,nk3=1,prefix="toy")
            def write_u(name, array):
                with (base/name).open("w") as stream:
                    stream.write(f"fixture\n2 {array.shape[1]} {array.shape[0]}\n")
                    for k in range(2):
                        stream.write(f"\n{k/2} 0 0\n")
                        for z in array[:,:,k].flatten(order="F"):
                            stream.write(f"{z.real} {z.imag}\n")
            u=np.stack([np.eye(2),np.array([[0,1j],[1,0]])],axis=2).astype(complex)
            dis=np.arange(12).reshape((3,2,2),order="F")+0.1j
            write_u("toy_u.mat",u); write_u("toy_u_dis.mat",dis)
            (base/"toy.win").write_text("use_ws_distance = false\n# dis_win_min = -10\n")
            (base/"toy_hr.dat").write_text("fixture\n2\n1\n2\n0 0 0 1 1 1 0\n0 0 0 2 1 2 -1\n0 0 0 1 2 2 1\n0 0 0 2 2 3 0\n")
            output=base/"converted"
            prepare.prepare(base,output,p)
            with (output/"u_dis_mat_HYC").open("rb") as stream:
                np.testing.assert_array_equal(read_array(stream,"<i4",(3,)),[3,2,2])
                np.testing.assert_array_equal(read_array(stream,"<c16",(3,2,2)),dis)
            with (output/"hr_HYC_1").open("rb") as stream:
                np.testing.assert_array_equal(read_array(stream,"<c16",(2,2)),[[1,2+1j],[2-1j,3]])
            with self.assertRaises(FileExistsError): prepare.prepare(base,output,p)
            (base/"toy.win").write_text("dis_win_max=4\n")
            with self.assertRaises(ValueError): prepare.prepare(base,base/"invalid",p)

    def test_isolated_generation(self):
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory)
            configure(output)
            self.assertTrue((output/"1_scf/scf.in").is_file())
            self.assertTrue((output/".build/parameter_mod.f90").is_file())

    @unittest.skipUnless(shutil.which("gfortran"),"gfortran is unavailable")
    def test_build_fingerprints_and_atomic_failure(self):
        (ROOT/".build").mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(dir=ROOT/".build") as directory:
            base=Path(directory); source=base/"probe.f90"; config=base/"parameter"
            config.write_text((ROOT/"parameter").read_text())
            source.write_text("program probe\nuse parameter_mod\ncharacter(len=256)::argument\n"
                              "call get_command_argument(1,argument)\nif(len_trim(argument)>0)then\n"
                              "open(10,file='run-output.txt',status='replace')\n"
                              "write(10,'(A)')trim(argument)\nclose(10)\nendif\n"
                              "print *,p_nwan\nend program\n")
            with patch.dict(os.environ,{"ANSATZ_PARAMETER":str(config),"ANSATZ_BUILD_DIR":str(base/"build")}):
                build([source],"gfortran",["-O0"],[])
                binary=verify(source); stamp=binary.stat().st_mtime_ns
                build([source],"gfortran",["-O0"],[])
                self.assertEqual(binary.stat().st_mtime_ns,stamp)
                self.assertEqual(subprocess.check_output([str(binary)],text=True).strip(),"30")
                with tempfile.TemporaryDirectory() as run_directory:
                    subprocess.run([sys.executable,str(ROOT/'tools/configure.py'),'--check'],
                                   cwd=run_directory,check=True,capture_output=True)
                    subprocess.run([sys.executable,str(ROOT/'tools/build.py'),str(source),
                                    '--fc','gfortran','--flags=-O0','--libs='],
                                   cwd=run_directory,check=True,capture_output=True)
                    command=[sys.executable,str(ROOT/'tools/run.py'),str(source.relative_to(ROOT)),
                             'argument with spaces']
                    self.assertEqual(subprocess.check_output(command,cwd=run_directory,text=True).strip(),"30")
                    self.assertEqual((Path(run_directory)/'run-output.txt').read_text(),'argument with spaces\n')
                    config.write_text(config.read_text().replace("nr1 = 64","nr1 = 32"))
                    rejected=subprocess.run(command,cwd=run_directory,text=True,capture_output=True)
                    self.assertNotEqual(rejected.returncode,0)
                    self.assertIn('Rebuild',rejected.stderr)
                with self.assertRaises(RuntimeError): verify(source)
                build([source],"gfortran",["-O0"],[]); verify(source)
                original=binary.read_bytes()
                source.write_text("broken Fortran\n")
                with self.assertRaises(subprocess.CalledProcessError):
                    build([source],"gfortran",["-O0"],[])
                self.assertEqual(binary.read_bytes(),original)
                with self.assertRaises(RuntimeError): verify(source)

    @unittest.skipUnless(shutil.which("gfortran"),"gfortran is unavailable")
    def test_qe_wavefunction_export_padding(self):
        source=(ROOT/"0_QE/common/PP/src/wfck2r.f90").read_text()
        helper=source[source.index("  subroutine write_hyc_wavefunction"):source.index("end PROGRAM wfck2r")]
        with tempfile.TemporaryDirectory() as directory:
            work=Path(directory)
            code="""program fixture
integer,parameter::DP=kind(0d0)
complex(DP)::values(9,2)
integer::j
values=reshape([(cmplx(j,-j,kind=DP),j=1,18)],[9,2])
call write_hyc_wavefunction('wfn',values,2,2,1,3,3,2)
contains
"""+helper+"end program fixture\n"
            (work/"fixture.f90").write_text(code)
            subprocess.run(["gfortran","-fcheck=all","fixture.f90","-o","fixture"],cwd=work,check=True,capture_output=True)
            subprocess.run([str(work/"fixture")],cwd=work,check=True)
            with (work/"wfn").open("rb") as stream: wave=read_array(stream,"<c16",(4,2))
            np.testing.assert_array_equal(wave[:,0],np.array([1,2,4,5])*(1-1j))
            np.testing.assert_array_equal(wave[:,1],np.array([10,11,13,14])*(1-1j))


if __name__ == "__main__": unittest.main(verbosity=2)
