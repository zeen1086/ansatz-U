from pathlib import Path
import contextlib
import importlib.util
import hashlib
import io
import json
import os
import shutil
import shlex
import struct
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from build import DEFAULT_LIBS, build, locations, verify
from configure import configure
from fortran_io import read_array, read_record, write_record
from parameters import ROOT, fortran_module, load


def module(name, relative):
    spec = importlib.util.spec_from_file_location(name, ROOT/relative)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result

prepare = module("prepare_wannier", "5_wannier/prepare.py")
rvr = module("prepare_rvr", "5_wannier/make_wan/prepare_rvr.py")


class WorkflowTests(unittest.TestCase):
    def test_qe_installer_rejects_modified_tree(self):
        installer=module("qe_install", "0_QE/install.py")
        with tempfile.TemporaryDirectory() as directory:
            base=Path(directory); overlay=base/"overlay"; target=base/"qe"
            (overlay/"common").mkdir(parents=True)
            (overlay/"magnetic").mkdir()
            target.mkdir()
            (overlay/"common/source.f90").write_text("modified source\n")
            (overlay/"upstream.json").write_text(json.dumps({"files":{"source.f90":hashlib.sha256(b"base\n").hexdigest()}}))
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

    def test_wannier_to_rvr_grid_and_spin_order(self):
        with tempfile.TemporaryDirectory() as directory:
            base=Path(directory); raw=base/"raw"; raw.mkdir()
            p={"nscf_dir":str(base),"func_dir":str(base),"convert_size_x":0,"convert_size_y":0,
               "convert_size_z":0,"nspin":2,"nwan":1}
            with (base/"V_HYC_dp").open("wb") as stream:
                write_record(stream,np.array([2,2,1],dtype="<i4"))
            vectors=np.array([[1,1,0],[0,0,0],[1,0,0],[0,1,0]],dtype="<i4")
            with (base/"R_vec").open("wb") as stream:
                write_record(stream,np.array([4],dtype="<i4")); write_record(stream,vectors)
            wave=(np.arange(8)+1j*np.arange(8)[::-1]).reshape((1,4,2),order="F").astype("<c16")
            with (raw/"wan_1").open("wb") as stream: write_record(stream,wave)
            rvr.prepare(raw,base/"out",p)
            got=np.fromfile(base/"out/wan_data/wan_1",dtype="<c16").reshape((1,4,2),order="F")
            np.testing.assert_array_equal(got,wave[:,[3,0,1,2],:])
            self.assertAlmostEqual(float(np.sum(abs(got)**2)),float(np.sum(abs(wave)**2)))
            self.assertEqual(np.loadtxt(base/"out/cal_TF_1"),1)

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
            source.write_text("program probe\nuse parameter_mod\nprint *,p_nwan\nend program\n")
            with patch.dict(os.environ,{"ANSATZ_PARAMETER":str(config),"ANSATZ_BUILD_DIR":str(base/"build")}):
                build([source],"gfortran",["-O0"],[])
                binary=verify(source); stamp=binary.stat().st_mtime_ns
                build([source],"gfortran",["-O0"],[])
                self.assertEqual(binary.stat().st_mtime_ns,stamp)
                self.assertEqual(subprocess.check_output([str(binary)],text=True).strip(),"44")
                config.write_text(config.read_text().replace("lambda_high = 0.0","lambda_high = 0.5"))
                with self.assertRaises(RuntimeError): verify(source)
                build([source],"gfortran",["-O0"],[]); verify(source)
                original=binary.read_bytes()
                source.write_text("broken Fortran\n")
                with self.assertRaises(subprocess.CalledProcessError):
                    build([source],"gfortran",["-O0"],[])
                self.assertEqual(binary.read_bytes(),original)
                with self.assertRaises(RuntimeError): verify(source)

    @unittest.skipUnless(shutil.which("gfortran"),"gfortran is unavailable")
    def test_fortran_spin_and_rvr_known_answers(self):
        with tempfile.TemporaryDirectory() as directory:
            work=Path(directory)
            (work/"parameter_mod.f90").write_text(fortran_module(load()))
            for test, dependency in (("test_spin_channels","7_spin/spin_channels_mod.f90"),
                                     ("test_rvr_kernel","6_RVR/speedup/rvr_kernel_mod.f90")):
                command=["gfortran","-O2","-fcheck=all",str(work/"parameter_mod.f90"),str(ROOT/dependency),
                         str(ROOT/f"tests/{test}.f90"),"-o",str(work/test)] + shlex.split(os.environ.get("LIBS", DEFAULT_LIBS))
                subprocess.run(command,cwd=work,check=True,capture_output=True)
                run=subprocess.run([str(work/test)],cwd=work,check=True,capture_output=True,text=True)
                self.assertIn("PASS",run.stdout)

    @unittest.skipUnless(shutil.which("gfortran"),"gfortran is unavailable")
    def test_hamiltonian_trace_and_supercell_fold(self):
        with tempfile.TemporaryDirectory() as directory:
            work=Path(directory)
            primitive=work/"primitive"; primitive.mkdir()
            folded=work/"super"; folded.mkdir()
            vectors=np.array([[-1,0,0],[0,0,0],[1,0,0]])
            np.savetxt(primitive/"ir_vec_HYC",vectors,fmt="%8d%8d%8d")
            np.savetxt(primitive/"H_deg_HYC",[[1,2],[2,1],[3,2]],fmt="%8d%8d")
            matrix=np.array([[1,2+1j],[2-1j,3]],dtype="<c16")
            for row in range(3):
                with (primitive/f"hr_HYC_{row+1}").open("wb") as stream:
                    write_record(stream,matrix*(row+1))
            p=load(); p.update(nwan=2,wan_max=2,nscf_dir=str(primitive),super_h_dir=str(folded),
                               period_x=2,period_y=1,period_z=1)
            (work/"parameter_mod.f90").write_text(fortran_module(p))
            subprocess.run(["gfortran","-fcheck=all",str(work/"parameter_mod.f90"),str(ROOT/"6_RVR/rearrange_HR.f90"),
                            "-o",str(work/"fold")],cwd=work,check=True,capture_output=True)
            subprocess.run([str(work/"fold")],cwd=work,check=True,capture_output=True)
            output_vectors=np.loadtxt(folded/"ir_vec_HYC",dtype=int).reshape((-1,3))
            for block, vector in enumerate(output_vectors):
                expected=np.zeros((4,4),dtype=complex)
                for cell in range(2):
                    for row, displacement in enumerate((-1,0,1)):
                        target=cell+displacement
                        if target//2 != vector[0]: continue
                        expected[2*cell:2*cell+2,2*(target%2):2*(target%2)+2]=matrix*(row+1)/[2,1,2][row]
                with (folded/f"hr_HYC_{block+1}").open("rb") as stream:
                    np.testing.assert_array_equal(read_array(stream,"<c16",(4,4)),expected)
            # Compare the actual normal and trace entry points with non-default controls.
            np.savetxt(folded/"ir_vec_HYC",[[0,0,0]],fmt="%8d%8d%8d")
            (work/"spin.in").write_text("30 80\n20 110\n")
            p.update(period_x=1,rvr_dir=str(work),nx=1,ny=1,nz=1,convert_size_x=0,
                     convert_size_y=0,convert_size_z=0,trace_i=1,trace_j=2,
                     channel_mask=[0,1,-1,1,0,1],channel_source=[6,4,2,1,3,3],lambda_high=0.7)
            for nf in (18,36):
                p["n_func"]=nf
                (work/"parameter_mod.f90").write_text(fortran_module(p))
                with (work/"B_general").open("wb") as stream:
                    write_record(stream,(np.arange(nf*16).reshape((nf,4,4),order="F")*0.01j+0.02).astype("<c16"))
                for source in ("cal_Hsk","cal_Hsk_trace"):
                    subprocess.run(["gfortran","-fcheck=all",str(work/"parameter_mod.f90"),str(ROOT/"7_spin/spin_channels_mod.f90"),
                                    str(ROOT/f"7_spin/{source}.f90"),"-o",str(work/source)],cwd=work,check=True,capture_output=True)
                    subprocess.run([str(work/source)],cwd=work,check=True,capture_output=True)
                with (work/"data/Hsk_1").open("rb") as stream: normal=read_array(stream,"<c16",(4,4))
                with (work/"data_trace/Hsk_1").open("rb") as stream: trace=read_array(stream,"<c16",(4,4))
                np.testing.assert_array_equal(normal,trace)

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
