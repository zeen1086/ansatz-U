"""Small offline checks: shared settings, safe generation and Fortran record IO."""
from pathlib import Path
import ast
import shutil
import struct
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from configure import configure, render
from parameters import ROOT, fortran_module, load


def write_records(path: Path, *records: bytes) -> None:
    with path.open("wb") as stream:
        for record in records:
            marker = struct.pack("<i", len(record))
            stream.write(marker + record + marker)


def read_records(path: Path) -> list[bytes]:
    data, records, offset = path.read_bytes(), [], 0
    while offset < len(data):
        size = struct.unpack_from("<i", data, offset)[0]
        end = offset + 4 + size
        if size < 0 or struct.unpack_from("<i", data, end)[0] != size:
            raise ValueError("invalid record markers")
        records.append(data[offset+4:end])
        offset = end + 4
    return records


class ConfigurationTests(unittest.TestCase):
    def altered(self, before: str, after: str) -> dict:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "parameter"
            path.write_text((ROOT / "parameter").read_text().replace(before, after))
            return load(path)

    def test_shared_grid_and_lattice_reach_all_inputs(self):
        p = self.altered("nk1 = 12", "nk1 = 6")
        output = render(p)
        self.assertIn("nk1 = 6", fortran_module(p))
        self.assertIn("6 12 1 0 0 0", output[ROOT / "1_scf/scf.in"])
        self.assertIn("mp_grid = 6 12 1", output[ROOT / "5_wannier/CrI3.win"])
        self.assertIn("last_k=72", output[ROOT / "5_wannier/wfck2r.in"])
        self.assertIn("write_v_qe=.true.", output[ROOT / "1_scf/print_V.in"])
        no_raw = render(self.altered("write_v_qe = true", "write_v_qe = false"))
        self.assertIn("write_v_qe=.false.", no_raw[ROOT / "1_scf/print_V.in"])
        krows = output[ROOT / "5_wannier/nscf.in"].split("K_POINTS crystal\n")[1].splitlines()
        self.assertEqual(int(krows[0]), 72)
        self.assertAlmostEqual(sum(float(row.split()[3]) for row in krows[1:]), 1)
        self.assertIn("nline=  3", output[ROOT / "7_spin/plt_path.in"])
        self.assertIn("npt= 20", output[ROOT / "7_spin/plt_path.in"])
        spinor = render(self.altered("nspin = 1", "nspin = 2"))
        self.assertIn("spinors = true", spinor[ROOT / "5_wannier/CrI3.win"])
        self.assertIn("noncolin=.true.", spinor[ROOT / "5_wannier/nscf.in"])
        changed = render(self.altered("a0 = 13.17139000", "a0 = 15.0"))
        self.assertIn("celldm(1)=15.0", changed[ROOT / "2_band/band.in"])
        self.assertNotEqual(changed[ROOT / "5_wannier/CrI3.win"], output[ROOT / "5_wannier/CrI3.win"])

    def test_rejects_inconsistent_parameters(self):
        for before, after in [("nwan = 44", "nwan = 20"), ("period_x = 5", "period_x = 7"),
                              ("nk3 = 1", "nk3 = 0"), ("n_func = 36", "n_func = 42"),
                              ("selected_eta = 1", "selected_eta = 2"),
                              ("eta = 0.002", "eta = nan"),
                              ('export_grid = "smooth"', 'export_grid = "invalid"'),
                              ("write_v_qe = true", 'write_v_qe = "true"'),
                              ("num_bands = 80", "num_bands = 79")]:
            with self.subTest(after=after), self.assertRaises(ValueError):
                self.altered(before, after)

    def test_generation_does_not_overwrite_unmanaged_input(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "scf.in"
            path.write_text("user input\n")
            with patch("configure.render", return_value={path: "new\n"}):
                with self.assertRaises(FileExistsError):
                    configure()
            self.assertEqual(path.read_text(), "user input\n")

    def test_all_python_sources_parse_and_no_removed_dependencies(self):
        for path in ROOT.rglob("*.py"):
            if ".git" not in path.parts:
                ast.parse(path.read_text(), filename=str(path), feature_version=(3, 10))
        self.assertFalse((ROOT / "frozen_magnon").exists())
        self.assertFalse((ROOT / "wannier").exists())
        self.assertFalse((ROOT / "7_spin/native_pw").exists())
        for stage in ("0_QE", "1_scf", "2_band", "3_band_no_mag", "4_potential", "5_wannier", "6_RVR", "7_spin"):
            self.assertTrue((ROOT / stage).is_dir())
            for path in (ROOT / stage).rglob("*"):
                if path.suffix in (".py", ".f90", ".sh", ".pbs"):
                    if not path.stem.endswith("_lanio2"):
                        self.assertNotIn("THE.in", path.read_text(), str(path))
                    self.assertNotIn("frozen_magnon", path.read_text(), str(path))

    @unittest.skipUnless(shutil.which("gfortran") or shutil.which("ifort"), "Fortran compiler is unavailable")
    def test_generated_fortran_values_and_row_validation(self):
        with tempfile.TemporaryDirectory() as directory:
            work = Path(directory)
            (work / "parameter_mod.f90").write_text(fortran_module(load()))
            (work / "probe.f90").write_text("""program probe
use parameter_mod
implicit none
integer :: row
call read_row_index(row)
write(*,*) p_nwan, p_n_func, p_nr1, row
end program
""")
            compiler = shutil.which("gfortran") or shutil.which("ifort")
            subprocess.run([compiler, "parameter_mod.f90", "probe.f90", "-o", "probe"], cwd=work, check=True)
            self.assertEqual(subprocess.check_output([str(work/"probe"), "7"], text=True).split(), ["44", "36", "81", "7"])
            for argument in ("0", "bad", "89"):
                run = subprocess.run([str(work/"probe"), argument], capture_output=True)
                self.assertNotEqual(run.returncode, 0)

    @unittest.skipUnless(shutil.which("gfortran"), "gfortran is unavailable")
    def test_small_fortran_record_kernels(self):
        names = ("4_potential/make_R_vec.x", "4_potential/V_split.x", "5_wannier/make_wan/make_w_density.x")
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            scf, potential, wan = base/"1_scf", base/"4_potential", base/"5_wannier/make_wan"
            for path in (scf, potential, wan): path.mkdir(parents=True)
            p = load()
            p.update(scf_dir=str(scf), func_dir=str(potential), split_super_x=1, split_super_y=1, split_super_z=1)
            (base/"parameter_mod.f90").write_text(fortran_module(p))
            binaries = []
            for name in names:
                binary = base/Path(name).name
                subprocess.run(["gfortran", "-ffree-line-length-none", str(base/"parameter_mod.f90"),
                                str(ROOT/Path(name).with_suffix(".f90")), "-o", str(binary)], cwd=base, check=True)
                binaries.append(binary)

            nx, ny, nz = 2, 2, 3
            grid = struct.pack("<3i", nx, ny, nz)
            values = np.arange(nx*ny*nz*4, dtype="<f8").reshape((-1,4), order="F")
            write_records(scf/"V_HYC_dp", grid, values.tobytes(order="F"))
            write_records(scf/"V_QE", grid)
            extended = np.arange((nx+1)*(ny+1)*(nz+1)*3*(p["n_func"]//3)*3, dtype="<f8").reshape((-1,3,p["n_func"]//3,3), order="F")
            write_records(potential/"V_func_tmp_tmp", b"", b"", b"", extended.tobytes(order="F"), b"")
            subprocess.run([str(binaries[0])], cwd=potential, check=True, stdout=subprocess.DEVNULL)
            got = np.frombuffer(read_records(potential/"V_func_tmp")[3], dtype="<f8").reshape((-1,3,p["n_func"]//3,3), order="F")
            expected = extended.reshape((nx+1,ny+1,nz+1,3,p["n_func"]//3,3), order="F")[:nx,:ny,:nz].reshape((-1,3,p["n_func"]//3,3), order="F")
            np.testing.assert_array_equal(got, expected)
            shutil.copyfile(scf/"V_HYC_dp", potential/"V_HYC_dp")
            (potential/"super_atom").write_text("".join(f"{'Cr':14s}{0.:20.10f}{0.:20.10f}{0.:20.10f}\n" for _ in range(2)))
            (potential/"polar2").write_text("3 0 0\n3 180 0\n")
            subprocess.run([str(binaries[1])], cwd=potential, check=True, stdout=subprocess.DEVNULL)
            got = np.frombuffer(read_records(potential/"V1")[1], dtype="<f8").reshape((-1,4),order="F")
            ix = [x%nx+nx*(y%ny+ny*(z%nz)) for z in range(nz+1) for y in range(ny+1) for x in range(nx+1)]
            np.testing.assert_array_equal(got, values[ix])
            for band in range(1,45):
                wf = np.full((9,nx*ny*nz,1), complex(band,.5), dtype="<c16")
                write_records(wan/f"wan_{band}", wf.tobytes(order="F"))
            out = wan/"weight"
            subprocess.run([str(binaries[2]), str(out)], cwd=wan, check=True, stdout=subprocess.DEVNULL)
            records = read_records(out)
            self.assertTrue(records[0].startswith(b"CRI3_WANNIER_WEIGHT_V1"))
            np.testing.assert_allclose(np.frombuffer(records[3],dtype="<f8"), 1, atol=1e-14)


if __name__ == "__main__":
    unittest.main(verbosity=2)
