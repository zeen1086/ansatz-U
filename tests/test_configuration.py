"""Offline input, parameter propagation and build-module checks."""
from pathlib import Path
import ast
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch
from configure import configure, render
from parameters import ROOT, fortran_module, load


class ConfigurationTests(unittest.TestCase):
    def altered(self, before, after):
        text=(ROOT/'parameter').read_text()
        self.assertIn(before,text)
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'parameter'
            path.write_text(text.replace(before,after))
            return load(path)

    def test_parameters_reach_inputs(self):
        p=self.altered('nk1 = 8','nk1 = 6')
        out=render(p)
        self.assertIn('p_nk1 = 6',fortran_module(p))
        self.assertIn('6 8 8 0 0 0',out[ROOT/'1_scf/scf.in'])
        self.assertIn('mp_grid = 6 8 8',out[ROOT/'5_wannier/LaNiO2.win'])
        self.assertIn('last_k=384',out[ROOT/'5_wannier/wfck2r.in'])
        rows=out[ROOT/'5_wannier/nscf.in'].split('K_POINTS crystal\n')[1].splitlines()
        self.assertEqual(int(rows[0]),384)
        self.assertAlmostEqual(sum(float(row.split()[3]) for row in rows[1:]),1)
        p=self.altered('a0 = 7.48180367','a0 = 9.0')
        self.assertIn('celldm(1)=9.0',render(p)[ROOT/'2_band/band.in'])
        self.assertIn('0 4.762594898127 0',render(p)[ROOT/'5_wannier/LaNiO2.win'])
        p=self.altered('ni_magnetization = 1.0','ni_magnetization = 0.37')
        self.assertIn('starting_magnetization(2)=0.37',render(p)[ROOT/'1_scf/scf.in'])
        p=self.altered('nspin = 1','nspin = 2')
        self.assertIn('spinors = true',render(p)[ROOT/'5_wannier/LaNiO2.win'])
        p=self.altered('write_v_qe = true','write_v_qe = false')
        self.assertIn('write_v_qe=.false.',render(p)[ROOT/'1_scf/print_V.in'])

    def test_rejects_inconsistent_parameters(self):
        for old,new in [('nwan = 30','nwan = 20'),('convert_size_x = 1','convert_size_x = 4'),
                        ('nk3 = 8','nk3 = 0'),('nspin = 1','nspin = 3'),
                        ('split_nspinor = 1','split_nspinor = 0'),('ecutwfc = 100.0','ecutwfc = nan'),
                        ('export_grid = "smooth"','export_grid = "invalid"'),
                        ('num_bands = 50','num_bands = 49')]:
            with self.subTest(new=new),self.assertRaises(ValueError):self.altered(old,new)

    def test_generation_preserves_unmanaged_input(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'scf.in';path.write_text('user input\n')
            with patch('configure.render',return_value={path:'new\n'}):
                with self.assertRaises(FileExistsError):configure()
            self.assertEqual(path.read_text(),'user input\n')

    def test_all_python_sources_parse(self):
        for path in ROOT.rglob('*.py'):
            if '.git' not in path.parts:ast.parse(path.read_text(),filename=str(path),feature_version=(3,10))

    @unittest.skipUnless(shutil.which('gfortran'),'gfortran is unavailable')
    def test_generated_values_and_row_bounds(self):
        with tempfile.TemporaryDirectory() as directory:
            work=Path(directory)
            (work/'parameter_mod.f90').write_text(fortran_module(load()))
            (work/'probe.f90').write_text('program probe\nuse parameter_mod\ninteger::row\ncall read_row_index(row)\nprint *,p_nwan,p_split_nspinor,p_nr1,row\nend\n')
            subprocess.run(['gfortran','parameter_mod.f90','probe.f90','-o','probe'],cwd=work,check=True)
            self.assertEqual(subprocess.check_output([str(work/'probe'),'7'],text=True).split(),['30','1','64','7'])
            for arg in ('0','bad','61'):
                self.assertNotEqual(subprocess.run([str(work/'probe'),arg],capture_output=True).returncode,0)


if __name__=='__main__':unittest.main(verbosity=2)
