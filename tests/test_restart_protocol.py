"""Bounded restart faults: input drift, torn publication, stale barriers and QE retry."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from concurrent.futures import ThreadPoolExecutor
from unittest.mock import patch
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'7_spin'))
import common
import run_material as material
import complete_energy as energy


class RestartProtocolTests(unittest.TestCase):
    def test_completed_density_assembly_preserves_metadata_on_reentry(self):
        with tempfile.TemporaryDirectory() as tmp,patch.object(energy,'R',Path(tmp)),patch.object(energy,'NG',2):
            out=energy.folder(1,0);out.mkdir(parents=True)
            common.atomic_npz(out/'DENSITY.npz',N=np.zeros((1,10,10)))
            common.save(out/'ELECTRONIC.json',dict(F=0.,band=0.,bare=0.,count=0.))
            (out/'DENSITY.bin').write_bytes(bytes(28+32*2))
            before=common.input_digest(out/'DENSITY.bin')
            with patch.object(energy,'write_density',side_effect=AssertionError('rewritten')):
                energy.assemble_density(1,0,1.)
            self.assertEqual(common.input_digest(out/'DENSITY.bin'),before)

    def test_eight_workers_rendezvous(self):
        with tempfile.TemporaryDirectory() as tmp,patch.object(energy,'ATTEMPT',Path(tmp)):
            with ThreadPoolExecutor(8) as pool:
                list(pool.map(lambda z:energy.barrier('current',z),range(8)))
            self.assertEqual(len(list((Path(tmp)/'barriers/current').glob('*.json'))),8)

    def test_truncated_density_recomputed_then_reused(self):
        with tempfile.TemporaryDirectory() as tmp,patch.object(energy,'R',Path(tmp)),patch.object(energy,'NG',2),patch.object(energy,'points',return_value=[]):
            out=energy.folder(1,0);out.mkdir(parents=True)
            common.save(out/'MU.json',dict(mu=0.))
            (out/'density_z0.npz').write_bytes(b'partial')
            energy.density_stage(1,0,0,{},1.,{})
            self.assertTrue(common.valid_npz(out/'density_z0.npz',{'density':(1,2,4),'N':(1,10,10)}))
            with patch.object(energy.time,'monotonic',side_effect=AssertionError('recomputed')):
                energy.density_stage(1,0,0,{},1.,{})

    def test_atomic_npz_failure_keeps_previous_and_detects_truncation(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'result.npz';common.atomic_npz(path,e=np.arange(3))
            original=path.read_bytes()
            def fail(f,**kw):f.write(b'partial');raise RuntimeError('interrupted')
            with patch.object(common.np,'savez',side_effect=fail):
                with self.assertRaises(RuntimeError):common.atomic_npz(path,e=np.arange(4))
            self.assertEqual(path.read_bytes(),original)
            self.assertTrue(common.valid_npz(path,['e']))
            path.write_bytes(original[:40])
            self.assertFalse(common.valid_npz(path,['e']))
            self.assertFalse(list(Path(tmp).glob('*.partial')))

    def test_spin_binding_precedes_completed_result_reuse(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'inputs/L02').mkdir(parents=True)
            np.savez(root/'inputs/L02/structure.npz',spins=np.array([[0.,0.,1.]]))
            case={'perturb_site_native_1based':None,'angle_deg':0}
            out=root/'results'
            with patch.object(material,'ROOT',root),patch.object(material,'input_digest',return_value='fixture'),patch.object(material,'prediction_identity',return_value={}):
                material.bind_case(2,{},case,out,2,.001)
                (out/'TORQUE.json').write_text('{"status":"COMPLETE"}')
                material.bind_case(2,{},case,out,2,.001)
                np.savez(root/'inputs/L02/structure.npz',spins=np.array([[1.,0.,0.]]))
                with self.assertRaisesRegex(ValueError,'inputs changed'):
                    material.bind_case(2,{},case,out,2,.001)
                (out/'INPUTS.json').unlink()
                with self.assertRaisesRegex(ValueError,'unknown prediction rule'):
                    material.bind_case(2,{},case,out,2,.001)
            self.assertTrue((out/'TORQUE.json').exists())

    def test_coefficients_restart_and_external_validation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);model=root/'V_func_tmp';model.mkdir();(model/'V_FUNC.json').write_text('{}')
            source=root/'source';source.mkdir();(source/'local_vfunc.py').write_text('fixture')
            structure=root/'structure.npz';np.savez(structure,spins=np.array([[0.,0.,1.]]))
            out=root/'coeff';out.mkdir();(out/'value.bin').write_bytes(b'partial')
            def evaluate(*args):
                stage=Path(args[args.index('--output')+1]);stage.mkdir()
                (stage/'PREDICTION.json').write_text('{"status":"COMPLETE"}')
                (stage/'value.bin').write_bytes(b'value');(stage/'dsite_0001.bin').write_bytes(b'derivative')
            with patch.object(material,'INTERFACE',root),patch.object(material,'SOURCE',source),patch.object(material,'prediction_identity',return_value={}),patch.object(material,'run',side_effect=evaluate) as run:
                material.coefficients(structure,out,[1],.001)
                material.coefficients(structure,out,[1],.001)
                self.assertEqual(run.call_count,1)
                (out/'dsite_0001.bin').write_bytes(b'broken')
                material.coefficients(structure,out,[1],.001)
                self.assertEqual(run.call_count,2)
                with np.load(structure) as a:identity=material.coefficient_identity(dict(a),[1],.001)
                self.assertTrue(material.valid_coefficients(out,identity))
                with patch.object(material,'prediction_identity',return_value={'model.py':'new revision'}):
                    with self.assertRaisesRegex(ValueError,'inputs changed'):
                        material.coefficients(structure,out,[1],.001)
                identity['structure']='different spin'
                with self.assertRaisesRegex(ValueError,'inputs changed'):material.valid_coefficients(out,identity)

    def test_stale_barriers_and_failure_do_not_release_worker(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);old=root/'attempts/old';new=root/'attempts/new';new.mkdir(parents=True)
            for z in range(8):common.save(old/f'barriers/phase/{z}.json',{'done':True})
            common.save(old/'FAILURE.json',{'error':'old failure'})
            with patch.object(energy,'R',root),patch.object(energy,'ATTEMPT',new),patch.object(energy.time,'sleep',side_effect=RuntimeError('still waiting')):
                with self.assertRaisesRegex(RuntimeError,'still waiting'):energy.barrier('phase',0)
                for z in range(1,8):common.save(new/f'barriers/phase/{z}.json',{'done':True})
                energy.barrier('phase',0)
                common.save(new/'FAILURE.json',{'error':'new failure'})
                with self.assertRaisesRegex(RuntimeError,'still waiting'):
                    energy.wait(new/'barriers/phase/missing.json')

    def test_launch_id_can_be_reused(self):
        with tempfile.TemporaryDirectory() as tmp,patch.object(energy,'R',Path(tmp)),patch.object(energy,'barrier'),patch.dict(os.environ,{'ANSATZ_RUN_ID':'test-launch'}):
            try:
                energy.start_attempt(0)
                energy._WORKER_LOCK.close()
                energy.start_attempt(0)
                self.assertTrue((energy.ATTEMPT/'worker_0.started').is_file())
            finally:
                if energy._WORKER_LOCK:energy._WORKER_LOCK.close()
                energy.ATTEMPT=None

    def test_automatic_launch_id_shared_by_workers(self):
        with tempfile.TemporaryDirectory() as tmp,patch.object(energy,'R',Path(tmp)),patch.object(energy,'barrier'),patch.dict(os.environ,{},clear=True):
            try:
                energy.start_attempt(0);first=energy.ATTEMPT
                energy._WORKER_LOCK.close()
                energy.start_attempt(1)
                self.assertEqual(energy.ATTEMPT,first)
                self.assertEqual((Path(tmp)/'run_ids'/f'local-{os.getppid()}').read_text().strip(),first.name)
            finally:
                if energy._WORKER_LOCK:energy._WORKER_LOCK.close()
                energy.ATTEMPT=None

    def test_qe_partial_retry_and_summary_reentry(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);binary=root/'qe';binary.write_bytes(b'qe fixture')
            for i in range(11):
                p=root/f'L04/case_{i:02}';p.mkdir(parents=True)
                (p/'DENSITY.bin').write_bytes(b'density')
                np.savez(p/'DENSITY.npz',N=np.eye(10)[None]*.3)
                common.save(p/'ELECTRONIC.json',dict(F=1.,band=2.,bare=3.))
            calls=[]
            def run(*args,**kw):
                cwd=kw['cwd'];self.assertFalse((cwd/'CORE_INPUT.bin').exists())
                (cwd/'CORE_INPUT.bin').write_bytes(b'partial')
                calls.append(cwd)
                if len(calls)==1:raise subprocess.CalledProcessError(1,args[0])
                for i in range(11):
                    np.savetxt(cwd/f'case_{i:02}_ENERGY.dat',np.arange(8.))
                    (cwd/f'case_{i:02}_TRUE.bin').write_bytes(b'field')
                kw['stdout'].write('GGA_BRANCH_REEXPORT_COMPLETE');kw['stdout'].flush()
            with patch.object(energy,'R',root),patch.dict(energy.ENERGY,{'native_binary':str(binary)}),patch.object(energy,'input_file'),patch.object(energy,'read_density',return_value=np.zeros((1,4))),patch.object(energy.subprocess,'run',side_effect=run):
                ueff=2.75
                with self.assertRaises(subprocess.CalledProcessError):energy.native(4,ueff)
                self.assertFalse((root/'L04/native/COMPLETE.json').exists())
                energy.native(4,ueff)
                (root/'L04/case_00/ENERGY.json').unlink()
                energy.native(4,ueff)
            self.assertEqual(len(calls),2)
            restored=json.loads((root/'L04/case_00/ENERGY.json').read_text())
            self.assertEqual(restored['Ueff_eV'],ueff)
            self.assertAlmostEqual(restored['terms']['EU'],.5*ueff*10*(.3-.3**2))
            self.assertEqual(len(list((root/'L04').glob('native.incomplete-*'))),1)


if __name__=='__main__':unittest.main()
