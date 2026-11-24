"""Analytic Heisenberg oracle and synthetic K4 interfaces; no material jobs."""
from dataclasses import dataclass
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO/'7_spin'))
from tools.spin_coordinates import CHART_KEYS, chart_spins, curvature_cases, derivative_files, tangent_basis
import local_vfunc
import run_material
import complete_energy
import summarize_energy
from torque_contract import derivative


def chart_jacobian(base, basis, eta):
    """Analytic differential, independent of the production finite differences."""
    r = np.linalg.norm(eta)
    if r == 0:
        return basis
    return (np.sin(r)/r*basis - np.outer(base, np.sin(r)/r*eta)
            + np.outer(basis@eta, (r*np.cos(r)-np.sin(r))/r**3*eta))


@dataclass
class Structure:
    spins: np.ndarray
    ni: np.ndarray


class CurvatureTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(245)
        self.base = rng.normal(size=(4,3))
        self.base /= np.linalg.norm(self.base,axis=1)[:,None]
        self.sites = [2,4]
        self.basis = tangent_basis(self.base[np.array(self.sites)-1])
        self.J = np.array([[0.,.8,-.2,.4],[.8,0.,.3,1.2],[-.2,.3,0.,-.6],[.4,1.2,-.6,0.]])
        self.steps = [.004,.008,.016]
        self.definition = dict(kind='K4',sites=self.sites,outer_steps_rad=self.steps,inner_step_rad=.001,
                               angle_configurations=curvature_cases(self.sites,self.steps))

    def evaluate(self, eta):
        spins = chart_spins(self.base,np.array(self.sites),self.basis,eta)
        gradient = self.J@spins
        tau = [-chart_jacobian(self.base[s-1],self.basis[k],eta[k]).T@gradient[s-1]
               for k,s in enumerate(self.sites)]
        return .5*np.sum(spins*gradient), np.array(tau).ravel()

    def analytic(self):
        gradient = self.J@self.base
        result = np.zeros((4,4))
        for a,i in enumerate(self.sites):
            for b,j in enumerate(self.sites):
                result[2*a:2*a+2,2*b:2*b+2] = (-np.dot(self.base[i-1],gradient[i-1])*np.eye(2)
                    if i==j else self.J[i-1,j-1]*self.basis[a].T@self.basis[b])
        return result

    def samples(self):
        rows = [self.evaluate(np.array(c['eta'])) for c in self.definition['angle_configurations']]
        return [r[0] for r in rows], [r[1] for r in rows]

    def test_nonstationary_heisenberg_all_blocks_and_energy_stencil(self):
        energies, torques = self.samples()
        result = summarize_energy.curvature_result(self.definition,energies,torques)
        expected = self.analytic()
        self.assertGreater(np.linalg.norm(result['tau0_eV_per_rad']),.1)
        np.testing.assert_allclose(result['K_eV_per_rad2'],expected,rtol=0,atol=2e-9)
        np.testing.assert_allclose(result['energy_diagonal_eV_per_rad2'],np.diag(expected),rtol=0,atol=2e-9)
        self.assertLess(result['max_reciprocity_error_eV_per_rad2'],2e-9)
        self.assertFalse(result['symmetrized'])
        self.assertEqual(result['modes'],[[2,1],[2,2],[4,1],[4,2]])
        # Independent energy differences include same-site mixed and intersite components.
        h=1e-4
        for a in range(4):
            for b in range(a+1,4):
                value=0.
                for sa,sb in ((1,1),(1,-1),(-1,1),(-1,-1)):
                    eta=np.zeros(4);eta[a]=sa*h;eta[b]=sb*h
                    value+=sa*sb*self.evaluate(eta.reshape(2,2))[0]/(4*h*h)
                self.assertAlmostEqual(value,expected[a,b],delta=1e-7)

    def test_steps_reduce_error_and_geometry_changes_K_with_fixed_J(self):
        result=summarize_energy.curvature_result(self.definition,*self.samples())
        errors=[np.max(abs(np.array(k)-self.analytic())) for k in result['K_steps_eV_per_rad2']]
        self.assertGreater(errors[2],3.5*errors[1])
        self.assertGreater(errors[1],3.5*errors[0])
        original=self.analytic()
        self.base=chart_spins(self.base,np.array(self.sites),self.basis,np.array([[.2,.1],[-.1,.15]]))
        self.basis=tangent_basis(self.base[np.array(self.sites)-1])
        self.assertGreater(np.max(abs(self.analytic()-original)),.01)

    def test_collinear_signs_and_two_rotation_zero_modes(self):
        self.J=np.array([[0.,1.3],[1.3,0.]])
        self.sites=[1,2]
        for sign in (1,-1):
            self.base=np.array([[0.,0.,1.],[0.,0.,float(sign)]])
            self.basis=tangent_basis(self.base)
            self.definition.update(sites=self.sites,angle_configurations=curvature_cases(self.sites,self.steps))
            result=summarize_energy.curvature_result(self.definition,*self.samples())
            np.testing.assert_allclose(result['K_eV_per_rad2'],self.analytic(),atol=1e-9)
            eigenvalues=np.linalg.eigvalsh(result['K_eV_per_rad2'])
            self.assertEqual(np.count_nonzero(abs(eigenvalues)<1e-8),2)
            np.testing.assert_allclose(np.diag(self.analytic()),-sign*1.3)

    def test_fixed_chart_norm_global_rotation_and_basis_covariance(self):
        eta=np.array([[.1,-.2],[.05,.3]])
        spins=chart_spins(self.base,np.array(self.sites),self.basis,eta)
        np.testing.assert_allclose(np.linalg.norm(spins,axis=1),1.,atol=1e-14)
        np.testing.assert_array_equal(spins[[0,2]],self.base[[0,2]])
        rotation=np.array([[0.,0.,1.],[1.,0.,0.],[0.,1.,0.]])
        rotated=chart_spins(self.base@rotation.T,np.array(self.sites),rotation@self.basis,eta)
        np.testing.assert_allclose(rotated,spins@rotation.T,atol=1e-14)
        before=self.analytic()
        q=np.array([[.6,-.8],[.8,.6]])
        self.basis=self.basis@q
        block=np.kron(np.eye(2),q)
        np.testing.assert_allclose(self.analytic(),block.T@before@block,atol=1e-14)

    def test_invalid_inputs_and_incomplete_results_are_rejected(self):
        for sites,steps in (([2,2],self.steps),([0,4],self.steps),([2.,4.],self.steps),
                            (self.sites,[.01,.02,.05]),(self.sites,[0.,0.,0.])):
            with self.assertRaises(ValueError):curvature_cases(sites,steps)
        with self.assertRaises(ValueError):
            chart_spins(self.base,np.array(self.sites),self.basis*2,np.zeros((2,2)))
        with self.assertRaises(ValueError):summarize_energy.curvature_result(self.definition,[0.]*24,np.zeros((24,4)))
        with self.assertRaises(ValueError):summarize_energy.curvature_result(self.definition,[0.]*25,np.full((25,4),np.nan))
        bad=dict(self.definition,angle_configurations=self.definition['angle_configurations'][::-1])
        with self.assertRaises(ValueError):summarize_energy.curvature_result(bad,*self.samples())

    def test_evaluate_exports_both_directions_at_displaced_chart_point(self):
        eta=np.array([[.12,-.07],[.03,.09]])
        spins=chart_spins(self.base,np.array(self.sites),self.basis,eta)
        structure=Structure(spins,np.arange(4))
        def predict(s):
            p=np.zeros((4,1,4));p[:,0,0]=s.spins[:,0]+2*s.spins[:,1];p[:,0,1:]=s.spins
            n=np.tile(np.eye(10,dtype=complex)*.4,(4,1,1))
            n[:,0,1]=.02*(s.spins[:,0]+1j*s.spins[:,1]);n[:,1,0]=n[:,0,1].conj()
            return SimpleNamespace(p=p,N=n)
        model=SimpleNamespace(predict=predict,physical_p=lambda p:p.p)
        runtime=SimpleNamespace(Bundle=lambda path:SimpleNamespace(model=model))
        meta=dict(schema=local_vfunc.SCHEMA,ncoeff=1,model_id='a'*64,Ueff_eV=1.,scalar_offset_eV=0.)
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'s.npz'
            np.savez(source,spins=spins,**dict(zip(CHART_KEYS,(self.base,np.array(self.sites),self.basis,eta))))
            with patch.object(local_vfunc,'metadata',return_value=meta),patch.object(local_vfunc,'runtime',return_value=(runtime,None)), \
                 patch.object(local_vfunc,'read_structure',return_value=(structure,2,np.arange(4))):
                report=local_vfunc.evaluate(root,source,root/'out',self.sites,.001,tangent=True)
                with self.assertRaisesRegex(ValueError,'--tangent'):
                    local_vfunc.evaluate(root,source,root/'y_axis',self.sites,.001)
            self.assertEqual(report['derivative_modes'],[(2,1),(2,2),(4,1),(4,2)])
            for k,site in enumerate(self.sites):
                jac=chart_jacobian(self.base[site-1],self.basis[k],eta[k])
                for alpha in (1,2):
                    data=derivative(root/'out'/f'dsite_{site:04}_a{alpha}.bin')
                    expected=np.zeros((1,4,4))
                    expected[0,site-1,0]=jac[0,alpha-1]+2*jac[1,alpha-1]
                    expected[0,site-1,1:]=jac[:,alpha-1]
                    np.testing.assert_allclose(data['p'],expected,rtol=0,atol=2e-10)
                    self.assertAlmostEqual(data['N'][site-1,0,1],.02*(jac[0,alpha-1]+1j*jac[1,alpha-1]),delta=2e-11)

    def test_integer_chart_origin_preserves_small_derivatives(self):
        def predict(s):
            p=np.zeros((4,1,4));p[:,0,1:]=s.spins
            return SimpleNamespace(p=p,N=np.zeros((4,10,10),dtype=complex))
        model=SimpleNamespace(predict=predict,physical_p=lambda p:p.p)
        runtime=SimpleNamespace(Bundle=lambda path:SimpleNamespace(model=model))
        meta=dict(schema=local_vfunc.SCHEMA,ncoeff=1,model_id='a'*64,Ueff_eV=1.,scalar_offset_eV=0.)
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'s.npz'
            np.savez(source,spins=self.base,**dict(zip(CHART_KEYS,
                (self.base,np.array(self.sites),self.basis,np.zeros((2,2),dtype=int)))))
            with patch.object(local_vfunc,'metadata',return_value=meta), \
                 patch.object(local_vfunc,'runtime',return_value=(runtime,None)), \
                 patch.object(local_vfunc,'read_structure',return_value=(Structure(self.base,np.arange(4)),2,np.arange(4))):
                local_vfunc.evaluate(root,source,root/'out',self.sites,.001,tangent=True)
            for k,site in enumerate(self.sites):
                for alpha in (1,2):
                    data=derivative(root/'out'/f'dsite_{site:04}_a{alpha}.bin')
                    np.testing.assert_allclose(data['p'][0,site-1,1:],self.basis[k,:,alpha-1],atol=1e-11)

    def test_native_input_reordering_preserves_chart_site_meaning(self):
        cell=np.diag([2.,2.,1.]);positions=np.array([[0.,0.,0.],[0.,1.,0.],[1.,0.,0.],[1.,1.,0.]])
        permutation=np.array([3,0,2,1])
        meta=dict(primitive_cell_A=np.eye(3).tolist(),primitive_motif=[[28,0.,0.,0.]],Ni_fractional=[0.,0.,0.])
        @dataclass
        class Geometry:
            cell: np.ndarray
            positions: np.ndarray
            species: np.ndarray
            ni: np.ndarray
            spins: np.ndarray
            def validate(self):pass
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'inputs/L02';source.mkdir(parents=True)
            np.savez(source/'structure.npz',cell=cell,positions=positions,species=np.full(4,28),
                     ni=permutation,spins=self.base[permutation])
            with patch.object(run_material,'ROOT',root),patch.object(local_vfunc,'metadata',return_value=meta), \
                 patch.object(local_vfunc,'runtime',return_value=(None,SimpleNamespace(Structure=Geometry))):
                result=run_material.spin_structure(2,self.definition['angle_configurations'][1])
            np.testing.assert_array_equal(result['chart_base_spins'],self.base)
            np.testing.assert_array_equal(result['ni'],np.arange(4))
            np.testing.assert_allclose(result['spins'],chart_spins(self.base,np.array(self.sites),self.basis,
                                                                 np.array(self.definition['angle_configurations'][1]['eta'])))

    def test_checkpoint_binds_chart_and_mode_and_reuses_only_matching_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'V_func_tmp').mkdir();(root/'V_func_tmp/V_FUNC.json').write_text('{}')
            source=root/'s.npz'
            chart=dict(zip(CHART_KEYS,(self.base,np.array(self.sites),self.basis,np.zeros((2,2)))))
            np.savez(source,spins=self.base,**chart)
            def run(*args):
                self.assertIn('--tangent',args)
                out=Path(args[args.index('--output')+1]);out.mkdir()
                for name in ['PREDICTION.json','value.bin',*derivative_files(self.sites,True)]:
                    (out/name).write_text('fixture')
            with patch.object(run_material,'INTERFACE',root),patch.object(run_material,'prediction_identity',return_value={}), \
                 patch.object(run_material,'run',side_effect=run) as call:
                run_material.coefficients(source,root/'coeff',self.sites,.001,tangent=True)
                run_material.coefficients(source,root/'coeff',self.sites,.001,tangent=True)
                self.assertEqual(call.call_count,1)
                with self.assertRaisesRegex(ValueError,'inputs changed'):
                    run_material.coefficients(source,root/'coeff',self.sites,.001)
                chart['chart_basis']=-self.basis
                np.savez(source,spins=self.base,**chart)
                with self.assertRaisesRegex(ValueError,'inputs changed'):
                    run_material.coefficients(source,root/'coeff',self.sites,.001,tangent=True)

    def test_collect_and_csv_use_all_modes_and_reject_swapped_modes(self):
        energies,torques=self.samples()
        labels=[[s,a] for s in self.sites for a in (1,2)]
        spec=dict(sites=self.sites,outer_steps_rad=self.steps,inner_step_rad=.001)
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            for case,(energy,torque) in enumerate(zip(energies,torques)):
                path=root/f'L02/case_{case:02}';path.mkdir(parents=True)
                (path/'ENERGY.json').write_text(json.dumps(dict(energy=energy)))
                for k in range(32):np.savez(path/f'torque_k{k:04}.npz',sites=labels,torque=torque,imaginary=0.)
            with patch.object(complete_energy,'R',root),patch.dict(complete_energy.ENERGY,{'curvature':spec}), \
                 patch.object(complete_energy,'spin_structure',return_value={'chart_base_spins':self.base,'chart_basis':self.basis}):
                complete_energy.collect(2)
                result=json.loads((root/'L02/RESULT.json').read_text())
                np.testing.assert_allclose(result['K_eV_per_rad2'],self.analytic(),atol=2e-9)
                self.assertFalse(result['physical_validation'])
                with patch.object(summarize_energy,'ROOT',root),patch.dict(summarize_energy.CONFIG['complete_energy'],{'sizes':[2]}),patch('builtins.print'):
                    summarize_energy.main()
                np.testing.assert_allclose(np.loadtxt(root/'L02/K4_eV_per_rad2.csv',delimiter=','),self.analytic(),atol=2e-9)
                np.savez(root/'L02/case_00/torque_k0000.npz',sites=labels[::-1],torque=torques[0],imaginary=0.)
                with self.assertRaisesRegex(ValueError,'Torque modes'):
                    complete_energy.collect(2)

    def test_torque_stage_batches_four_modes_and_reuses_completed_points(self):
        _,torques=self.samples()
        spec=dict(sites=self.sites,outer_steps_rad=self.steps)
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);ops={}
            for k in (0,8,16,24):
                ops[k]=root/f'op{k}';ops[k].mkdir()
                np.save(ops[k]/'native_projectors.npy',np.zeros((1,1,1)))
            for case in range(25):
                path=root/f'L02/case_{case:02}';path.mkdir(parents=True)
                (path/'MU.json').write_text('{"mu":0.0}')
                for z in range(8):(path/f'alpha_z{z}.json').write_text('{"numerator":0.0,"denominator":-1.0}')
                for k in ops:
                    np.savez(path/f'k_{k:04}.npz',e=np.array([-.1,.1]),Q=np.eye(2))
                    np.savez(path/f'adjoint_{k:04}.npz',Ae=np.eye(2))
            def load(path):
                case=int(path.parent.parent.name.split('_')[1])
                return dict(Ueff=1.,case=case,name=path.name)
            batches=[]
            def contract(rvr,background,points,densities,derivatives,**kwargs):
                batches.append(len(points))
                answer=[]
                for ds in kwargs['derivatives_by_point']:
                    self.assertEqual([d['name'] for d in ds],derivative_files(self.sites,True))
                    answer.append(torques[ds[0]['case']])
                return answer,np.zeros(len(points)),{}
            with patch.object(complete_energy,'R',root),patch.dict(complete_energy.ENERGY,{'curvature':spec,'reuse_cases':{}}), \
                 patch.object(complete_energy,'derivative',side_effect=load), \
                 patch.object(complete_energy,'fold_unitary',return_value=np.eye(2)), \
                 patch.object(complete_energy,'contract_batch',side_effect=contract) as called:
                complete_energy.torque_stage(2,0,root/'rvr',ops,1.)
                self.assertEqual(batches,[16]*6+[4])
                with np.load(root/'L02/case_24/torque_k0024.npz') as data:
                    np.testing.assert_array_equal(data['sites'],[[2,1],[2,2],[4,1],[4,2]])
                    np.testing.assert_allclose(data['torque'],torques[24])
                called.reset_mock()
                complete_energy.torque_stage(2,0,root/'rvr',ops,1.)
                called.assert_not_called()

    def test_native_adapter_receives_variable_case_count(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);binary=root/'native';binary.write_text('fixture')
            for case in range(25):
                path=root/f'L02/case_{case:02}';path.mkdir(parents=True)
                (path/'DENSITY.bin').write_bytes(b'fixture')
                np.savez(path/'DENSITY.npz',N=np.eye(10)[None]*.4)
                (path/'ELECTRONIC.json').write_text('{"F":1.0,"band":2.0,"bare":3.0}')
            def run(*args,**kwargs):
                stage=kwargs['cwd'];lines=(stage/'input_list.txt').read_text().splitlines()
                self.assertEqual(lines[0],'25');self.assertEqual(len(lines),26)
                for case in range(25):
                    self.assertTrue(lines[case+1].startswith(f'case_{case:02} '))
                    np.savetxt(stage/f'case_{case:02}_ENERGY.dat',np.zeros(8))
                    (stage/f'case_{case:02}_TRUE.bin').write_bytes(b'field')
                kwargs['stdout'].write('GGA_BRANCH_REEXPORT_COMPLETE')
            with patch.object(complete_energy,'R',root),patch.dict(complete_energy.ENERGY,{'native_binary':str(binary)}), \
                 patch.object(complete_energy,'input_file'),patch.object(complete_energy,'read_density',return_value=np.zeros((1,4))), \
                 patch.object(complete_energy.subprocess,'run',side_effect=run) as called:
                complete_energy.native(2,1.,25)
                complete_energy.native(2,1.,25)
                self.assertEqual(called.call_count,1)
            result=json.loads((root/'L02/case_24/ENERGY.json').read_text())
            self.assertAlmostEqual(result['energy'],3.2)


if __name__=='__main__':unittest.main()
