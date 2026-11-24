"""Reviewed native H, common-mu HF torque and shell-J execution; no scheduler calls."""
import argparse
import csv
import json
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import numpy as np
from scipy.spatial.transform import Rotation
from scipy.special import erfc
from common import CONFIG, INTERFACE, BUILD, REF, ROOT, WORK, SPINS, PYTHON, SIGMA, SOURCE
from common import compute_node, run, save, load_h, chemical_potential, check_regression, quality_check
from common import atomic_npz, valid_npz, valid_json, archive_incomplete, bind_inputs, input_digest, arrays_digest
from torque_contract import derivative,contract_batch
sys.path.insert(0, str(SOURCE))
from tools.spin_coordinates import chart_spins, tangent_basis, derivative_files


def binding() -> None:
    sys.path.insert(0,str(SOURCE))
    from local_vfunc import metadata
    meta=metadata(INTERFACE/'V_func_tmp')
    basis=json.loads((WORK/'RVR/PROJECTION.json').read_text())
    if (meta['model_id']!=CONFIG['model_id'] or basis['model_id']!=meta['model_id']):
        raise ValueError('Model/background provenance changed')
    with (WORK/'RVR/B_general').open('rb') as f: head=f.read(176)
    if head[:16]!=b'LANIO2_LRVR_V2  ' or head[48:112].decode()!=meta['model_id'] or head[112:176].decode()!=basis['basis_id']:
        raise ValueError('RVR identity mismatch')


def run_native(coeff: Path, k: np.ndarray, out: Path, components: bool, *, rvr: Path | None=None, background_files: list[Path] | None=None) -> float:
    packed=np.column_stack((k[:,:3],np.full(len(k),1/len(k))))
    np.savetxt(out/'k.dat',packed,header=str(len(k)),comments='',fmt='%.17g')
    values=dict(rvr_path=rvr or WORK/'RVR/B_general',background_path=WORK/'basis/pure_H0.au',
                coefficient_path=coeff,kpoints_path=out/'k.dat',output_prefix=out/'H')
    if background_files is not None:
        if len(background_files)!=len(k):raise ValueError('One background cache required per k point')
        listing=out/'background_cache.txt'
        listing.write_text(''.join(str(p.resolve())+'\n' for p in background_files))
        values['background_cache_list']=listing
    lines=['&local_hamiltonian']
    lines += [f" {key}='{str(value.resolve()).replace(chr(39),chr(39)*2)}'," for key,value in values.items()]
    lines += [" background_kind='pure_H0_no_local_fields_no_HU',",
              f" dump_components={'.true.' if components else '.false.'}",'/']
    (out/'run.nml').write_text('\n'.join(lines)+'\n');start=time.monotonic()
    env=os.environ.copy()
    workers=int(env.get('V1_K_WORKERS','1'))
    if len(k)>1 and workers>1:
        env.update(OMP_NUM_THREADS=str(workers),MKL_NUM_THREADS=str(workers))
    subprocess.run([BUILD/'local/cal_Hsk',out/'run.nml'],env=env,check=True)
    return time.monotonic()-start


def native(coeff: Path, k: np.ndarray, out: Path, components: bool=False) -> Path:
    # The derivative header carries the same supercell size as the value coefficients.
    with coeff.open('rb') as f:
        f.seek(16); length=int(np.fromfile(f,'<i4',4)[2])
    out.parent.mkdir(parents=True,exist_ok=True)
    native_batch(coeff,k[None],[out],out.parent,60*length*length,components)
    return out/'H.k000001.h'


def native_batch(coeff: Path, k: np.ndarray, directories: list[Path], base: Path,
                 dimension: int, components: bool=False, *, rvr: Path | None=None, background_files: list[Path] | None=None) -> None:
    """Assemble one kz group and publish per-point checkpoints atomically."""
    if len(k)!=len(directories) or not 1<=len(k)<=16 or np.ptp(k[:,2])>1e-10:
        raise ValueError('Expected one kz group of up to 16 points and matching output paths')
    if background_files is not None and len(background_files)!=len(k):
        raise ValueError('One background cache required per k point')
    missing=[]
    suffixes=['h','h0','hv','hb','hu'] if components else ['h']
    expected=16*dimension**2
    basis=json.loads((WORK/'RVR/PROJECTION.json').read_text())['basis_id']
    ids=[CONFIG['model_id'],basis]
    identity=dict(coeff=input_digest(coeff),rvr=input_digest(rvr or WORK/'RVR/B_general'),
                  background=input_digest(WORK/'basis/pure_H0.au'),binary=input_digest(BUILD/'local/cal_Hsk'),
                  dimension=dimension,components=components)
    identities=[dict(identity,k=list(map(float,point[:3])),background_cache=None if background_files is None else input_digest(background_files[i])) for i,point in enumerate(k)]
    for index,out in enumerate(directories):
        if (out/'COMPLETE.json').is_file():
            try:
                recorded=json.loads((out/'COMPLETE.json').read_text())
            except (OSError,ValueError):recorded={}
            if recorded and recorded.get('inputs') != identities[index]:
                raise ValueError(f'Completed Hamiltonian inputs changed or unbound: {out}')
            names=['H.index',*[f'H.k000001.{suffix}' for suffix in suffixes]]
            if recorded.get('files') and set(recorded['files'])==set(names) and all(
                    (out/name).is_file() and input_digest(out/name)==recorded['files'][name] for name in names):
                continue
        archive_incomplete(out)
        missing.append(index)
    if not missing:
        return
    with tempfile.TemporaryDirectory(prefix='.cal_Hsk_batch_',dir=base) as tmp:
        staged=Path(tmp)
        if background_files is not None:
            seconds=run_native(coeff,k[missing],staged,components,rvr=rvr,
                               background_files=[background_files[i] for i in missing])
        elif rvr is None:
            seconds=run_native(coeff,k[missing],staged,components)
        else:
            seconds=run_native(coeff,k[missing],staged,components,rvr=rvr)
        metadata=(staged/'H.index').read_text().splitlines()
        if len(metadata)!=4+len(missing) or [s.strip() for s in metadata[2:4]]!=ids:
            raise ValueError('Incomplete cal_Hsk batch index')
        for slot,index in enumerate(missing,1):
            row=metadata[3+slot].split()
            if int(row[0])!=slot or not np.allclose(np.array(row[1:4],float),k[index,:3],rtol=0,atol=1e-12):
                raise ValueError('cal_Hsk batch point order changed')
            point=staged/f'point_{slot}';point.mkdir()
            for suffix in suffixes:
                part=staged/f'H.k{slot:06}.{suffix}'
                if part.stat().st_size!=expected:
                    raise ValueError(f'Truncated cal_Hsk batch output: {part}')
                part.rename(point/f'H.k000001.{suffix}')
            head=metadata[1].split();head[2]='1';row[0]='1';row[4]='1.0'
            (point/'H.index').write_text('\n'.join([metadata[0],' '.join(head),*metadata[2:4],' '.join(row)])+'\n')
            np.savetxt(point/'k.dat',np.r_[k[index,:3],1.][None],header='1',comments='',fmt='%.17g')
            save(point/'COMPLETE.json',dict(inputs=identities[index],files={name:input_digest(point/name) for name in ['H.index',*[f'H.k000001.{suffix}' for suffix in suffixes]]},seconds=seconds/len(missing),batch_seconds=seconds,batch_size=len(missing)))
        for slot,index in enumerate(missing,1):
            out=directories[index];out.parent.mkdir(parents=True,exist_ok=True)
            (staged/f'point_{slot}').rename(out)


def kz_groups(k: np.ndarray):
    for kz in np.unique(k[:,2]):
        members=np.flatnonzero(k[:,2]==kz)
        for start in range(0,len(members),16):
            yield members[start:start+16]


def prediction_identity():
    if str(SOURCE) not in sys.path:
        sys.path.insert(0,str(SOURCE))
    from local_vfunc import runtime_names, runtime_path
    model = INTERFACE/'V_func_tmp'
    meta = json.loads((model/'V_FUNC.json').read_text())
    code = runtime_path(model, meta)
    return {f:input_digest(code/f) for f in runtime_names(meta)}


def coefficient_identity(structure, sites, step, *, tangent=False):
    identity = dict(structure=arrays_digest(structure),sites=list(sites),step=step,
                model=input_digest(INTERFACE/'V_func_tmp/V_FUNC.json'),
                evaluator=input_digest(SOURCE/'local_vfunc.py'), runtime=prediction_identity())
    if tangent:
        identity.update(tangent=True, coordinates=input_digest(SOURCE/'tools/spin_coordinates.py'))
    return identity


def valid_coefficients(directory, identity):
    try:
        record=json.loads((directory/'CHECKPOINT.json').read_text())
        if record['inputs']!=identity:raise ValueError(f'Coefficient inputs changed: {directory}')
        expected=['PREDICTION.json','value.bin',*derivative_files(identity['sites'], identity.get('tangent', False))]
        return set(record['files'])==set(expected) and all(input_digest(directory/name)==record['files'][name] for name in expected)
    except (OSError,KeyError,json.JSONDecodeError):return False


def coefficients(structure: Path, directory: Path, sites: list[int], step: float, *, tangent=False) -> None:
    with np.load(structure) as data:identity=coefficient_identity(dict(data),sites,step,tangent=tangent)
    import fcntl
    directory.parent.mkdir(parents=True,exist_ok=True)
    with (directory.parent/(directory.name+'.lock')).open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if valid_coefficients(directory,identity):return
        archive_incomplete(directory)
        with tempfile.TemporaryDirectory(prefix='.coeff_',dir=directory.parent) as tmp:
            staged=Path(tmp)/'coeff'
            run(PYTHON,SOURCE/'local_vfunc.py','evaluate','--model',INTERFACE/'V_func_tmp',
                '--structure',structure,'--output',staged,'--derivative-sites',*sites,'--step',step,
                *(['--tangent'] if tangent else []))
            names=['PREDICTION.json','value.bin',*derivative_files(sites, tangent)]
            save(staged/'CHECKPOINT.json',dict(inputs=identity,files={name:input_digest(staged/name) for name in names}))
            staged.rename(directory)


def global_spin_rotation(case: dict) -> Rotation:
    """Active rotation after the local perturbation, in reference Cartesian axes."""
    degrees = np.asarray(case.get('global_spin_rotation_deg', [0., 0., 0.]), dtype=float)
    if degrees.shape != (3,) or not np.isfinite(degrees).all():
        raise ValueError('global_spin_rotation_deg must be a finite three-component rotation vector')
    return Rotation.from_rotvec(np.deg2rad(degrees))


def spin_structure(length, case):
    with np.load(ROOT/f'inputs/L{length:02}/structure.npz') as data:s=dict(data)
    if 'tangent_sites' in case:
        if 'global_spin_rotation_deg' in case:
            raise ValueError('K4 uses a fixed complement frame; rotate the entire input/reference package separately')
        from local_vfunc import metadata, read_structure, runtime
        model = INTERFACE/'V_func_tmp'; meta = metadata(model)
        _, features = runtime(model, meta)
        base, cell, _ = read_structure(ROOT/f'inputs/L{length:02}/structure.npz', meta, features)
        if cell != length:
            raise ValueError('Input structure and requested cell size differ')
        sites = np.asarray(case['tangent_sites'])
        if sites.shape != (2,) or sites.dtype.kind not in 'iu' or np.any(sites < 1) or np.any(sites > len(base.spins)):
            raise ValueError('Invalid native tangent sites')
        basis = tangent_basis(base.spins[sites-1]); eta = np.asarray(case['eta'], dtype=float)
        s.update(ni=base.ni, spins=chart_spins(base.spins, sites, basis, eta),
                 chart_base_spins=base.spins, chart_sites=sites, chart_basis=basis, chart_eta=eta)
        return s
    site=case['perturb_site_native_1based']
    if site is not None:
        angle=np.deg2rad(case['angle_deg']);v=s['spins'][site-1].copy()
        s['spins'][site-1]=[np.cos(angle)*v[0]+np.sin(angle)*v[2],v[1],-np.sin(angle)*v[0]+np.cos(angle)*v[2]]
    if 'global_spin_rotation_deg' in case:
        s['spins'] = global_spin_rotation(case).apply(s['spins'])
    return s


def bind_case(length, definition, case, out, mesh, step):
    if not (out/'INPUTS.json').exists() and any(
            any(out.glob(pattern)) for pattern in ('k_*','SPECTRA.*','TORQUE.json')):
        raise ValueError(f'Unbound spin results have unknown prediction rule: {out}; use a new output directory')
    files=[INTERFACE/'V_func_tmp/V_FUNC.json',WORK/'RVR/B_general',WORK/'RVR/PROJECTION.json',WORK/'basis/pure_H0.au',BUILD/'local/cal_Hsk']
    files += [Path(__file__).with_name(name) for name in ('run_material.py','common.py','kz_production.py','torque_contract.py')]
    bind_inputs(out,dict(structure=arrays_digest(spin_structure(length,case)),definition=definition,
                        case=case,mesh=mesh,step=step,runtime=prediction_identity(),
                        files={str(f):input_digest(f) for f in files}))


def regression_assemble() -> None:
    binding();data=np.load(REF/'data.npz'); base=dict(np.load(ROOT/'inputs/L02/structure.npz'))
    k=data['K'];out=SPINS/'regression';out.mkdir(parents=True,exist_ok=True)
    for name in ('Neel_anchor','fm'):
        row=out/name;row.mkdir(exist_ok=True)
        base['spins']=data['u'][data['names'].tolist().index(name)][[0,2,1,3]]
        atomic_npz(row/'structure.npz',**base)
        coefficients(row/'structure.npz',row/'coeff',[],.001)
        for group in kz_groups(k):
            directories=[row/f'k_{ki:04}' for ki in group]
            native_batch(row/'coeff/value.bin',k[group],directories,row,len(base['spins'])*60,components=True)
            print('REGRESSION_H_BATCH',name,group.tolist(),flush=True)


def case_root(length: int, index: int, mesh: int, step: float):
    definition=json.loads((ROOT/f'inputs/L{length:02}/CASE.json').read_text())
    case=definition['angle_configurations'][index]
    work=SPINS/f'J/L{length:02}/mesh{mesh}_h{step:g}/case_{index:02}'
    return definition,case,work


def map_k(k: np.ndarray,worker) -> None:
    workers=int(os.environ.get('V1_K_WORKERS','1'))
    if not 1<=workers<=16:raise ValueError('Expected 1..16 allocated k workers')
    if workers==1:
        for ki in np.argsort(k[:,2],kind='stable'):worker(int(ki))
        return
    # Keep one kz slab hot while independent native processes/NumPy calls use one core each.
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for kz in np.unique(k[:,2]):list(pool.map(worker,map(int,np.flatnonzero(k[:,2]==kz))))


def spectra(length: int, index: int, mesh: int, step: float) -> None:
    check_regression();binding();definition,case,out=case_root(length,index,mesh,step)
    bind_case(length,definition,case,out,mesh,step)
    out.mkdir(parents=True,exist_ok=True);start=time.monotonic()
    s=spin_structure(length,case)
    atomic_npz(out/'structure.npz',**s)
    coefficients(out/'structure.npz',out/'coeff',case['readout_sites_native_1based'],step)
    k=np.array([[x/mesh,y/mesh,z/8,1/(mesh*mesh*8)] for x in range(mesh) for y in range(mesh) for z in range(8)])
    energies=[None]*len(k);timing=[0.]*len(k);dimension=definition['dimension']
    def spectrum_k(ki):
        target=out/f'k_{ki:04}';tick=time.monotonic()
        hfile=target/'value/H.k000001.h'
        if not valid_npz(target/'eigensystem.npz',{'e':(dimension,),'v':(dimension,dimension)}):
            h=load_h(hfile,dimension);e,v=np.linalg.eigh(h)
            atomic_npz(target/'eigensystem.npz',e=e,v=v);del h,v
        with np.load(target/'eigensystem.npz') as z: energies[ki]=z['e']
        timing[ki]=time.monotonic()-tick;print('SPECTRUM',length,index,ki,timing[ki],flush=True)
    for group in kz_groups(k):
        directories=[out/f'k_{ki:04}'/'value' for ki in group]
        native_batch(out/'coeff/value.bin',k[group],directories,out,dimension)
        map_k(k[group],lambda i:spectrum_k(int(group[i])))
    energies=np.array(energies);mu=chemical_potential(energies,k[:,3],definition['electrons'])
    atomic_npz(out/'SPECTRA.npz',k=k,e=energies,mu=mu)
    save(out/'SPECTRA.json',dict(status='COMPLETE',L=length,case=case,mesh=[mesh,mesh,8],
         derivative_step_rad=step,seconds=time.monotonic()-start,k_seconds=timing,mu_eV=mu,
         electrons=definition['electrons'],model_id=CONFIG['model_id']))


def torque(length: int,index: int,mesh: int,step: float) -> None:
    check_regression();binding();definition,case,out=case_root(length,index,mesh,step)
    bind_case(length,definition,case,out,mesh,step)
    if valid_json(out/'TORQUE.json',['torque_eV','F_eV','sites']): return
    spectra(length,index,mesh,step)
    with np.load(out/'SPECTRA.npz') as z: k,e,mu=z['k'],z['e'],float(z['mu'])
    sites=case['readout_sites_native_1based'];nd=definition['dimension'];start=time.monotonic()
    derivatives=[derivative(out/f'coeff/dsite_{site:04}.bin') for site in sites]
    values=[None]*len(k);imaginaries=[0.]*len(k);energies=[0.]*len(k);counts=[0.]*len(k)
    # The PBS shell capped the old per-point readers. Here only one thread reads RVR.
    workers=int(os.environ.get('K_WORKERS',os.environ.get('V1_K_WORKERS','1')))
    for group in kz_groups(k):
        missing=[int(ki) for ki in group if not valid_npz(out/f'k_{ki:04}/TORQUE.npz',['torque','imaginary','count','F'])]
        if missing:
            tick=time.monotonic()
            densities=np.empty((len(missing),nd,nd),dtype=np.complex128)
            point_stats=[]
            for slot,ki in enumerate(missing):
                with np.load(out/f'k_{ki:04}/eigensystem.npz') as z: ev,vec=z['e'],z['v']
                x=(ev-mu)/SIGMA
                occ=.5*erfc(x)
                densities[slot]=(vec*occ)@vec.conj().T
                point_stats.append(dict(count=occ.sum(),F=np.sum(ev*occ-SIGMA*np.exp(-x*x)/(2*np.sqrt(np.pi)))))
                del vec
            ts,im,measure=contract_batch(WORK/'RVR/B_general',WORK/'basis/pure_H0.au',k[missing,:3],
                                         densities,derivatives,workers=workers)
            if length==2 and index in (0,1) and 0 in missing:
                slot=missing.index(0)
                site=sites[0];directory=out/f'k_0000/dsite_{site:04}'
                dh=load_h(native(out/f'coeff/dsite_{site:04}.bin',k[0],directory),nd)
                direct=-np.einsum('ij,ji->',densities[slot],dh);del dh
                error=float(abs(ts[slot,0]-direct))
                if error>1e-6:raise ValueError(f'Actual material adjoint/direct derivative mismatch: {error}')
                save(out/'ADJOINT_MATERIAL_CHECK.json',dict(status='PASS',site=site,k_index=0,error_eV=error))
            del densities
            for slot,ki in enumerate(missing):
                part=out/f'k_{ki:04}/TORQUE.npz'
                atomic_npz(part,torque=ts[slot],imaginary=im[slot],**point_stats[slot])
            print('TORQUE_KZ_BATCH',json.dumps(dict(L=length,case=index,kz=float(k[missing[0],2]),
                  k_indices=missing,workers=min(workers,len(missing)),seconds=time.monotonic()-tick,**measure)),flush=True)
        for ki in group:
            with np.load(out/f'k_{ki:04}/TORQUE.npz') as z:
                values[ki]=z['torque'];imaginaries[ki]=float(z['imaginary'])
                counts[ki]=float(z['count']);energies[ki]=float(z['F'])
            print('TORQUE',length,index,ki,flush=True)
    result=dict(status='COMPLETE',L=length,case=case,sites=sites,
                torque_eV=(k[:,3]@np.array(values)).tolist(),F_eV=float(k[:,3]@energies),
                electron_error=float(abs(k[:,3]@counts-definition['electrons'])),
                imaginary_gradient_eV=max(imaginaries),seconds=time.monotonic()-start,
                mesh=[mesh,mesh,8],derivative_step_rad=step)
    quality_check(result['electron_error']<=1e-8, f"Electron count error: {result['electron_error']}")
    save(out/'TORQUE.json',result)
    # Cleanup only after the entire case is published, never while it is resumable.
    for ki in range(len(k)):
        for site in sites:(out/f'k_{ki:04}/dsite_{site:04}/H.k000001.h').unlink(missing_ok=True)
        if index!=0:
            (out/f'k_{ki:04}/value/H.k000001.h').unlink(missing_ok=True)
            (out/f'k_{ki:04}/eigensystem.npz').unlink(missing_ok=True)


def collect(length: int,mesh: int,step: float) -> None:
    check_regression();definition,_,out=case_root(length,0,mesh,step);out=out.parent
    cases=definition['angle_configurations']
    for i,case in enumerate(cases):bind_case(length,definition,case,out/f'case_{i:02}',mesh,step)
    records=[json.loads((out/f'case_{i:02}/TORQUE.json').read_text()) for i in range(len(cases))]
    def value(perturb,angle,target):
        i=next(i for i,c in enumerate(cases) if c['perturb_site_native_1based']==perturb and c['angle_deg']==angle)
        row=records[i];return row['torque_eV'][row['sites'].index(target)]
    rows=[]
    for shell in definition['selected_shells']:
        representative=min(p['site_native_1based'] for p in shell['pairs'])
        for pair in shell['pairs']:
            target=pair['site_native_1based'];sign=pair['neel_sign'];jj=[]
            for delta in (.25,.5,1.):
                jj.append(-1000*(value(1,delta,target)-value(1,-delta,target))/(2*np.deg2rad(delta)*sign))
            small=(4*jj[0]-jj[1])/3;large=(4*jj[1]-jj[2])/3;reverse=None
            if target==representative:
                back=[-1000*(value(target,d,1)-value(target,-d,1))/(2*np.deg2rad(d)*sign) for d in (.25,.5)]
                reverse=(4*back[0]-back[1])/3
            rows.append(dict(shell=shell['label'],target=target,distance_A=shell['distance_A'],
                 J_meV=small,step_error_meV=abs(small-large),reverse_J_meV=reverse,
                 reciprocity_error_meV=None if reverse is None else abs(small-reverse),steps_meV=jj,
                 image_multiplicity=pair['folded_shell_multiplicity']))
    shells=[]
    for shell in definition['selected_shells']:
        values=[r['J_meV'] for r in rows if r['shell']==shell['label']]
        shells.append(dict(shell=shell['label'],mean_meV=float(np.mean(values)),spread_meV=float(np.ptp(values))))
    for name,data in [('J_pairs.csv',rows),('J_shells.csv',shells)]:
        with (out/name).open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(data[0]));writer.writeheader();writer.writerows(data)
    save(out/'RESULT.json',dict(status='COMPLETE',L=length,mesh=[mesh,mesh,8],pairs=rows,shells=shells,
         stationarity_meV=max(abs(np.array(records[0]['torque_eV'])))*1000,
         convention='unit spins; E=sum_i<j J_ij e_i dot e_j; positive J AFM',
         interpretation='periodic sums; no division by image multiplicity',
         k_mesh_converged=False,physical_DFT_accuracy_validated=False))


def main() -> None:
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['regression','spectrum','torque','collect'])
    p.add_argument('--L',type=int,choices=[2,4,6,8,10]);p.add_argument('--case',type=int,default=0)
    p.add_argument('--mesh',type=int);p.add_argument('--step',type=float,default=.001)
    a=p.parse_args();compute_node()
    if a.stage=='regression': regression_assemble();return
    if a.L is None or a.step<=0 or not np.isfinite(a.step): raise ValueError('Valid L and inner step required')
    mesh=a.mesh if a.mesh is not None else (8+a.L-1)//a.L
    if mesh<1: raise ValueError('Mesh must be positive')
    if a.stage=='collect': collect(a.L,mesh,a.step)
    elif a.stage=='spectrum': spectra(a.L,a.case,mesh,a.step)
    else: torque(a.L,a.case,mesh,a.step)


if __name__=='__main__': main()
