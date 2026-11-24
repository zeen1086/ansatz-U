"""Eight kz workers: complete trial functional, adjoint J12 or two-site K4."""
from pathlib import Path
import os,sys,json,time,subprocess,struct
from itertools import product
import numpy as np
from scipy.optimize import brentq
from density_map import NG,ws_weights,fold_unitary,fold_density,fold_project,cells_to_grid,grid_to_cells
from complete_readout import static,operators,metric,fermi,loewner,hubbard,trace,RY,SIGMA
from native_input import input_file,write_density,read_density
from common import ROOT, REF, CONFIG, WORK, INTERFACE, SOURCE, compute_node, load_h, quality_check, check_regression, shared_bundle
from common import atomic_npz, valid_npz, valid_json, bind_inputs, archive_incomplete, input_digest
from run_material import native_batch, binding, spin_structure, global_spin_rotation, coefficient_identity, valid_coefficients, prediction_identity
from torque_contract import derivative, contract_batch
from kz_production import stage_slice, layout
from spin_coordinates import curvature_cases, derivative_files
ENERGY=CONFIG['complete_energy']
R=Path(ENERGY['work']);MIN=Path(ENERGY['reference_run'])
SHARED=Path(ENERGY['shared_work'])
SIZES=ENERGY['sizes']

def case_definition(L):
    if 'curvature' not in ENERGY:
        return json.loads((ROOT/f'inputs/L{L:02}/CASE.json').read_text())
    spec = ENERGY['curvature']
    if not isinstance(spec, dict) or set(spec)-{'sites', 'outer_steps_rad', 'inner_step_rad'}:
        raise ValueError('Invalid curvature specification')
    steps = spec.get('outer_steps_rad', np.deg2rad([.25, .5, 1.]).tolist())
    cases = curvature_cases(spec['sites'], steps)
    inner = spec.get('inner_step_rad', .001)
    if not np.isfinite(inner) or inner <= 0 or max(steps)+2*inner >= np.pi:
        raise ValueError('Invalid inner derivative step')
    if max(spec['sites']) > L*L:
        raise ValueError('Curvature sites exceed the supercell')
    return dict(kind='K4', sites=spec['sites'], outer_steps_rad=steps,
                inner_step_rad=inner, angle_configurations=cases)


def readout_sites(definition, case):
    if definition.get('kind') == 'K4':
        return definition['sites']
    targets = [p['site_native_1based'] for p in definition['selected_shells'][0]['pairs']]
    return [1] if case == 0 or case >= 7 else [1, *targets]

def frozen_metadata():
    sys.path.insert(0,str(SOURCE))
    from local_vfunc import metadata
    meta=metadata(INTERFACE/'V_func_tmp')
    if meta['model_id']!=CONFIG['model_id']:
        raise ValueError('Frozen V1 identity changed')
    return meta

from common import save

ATTEMPT = None
_WORKER_LOCK = None

def start_attempt(z):
    import fcntl,re,uuid
    global ATTEMPT,_WORKER_LOCK
    run_id=os.environ.get('ANSATZ_RUN_ID','')
    R.mkdir(parents=True,exist_ok=True)
    if not run_id:
        key=os.environ.get('PBS_JOBID') or f'local-{os.getppid()}'
        key=key.replace('/', '_')
        registry=R/'run_ids';registry.mkdir(exist_ok=True)
        with (registry/key).open('a+') as entry:
            fcntl.flock(entry,fcntl.LOCK_EX)
            entry.seek(0);run_id=entry.read().strip()
            if not run_id:
                run_id=uuid.uuid4().hex
                entry.write(run_id+'\n');entry.flush();os.fsync(entry.fileno())
        print(f'NOTE: using automatic ANSATZ_RUN_ID={run_id} for {key}',flush=True)
    if not re.fullmatch(r'[A-Za-z0-9_-]{1,100}',run_id):
        raise ValueError('ANSATZ_RUN_ID must contain 1..100 letters, digits, underscores or hyphens')
    R.mkdir(parents=True,exist_ok=True)
    locks=R/'worker_locks';locks.mkdir(exist_ok=True)
    lock=(locks/f'{z}.lock').open('a')
    try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:
        lock.close();raise RuntimeError(f'Worker {z} is already running in {R}')
    _WORKER_LOCK=lock
    candidate=R/'attempts'/run_id;candidate.mkdir(parents=True,exist_ok=True)
    if (candidate/f'worker_{z}.started').exists():
        print(f'NOTE: reusing worker {z} start marker in {candidate}',flush=True)
    try:
        with (candidate/f'worker_{z}.started').open('w') as f:f.write(str(os.getpid()))
    except BaseException:
        lock.close();raise
    ATTEMPT=candidate
    barrier('startup',z)


def wait(path):
    if ATTEMPT is None:raise RuntimeError('Run attempt not initialized')
    noted=False
    while True:
        if not noted and (ATTEMPT/'FAILURE.json').exists():
            print(f'NOTE: worker failure recorded in {ATTEMPT}/FAILURE.json; continuing to wait',flush=True);noted=True
        if path.exists():return
        time.sleep(2)


def barrier(name,z):
    if ATTEMPT is None:raise RuntimeError('Run attempt not initialized')
    p=ATTEMPT/'barriers'/name;save(p/f'{z}.json',dict(done=True))
    for i in range(8):wait(p/f'{i}.json')

def points(L,z):return [(8*(2*x+y)+z,np.array([x/2,y/2,z/8])) for x,y in product(range(2),repeat=2)]
def folder(L,case):return R/f'L{L:02}/case_{case:02}'
def coefficient_case(L,case):
    source=ENERGY.get('reuse_cases',{}).get(str(L))
    return Path(source)/f'case_{case:02}' if source else R/f'coefficients/L{L:02}/case_{case:02}'
def cache_identity():
    # Metadata and declared model IDs bind inputs without rereading large payloads.
    files = {'bloch': WORK/'basis/bloch.npy', 'background': WORK/'basis/pure_H0.au',
             'projection': WORK/'RVR/PROJECTION.json', 'rvr': WORK/'RVR/B_general', 'model': INTERFACE/'V_func_tmp/V_FUNC.json'}
    for folder, label in [(REF, 'wannier'), (MIN/'results', 'reference_results'),
                          (MIN/'inputs/reference_actual', 'reference_density')]:
        for path in sorted(folder.iterdir()):
            if path.is_file(): files[label+'/'+path.name] = path
    for name in ('complete_energy.py','complete_readout.py','density_map.py','cal_Hsk.f90','common.py','run_material.py','torque_contract.py','native_input.py','kz_production.py'):
        files['code/'+name] = Path(__file__).with_name(name)
    files['code/spin_coordinates.py'] = SOURCE/'spin_coordinates.py'
    files['code/summarize_energy.py'] = Path(__file__).with_name('summarize_energy.py')
    files['code/frozen_io.f90']=SOURCE/'6_RVR/frozen_io.f90'
    return dict(version=1, files={name:input_digest(path) for name,path in files.items()},
                model=CONFIG['model_id'],
                grid=[40,40,36], orbitals=30, mesh=[2,2,8], sigma=SIGMA)


def bind_cache(identity):
    # Adopt unbound reference data, but never override an existing identity.
    import fcntl
    with (SHARED/'INPUTS.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        path=SHARED/'INPUTS.json'
        if path.exists():
            if json.loads(path.read_text()) != identity:
                raise ValueError('Shared inputs changed; use a new complete_energy.shared_work directory')
        else:
            if (SHARED/'reference').exists() or (SHARED/'cache').exists():
                print(f'NOTE: adopting shared cache without INPUTS.json: {SHARED}')
            save(path,identity)


def stage_bloch(z):
    src=np.load(WORK/'basis/bloch.npy',mmap_mode='r')
    if src.shape != (8,8,8,NG,30) or src.dtype != np.dtype('complex128'):
        raise ValueError('Unexpected common Bloch grid or dtype')
    return src[:,:,z]


def cached_operators(L,k,identity,st):
    key='_'.join(float(v).hex() for v in k)
    path=SHARED/'cache'/f'L{L:02}'/('operators_'+key)
    def build(out):
        h0,aa=operators(L,k,*st)
        G=fold_unitary(L,k)
        native_h=G.conj().T@h0@G
        native_a=aa@G
        np.save(out/'h0.npy',h0);np.save(out/'projectors.npy',aa)
        np.save(out/'native_projectors.npy',native_a)
        # Fortran stream: explicit model/basis, dimensions, k, then native-order matrices.
        basis=json.loads((WORK/'RVR/PROJECTION.json').read_text())['basis_id']
        with (out/'native.bin').open('wb') as f:
            f.write(b'V1_BG_CACHE_V1  '.ljust(16))
            f.write(struct.pack('<2i',L,60*L*L));f.write(np.asarray(k,dtype='<f8').tobytes())
            f.write(CONFIG['model_id'].encode());f.write(basis.encode())
            f.write(np.asarray(native_h,dtype='<c16').tobytes(order='F'))
            for a in native_a:f.write(np.asarray(a,dtype='<c16').tobytes(order='F'))
    return shared_bundle(path,dict(inputs=identity,L=L,k=list(map(float,k)),kind='operators'),build)


def load_operators(path):
    return np.load(path/'h0.npy',mmap_mode='r'),np.load(path/'projectors.npy',mmap_mode='r')


def reference(z,slab,dv,identity,st,ueff):
    out=SHARED/'reference';out.mkdir(exist_ok=True)
    if valid_npz(out/f'z{z}.npz',{'density':(4,NG,4),'N':(4,10,10),'count':()}):return
    data=np.load(REF/'data.npz');i0=data['names'].tolist().index('Neel_anchor');mu=np.load(MIN/'results/BANDS.npz')['mu'][1,0]
    density=np.zeros((4,NG,4));N=np.zeros((4,10,10),complex);count=0.
    for ki,k in enumerate(data['K']):
        if abs(k[2]-z/8)>1e-10:continue
        a=np.load(REF/f'k_{ki+1:04}.npz');bits=a['fold_bits'];order=[next(i for i,b in enumerate(bits) if np.array_equal(b,c)) for c in product(range(2),repeat=2)]
        ix=np.concatenate([np.arange(i*60,(i+1)*60) for i in order]);h=a['H_B_eV'][i0][np.ix_(ix,ix)];e,u=np.linalg.eigh(h);f,_,_=fermi(e,mu);p=(u*f)@u.conj().T
        q=(np.array(list(product(range(2),repeat=2)))+k[:2])/2;psi=np.array([slab[tuple(np.rint(v*8).astype(int)%8)] for v in q])
        if ki==0:
            fields=np.load(MIN/'results/FIELDS.npz')['DFT'][0];nin=np.load(MIN/'results/FIELDS.npz')['NDFT'][0][[0,2,1,3]]
            h0,aa=load_operators(cached_operators(2,k,identity,st));_,vu=hubbard(nin,ueff)
            reconstructed=h0+fold_project(psi,grid_to_cells(fields,2),2)+sum(aa[j].conj().T@vu[j]@aa[j] for j in range(4))
            from density_map import Projection,indices,spin_matrix
            old=Projection.__new__(Projection);old.ix=indices(data['grid']);old.bits=bits;old.block=30;old.w=np.load(REF/f'wr_{ki+1:04}.npy',mmap_mode='r')
            direct=spin_matrix(old.general(fields.T))[np.ix_(ix,ix)]
            local=fold_project(psi,grid_to_cells(fields,2),2)
            error=float(np.max(abs(local-direct)))
            h0error=float(np.max(abs(h0-a['H0_eV'][np.ix_(ix,ix)])))
            delta=reconstructed-h;shift=float(np.trace(delta).real/240)
            save(SHARED/'reference/OPERATOR_PARITY.json',dict(local_projection_error_eV=error,H0_roundtrip_error_eV=h0error,DFT_H_uniform_gauge_shift_eV=shift,DFT_H_residual_after_uniform_alignment_eV=float(np.max(abs(delta-shift*np.eye(240))))))
            quality_check(max(error,h0error)<1e-6, f'Density map/native projection residuals: {error}, {h0error}')
        density+=fold_density(psi,p,2)/128/dv;count+=f.sum()/128
        for j,site in enumerate([0,2,1,3]):
            aa=np.zeros((10,240),complex);aa[:5,::2]=data['A'][ki,site];aa[5:,1::2]=data['A'][ki,site];aa=aa[:,ix]
            N[j]+=aa@p@aa.conj().T/128
    atomic_npz(out/f'z{z}.npz',density=density,N=N,count=count)

def background(z,dv):
    barrier('reference',z)
    if z==0 and not valid_npz(SHARED/'reference/background.npz',{'density':(4,NG,4),'N':(4,10,10)}):
        rows=[np.load(SHARED/f'reference/z{i}.npz') for i in range(8)];act=sum(v['density'] for v in rows);n=sum(v['N'] for v in rows)
        full=grid_to_cells(read_density(MIN/'inputs/reference_actual/DENSITY.bin',2),2);nt=np.load(MIN/'results/FIELDS.npz')['NDFT'][0][[0,2,1,3]]
        if any(not np.isfinite(v).all() for v in (act,n,full,nt)):raise ValueError('Nonfinite reference density or occupation')
        quality_check(abs(act[:,:,0].sum()*dv-84)<1e-7, f'Reference active electron count: {act[:,:,0].sum()*dv}; expected 84')
        atomic_npz(SHARED/'reference/background.npz',density=full-act,N=nt-n)
    barrier('background',z)

def wave_cache(L,z,identity,st,slab):
    paths={};ops={};checks=[];cells=np.array(list(product(range(L),repeat=2)))
    for ki,k in points(L,z):
        ops[ki]=cached_operators(L,k,identity,st)
        def build(out):
            import shutil
            required=16*(L*L*NG*30+2*(60*L*L)**2)
            if shutil.disk_usage(out).free < required*1.05:
                raise OSError(f'Insufficient shared space for L={L}, k={ki}: need at least {required} bytes')
            source_slab=np.ascontiguousarray(slab).reshape(64,NG,30)
            psi=np.lib.format.open_memmap(out/'wave.npy',mode='w+',dtype=complex,shape=(L*L,NG,30))
            for j,q in enumerate((cells+k[:2])/L):
                weights=ws_weights(np.r_[q,z/8])[:,:,z];ix=np.flatnonzero(abs(weights.reshape(-1))>1e-13)
                if len(ix)==1:psi[j]=source_slab[ix[0]]*weights.reshape(-1)[ix[0]]
                else:psi[j]=np.einsum('q,qgo->go',weights.reshape(-1)[ix],source_slab[ix],optimize=True)
            psi.flush();S,T,err=metric(psi)
            np.save(out/'overlap.npy',S);np.save(out/'transform.npy',T)
            save(out/'metric.json',dict(gram_max=err))
        path=SHARED/'cache'/f'L{L:02}'/f'k{ki:04}'
        path=shared_bundle(path,dict(inputs=identity,L=L,k=k.tolist(),kind='wave_metric'),build)
        paths[ki]=path/'wave.npy'
        checks.append(dict(k=ki,gram_max=json.loads((path/'metric.json').read_text())['gram_max']))
    save(SHARED/f'L{L:02}/metrics/z{z}.json',dict(points=checks))
    return paths,ops

def value_offset(path: Path, length: int, ueff: float) -> float:
    with path.open('rb') as stream:
        header = stream.read(128)
    if len(header) != 128:
        raise ValueError(f'Truncated coefficient header: {path}')
    magic, *fields, model_id, coefficient_ueff, offset = struct.unpack('<16s8i64s2d', header)
    version, ncoeff, cell, sites, derivative = fields[:5]
    if (magic != b'LANIO2_COEFF_V2 ' or version != 1 or ncoeff < 1 or cell != length
            or sites != length**2 or derivative != 0 or model_id != CONFIG['model_id'].encode()
            or not np.isfinite(offset)):
        raise ValueError(f'Invalid value coefficient header: {path}')
    if coefficient_ueff != ueff:
        raise ValueError(f'Coefficient Ueff {coefficient_ueff} differs from model Ueff {ueff}: {path}')
    return offset


def spectrum(L,case,z,paths,ops,rvr,ueff):
    p=folder(L,case);p.mkdir(parents=True,exist_ok=True)
    source=coefficient_case(L,case)
    batch=points(L,z)
    coeff=source/'coeff/value.bin';offset=value_offset(coeff,L,ueff)
    native_batch(coeff,np.array([k for _,k in batch]),[p/f'H_{ki:04}' for ki,_ in batch],p,60*L*L,rvr=rvr,background_files=[ops[ki]/'native.bin' for ki,_ in batch])
    for ki,k in batch:
        out=p/f'k_{ki:04}.npz'
        if valid_npz(out,{'e':(60*L*L,),'Q':(60*L*L,60*L*L),'H':(60*L*L,60*L*L)}):continue
        G=fold_unitary(L,k);T=np.load(paths[ki].parent/'transform.npy',mmap_mode='r');err=json.loads((paths[ki].parent/'metric.json').read_text())['gram_max']
        S=np.load(paths[ki].parent/'overlap.npy',mmap_mode='r')
        hpath=p/f'H_{ki:04}/H.k000001.h'
        if not hpath.exists():
            target=p/f'H_{ki:04}';native_batch(coeff,k[None],[target],p,60*L*L,rvr=rvr,background_files=[ops[ki]/'native.bin']);hpath=target/'H.k000001.h'
        h=G@load_h(hpath,60*L*L)@G.conj().T
        # Replace the native cI with the scalar operator cS in the folded basis.
        if offset:
            h+=offset*S
            h.flat[::len(h)+1]-=offset
        checkpoint=source/f'k_{ki:04}/eigensystem.npz'
        if err<1e-10 and valid_npz(checkpoint,{'e':(60*L*L,),'v':(60*L*L,60*L*L)}):
            with np.load(checkpoint) as b:
                e=b['e'];Q=T@(G@b['v'])
            if np.max(abs(h@Q-(S@Q)*e))>=1e-8:
                e,u=np.linalg.eigh(T@h@T);Q=T@u
        else:
            e,u=np.linalg.eigh(T@h@T);Q=T@u
        atomic_npz(out,e=e,Q=Q,H=h)
        print('SPECTRUM',L,case,ki,'metric',err,flush=True)
    barrier(f'spectrum_L{L}_c{case}',z)
    if z==0:
        e=np.array([np.load(p/f'k_{i:04}.npz')['e'] for i in range(32)]);target=21*L*L
        mu=brentq(lambda m:np.mean([fermi(a,m)[0].sum() for a in e])-target,e.min()-20*SIGMA,e.max()+20*SIGMA,xtol=1e-12)
        save(p/'MU.json',dict(mu=mu,electrons=target))
    barrier(f'mu_L{L}_c{case}',z)

def density_stage(L,case,z,paths,dv,ops):
    p=folder(L,case);out=p/f'density_z{z}.npz'
    if valid_npz(out,{'density':(L*L,NG,4),'N':(L*L,10,10),'bare':(),'F':(),'band':(),'count':()}):return
    start=time.monotonic();mu=json.loads((p/'MU.json').read_text())['mu'];acc=np.zeros((L*L,NG,4));N=np.zeros((L*L,10,10),complex);bare=0.;F=0.;band=0.;count=0.
    for ki,k in points(L,z):
        b=np.load(p/f'k_{ki:04}.npz');e,Q=b['e'],b['Q'];f,_,free=fermi(e,mu);P=(Q*f)@Q.conj().T;h0,aa=load_operators(ops[ki])
        psi=np.load(paths[ki],mmap_mode='r');acc+=fold_density(psi,P,L)/32/dv
        N+=np.array([a@P@a.conj().T for a in aa])/32;bare+=trace(P,h0)/32;F+=free/32;band+=float(e@f)/32;count+=float(f.sum())/32
        print('DENSITY_POINT',L,case,ki,flush=True)
    atomic_npz(out,density=acc,N=N,bare=bare,F=F,band=band,count=count,seconds=time.monotonic()-start)
    print('DENSITY_CASE',L,case,z,time.monotonic()-start,flush=True)

def reference_sites(length):
    x, y = np.indices((length, length))
    return (2*(x % 2) + y % 2).ravel()


def total_occupation(active, reference_n, length, case):
    if active.shape != (length**2,10,10) or reference_n.shape != (4,10,10):
        raise ValueError('Expected one active Ni shell per cell and four reference shells')
    if 'global_spin_rotation_deg' in case:
        x, y, z, w = global_spin_rotation(case).as_quat()
        # U = w I - i (x,y,z)·sigma, in the up-d1..5, down-d1..5 ordering.
        spinor = np.array([[w-1j*z, -y-1j*x], [y-1j*x, w+1j*z]])
        rotation = np.kron(spinor, np.eye(5))
        reference_n = rotation @ reference_n @ rotation.conj().T
    return active + reference_n[reference_sites(length)]


def assemble_density(L,case,dv):
    p=folder(L,case)
    if (valid_npz(p/'DENSITY.npz',{'N':(L*L,10,10)}) and valid_json(p/'ELECTRONIC.json',['F','band','bare','count'])
            and (p/'DENSITY.bin').is_file() and (p/'DENSITY.bin').stat().st_size==28+32*NG*L*L):return
    case_input=case_definition(L)['angle_configurations'][case]
    rows=[np.load(p/f'density_z{i}.npz') for i in range(8)]
    with np.load(SHARED/'reference/background.npz') as bg:
        reference_density=bg['density']
        if 'global_spin_rotation_deg' in case_input:
            reference_density[...,1:] = reference_density[...,1:] @ global_spin_rotation(case_input).as_matrix().T
        rho=sum(r['density'] for r in rows)+reference_density[reference_sites(L)]
        N=total_occupation(sum(r['N'] for r in rows),bg['N'],L,case_input)
    charge=float(rho[:,:,0].sum()*dv)
    if not np.isfinite(rho).all() or not np.isfinite(N).all():raise ValueError('Nonfinite density or occupation')
    quality_check(abs(charge-41*L*L)<1e-6, f'L={L}, case={case}: charge={charge}; expected {41*L*L}')
    atomic_npz(p/'DENSITY.npz',N=N);write_density(p/'DENSITY.bin',cells_to_grid(rho,L),L)
    save(p/'ELECTRONIC.json',dict(charge=charge,charge_error=abs(charge-41*L*L),**{key:float(sum(r[key] for r in rows)) for key in ['bare','F','band','count']}))

def native(L,ueff,case_count=11):
    import tempfile
    out=R/f'L{L:02}/native';out.parent.mkdir(parents=True,exist_ok=True)
    identity={str(folder(L,i)/'DENSITY.bin'):input_digest(folder(L,i)/'DENSITY.bin') for i in range(case_count)}
    identity['binary']=input_digest(Path(ENERGY['native_binary']))
    names=[f'case_{i:02}_{suffix}' for i in range(case_count) for suffix in ('ENERGY.dat','TRUE.bin')]
    complete=False
    try:
        report=json.loads((out/'COMPLETE.json').read_text())
        if report['inputs']!=identity:raise ValueError(f'Native density inputs changed: {out}')
        complete=set(report['files'])==set(names) and all(input_digest(out/name)==report['files'][name] for name in names)
    except (OSError,KeyError,json.JSONDecodeError):pass
    if not complete:
        archive_incomplete(out)
        with tempfile.TemporaryDirectory(prefix='.native_',dir=out.parent) as tmp:
            stage=Path(tmp)/'native';stage.mkdir()
            input_file(stage,L);(stage/'input_list.txt').write_text(f'{case_count}\n'+''.join(f'case_{i:02} "{folder(L,i)}/DENSITY.bin" 0 0 1\n' for i in range(case_count)))
            try:
                with (stage/'native.out').open('w') as log:
                    subprocess.run([str(Path(ENERGY['native_binary'])),'grid.in'],cwd=stage,stdout=log,stderr=subprocess.STDOUT,check=True)
            except BaseException:
                stage.rename(out)
                raise
            if 'GGA_BRANCH_REEXPORT_COMPLETE' not in (stage/'native.out').read_text():raise RuntimeError('Incomplete native evaluation')
            for i in range(case_count):
                ee=np.loadtxt(stage/f'case_{i:02}_ENERGY.dat')
                if ee.size<7 or not np.isfinite(ee).all():raise ValueError('Invalid native energy output')
                field=read_density(stage/f'case_{i:02}_TRUE.bin',L)
                if not np.isfinite(field).all():raise ValueError('Nonfinite native potential output')
            save(stage/'COMPLETE.json',dict(inputs=identity,files={name:input_digest(stage/name) for name in names}))
            stage.rename(out)
    # Recreate these cheap summaries even if interruption followed native publication.
    for case in range(case_count):
        p=folder(L,case);ee=np.loadtxt(out/f'case_{case:02}_ENERGY.dat');n=np.load(p/'DENSITY.npz')['N'];eu,_=hubbard(n,ueff);b=json.loads((p/'ELECTRONIC.json').read_text());terms=dict(band_free=b['F'],minus_active_potential=-b['band']+b['bare'],EH=ee[0]*RY,Exc=ee[1]*RY,minus_core=-ee[2]*RY,ionic_density=ee[6]*RY,EU=eu)
        save(p/'ENERGY.json',dict(terms=terms,energy=sum(terms.values()),Ueff_eV=ueff,constant_reference_excluded=True))

def potential_stage(L,case,z,paths,ops,ueff):
    p=folder(L,case);mu=json.loads((p/'MU.json').read_text())['mu'];field=grid_to_cells(read_density(R/f'L{L:02}/native/case_{case:02}_TRUE.bin',L),L)*RY
    nt=np.load(p/'DENSITY.npz')['N'];_,vu=hubbard(nt,ueff);numerator=0.;denominator=0.
    for ki,k in points(L,z):
        dest=p/f'adjoint_{ki:04}.npz'
        if valid_npz(dest,{'Ae':(60*L*L,60*L*L),'numerator':(),'denominator':()}):
            a=np.load(dest);numerator+=float(a['numerator']);denominator+=float(a['denominator']);continue
        b=np.load(p/f'k_{ki:04}.npz');e,Q,H=b['e'],b['Q'],b['H'];h0,aa=load_operators(ops[ki]);psi=np.load(paths[ki],mmap_mode='r')
        A=H-h0-fold_project(psi,field,L)-sum(a.conj().T@v@a for a,v in zip(aa,vu));ae=Q.conj().T@A@Q;f,fp,_=fermi(e,mu)
        num=float(np.dot(fp,np.diag(ae).real))/32;den=float(fp.sum())/32;atomic_npz(dest,Ae=ae,numerator=num,denominator=den);numerator+=num;denominator+=den
        print('NATIVE_PROJECTED',L,case,ki,flush=True)
    save(p/f'alpha_z{z}.json',dict(numerator=numerator,denominator=denominator))

def torque_stage(L,z,rvr,ops,ueff):
    definition=case_definition(L); tangent=definition.get('kind')=='K4'
    count=len(definition['angle_configurations'])
    for start in range(0,count,4):
        chunk=list(range(start,min(start+4,count)))
        entries=[];densities=[];ds=[]
        for case in chunk:
            p=folder(L,case);mu=json.loads((p/'MU.json').read_text())['mu'];al=[json.loads((p/f'alpha_z{i}.json').read_text()) for i in range(8)];alpha=sum(a['numerator'] for a in al)/sum(a['denominator'] for a in al)
            sites=readout_sites(definition,case)
            derivs=[derivative(coefficient_case(L,case)/'coeff'/name) for name in derivative_files(sites,tangent)]
            labels=[(s,a) for s in sites for a in (1,2)] if tangent else sites
            if any(d['Ueff'] != ueff for d in derivs):
                raise ValueError(f'Derivative coefficient Ueff differs from model Ueff {ueff}: {coefficient_case(L,case)}')
            for ki,k in points(L,z):
                dest=p/f'torque_k{ki:04}.npz'
                if valid_npz(dest,{'torque':(len(labels),),'sites':np.shape(labels),'imaginary':()}):continue
                b=np.load(p/f'k_{ki:04}.npz');e,Q=b['e'],b['Q'];ae=np.load(p/f'adjoint_{ki:04}.npz')['Ae'];ae.flat[::len(e)+1]-=alpha;f,fp,_=fermi(e,mu);Z=Q@(loewner(e,f,fp)*ae)@Q.conj().T;G=fold_unitary(L,k)
                densities.append(-G.conj().T@Z@G);entries.append((case,ki,k,dest,labels));ds.append(derivs)
        if entries:
            vals,im,stats=contract_batch(rvr,WORK/'basis/pure_H0.au',np.array([a[2] for a in entries]),np.array(densities),ds[0],workers=1,derivatives_by_point=ds,projectors_by_point=[np.load(ops[a[1]]/'native_projectors.npy',mmap_mode='r') for a in entries])
            for row,v,imag in zip(entries,vals,im):atomic_npz(row[3],torque=v,sites=row[4],imaginary=imag)
            save(R/f'L{L:02}/torque_batch_{z}_{chunk[0]}.json',stats)

def collect(L):
    import csv
    d=case_definition(L);rec=[]
    for case in range(len(d['angle_configurations'])):
        rows=[]
        for i in range(32):
            with np.load(folder(L,case)/f'torque_k{i:04}.npz') as data:
                rows.append(dict(data))
        expected=[(s,a) for s in d['sites'] for a in (1,2)] if d.get('kind')=='K4' else readout_sites(d,case)
        if any(not np.array_equal(r['sites'],expected) or r['torque'].shape!=(len(expected),)
               or not np.isfinite(r['torque']).all() for r in rows):
            raise ValueError('Torque modes or values do not match the requested case')
        rec.append(dict(sites=rows[0]['sites'].tolist(),torque=np.mean([r['torque'] for r in rows],axis=0).tolist()))
    if d.get('kind')=='K4':
        from summarize_energy import curvature_result
        energies=[json.loads((folder(L,i)/'ENERGY.json').read_text())['energy'] for i in range(len(rec))]
        result=curvature_result(d,energies,[r['torque'] for r in rec])
        base=spin_structure(L,d['angle_configurations'][0])
        result.update(L=L,base_spins=base['chart_base_spins'].tolist(),tangent_basis=base['chart_basis'].tolist(),
                      mesh=[2,2,8],sigma_eV=SIGMA,active_electrons=21*L*L,total_electrons=41*L*L,
                      functional='complete trial-density Gaussian free energy',
                      complement='fixed reference density and occupation in the input frame',
                      other_spins='fixed',physical_validation=False,
                      interpretation='periodic supercell Hessian; no image-multiplicity division',
                      provenance=str(R/'INPUTS.json'))
        save(R/f'L{L:02}/RESULT.json',result)
        return
    shell=d['selected_shells'][0]
    def t(case,site):return rec[case]['torque'][rec[case]['sites'].index(site)]
    rows=[]
    for p in shell['pairs']:
        site=p['site_native_1based'];sgn=p['neel_sign'];jj=[-1000*(t(2*j+1,site)-t(2*j+2,site))/(2*np.deg2rad(a)*sgn) for j,a in enumerate([.25,.5,1.])];val=(4*jj[0]-jj[1])/3;back=None
        if site==2:
            rev=[-1000*(t(7+2*j,1)-t(8+2*j,1))/(2*np.deg2rad(a)*sgn) for j,a in enumerate([.25,.5])];back=(4*rev[0]-rev[1])/3
        rows.append(dict(target=site,J12_meV=val,steps_meV=jj,step_indicator_meV=abs(val-(4*jj[1]-jj[2])/3),reverse_J12_meV=back))
    save(R/f'L{L:02}/RESULT.json',dict(L=L,J12_meV=float(np.mean([a['J12_meV'] for a in rows])),spread_meV=float(np.ptp([a['J12_meV'] for a in rows])),pairs=rows,torques=rec,mesh=[2,2,8],inner_step_rad=.001,outer_steps_deg=[.25,.5,1.],density='Hamiltonian active plus frozen DFT excluded reference',metric='common Wannier overlap retained in generalized diagonalization',physical_validation=False,V1_unchanged=True))

def prepare_coefficients(L,z):
    from run_material import coefficients
    definition=case_definition(L);tangent=definition.get('kind')=='K4'
    step=definition['inner_step_rad'] if tangent else .001
    if not tangent:
        targets=[p['site_native_1based'] for p in definition['selected_shells'][0]['pairs']]
        assert sorted(targets)==sorted([2,L,L+1,L*(L-1)+1])
        expected=[(None,0),(1,.25),(1,-.25),(1,.5),(1,-.5),(1,1),(1,-1),(2,.25),(2,-.25),(2,.5),(2,-.5)]
        assert [(c['perturb_site_native_1based'],c['angle_deg']) for c in definition['angle_configurations']]==expected
    for i,case in enumerate(definition['angle_configurations']):
        if i%8!=z:continue
        out=coefficient_case(L,i)
        base=spin_structure(L,case)
        sites=readout_sites(definition,i)
        if str(L) in ENERGY.get('reuse_cases',{}):
            if not valid_coefficients(out/'coeff',coefficient_identity(base,sites,step,tangent=tangent)):
                raise ValueError(f'External coefficients incomplete or unbound: {out}')
        else:
            out.mkdir(parents=True,exist_ok=True)
            atomic_npz(out/'structure.npz',**base)
            coefficients(out/'structure.npz',out/'coeff',sites,step,tangent=tangent)
    barrier(f'coefficients_L{L}',z)

def main(z):
    compute_node()
    if not 0<=z<8 or not SIZES or len(set(SIZES))!=len(SIZES) or any(L not in (4,6,8) for L in SIZES):
        raise ValueError('Require kz worker 0..7 and distinct sizes from 4,6,8')
    if not (R/'INPUTS.json').is_file() and any(
            path.is_file() for L in SIZES for path in (R/f'L{L:02}').glob('case_*/*')):
        raise ValueError(f'Unbound complete-energy checkpoints have unknown energy gauge: {R}; use a new work directory')
    definitions={L:case_definition(L) for L in SIZES}
    R.mkdir(parents=True,exist_ok=True)
    SHARED.mkdir(parents=True,exist_ok=True)
    start_attempt(z)
    check_regression();binding()
    meta=frozen_metadata()
    assert len(meta['files'])+1==13
    ueff=meta['Ueff_eV']
    if not np.isfinite(ueff) or ueff < 0:
        raise ValueError(f'Invalid model Ueff: {ueff}')
    identity=cache_identity();identity['Ueff_eV']=ueff;bind_cache(identity)
    input_names=('structure.npz',) if 'curvature' in ENERGY else ('structure.npz','CASE.json')
    spin_files={str(ROOT/f'inputs/L{L:02}'/name):input_digest(ROOT/f'inputs/L{L:02}'/name) for L in SIZES for name in input_names}
    spin_files['native_binary']=input_digest(Path(ENERGY['native_binary']))
    spin_files['prediction_runtime']=prediction_identity()
    for root in ENERGY.get('reuse_cases',{}).values():
        for path in sorted(Path(root).glob('case_*/coeff/*')):
            if path.is_file():spin_files[str(path)]=input_digest(path)
    for path in sorted(Path(ENERGY['pseudo_dir']).glob('*.upf')):spin_files[str(path)]=input_digest(path)
    bind_inputs(R,dict(shared=identity,spin=spin_files,energy=ENERGY),ignored=('attempts','worker_locks','run_ids'))
    barrier('inputs',z)
    if z==0:(R/'COMPLETE.json').unlink(missing_ok=True)
    for L in SIZES:
        prepare_coefficients(L,z)
    dv=abs(np.linalg.det(np.array(meta['primitive_cell_A'])/.529177210903))/NG
    slab=stage_bloch(z);st=static(WORK/'basis/pure_H0.au');reference(z,slab,dv,identity,st,ueff);background(z,dv)
    cached=CONFIG.get('rvr_slices',{}).get(str(z))
    if cached:
        rvr=Path(cached)
    else:
        cache=shared_bundle(SHARED/'cache'/f'z{z}',dict(inputs=identity,kz=z/8,kind='rvr_slice'),
                            lambda out:stage_slice(WORK/'RVR/B_general',out/'RVR.bin',z))
        rvr=cache/'RVR.bin'
    prefix,sectors,_,size=layout(rvr)
    original,_,_,full_size=layout(WORK/'RVR/B_general')
    if sectors!=(z/8,) or size!=full_size or prefix[:32]!=original[:32] or prefix[36:]!=original[36:]:
        raise ValueError('Cached RVR does not match the full source and kz')
    save(SHARED/f'reused_RVR_z{z}.json',dict(path=str(rvr),bytes=size))
    for L in SIZES:
        paths,ops=wave_cache(L,z,identity,st,slab)
        count=len(definitions[L]['angle_configurations'])
        for case in range(count):
            spectrum(L,case,z,paths,ops,rvr,ueff);density_stage(L,case,z,paths,dv,ops);barrier(f'density_L{L}_c{case}',z)
            if z==0:assemble_density(L,case,dv)
        barrier(f'all_density_L{L}',z)
        if z==0:native(L,ueff,count)
        barrier(f'native_L{L}',z)
        for case in range(count):potential_stage(L,case,z,paths,ops,ueff)
        barrier(f'potentials_L{L}',z);torque_stage(L,z,rvr,ops,ueff);barrier(f'torques_L{L}',z)
        if z==0:collect(L)
    if z==0:save(R/'COMPLETE.json',dict(V1_files_verified=len(frozen_metadata()['files'])+1,sizes=SIZES,
                                     **({'observable':'K4'} if 'curvature' in ENERGY else {'J':'J12 only'})))

if __name__=='__main__':
    try:main(int(sys.argv[1]))
    except BaseException as e:
        if ATTEMPT is not None:save(ATTEMPT/'FAILURE.json',dict(worker=int(sys.argv[1]),error=repr(e)))
        raise
