"""Eight kz workers: frozen V1, Hamiltonian density, native functional, complete adjoint J12."""
from pathlib import Path
import os,sys,json,time,subprocess
from itertools import product
import numpy as np
from scipy.optimize import brentq
from density_map import NG,ws_weights,fold_unitary,fold_density,fold_project,cells_to_grid,grid_to_cells
from complete_readout import static,operators,metric,fermi,loewner,hubbard,trace,RY,SIGMA
from native_input import input_file,write_density,read_density
from common import ROOT, REF, CONFIG, WORK, INTERFACE, SOURCE, compute_node, load_h, require_regression
from run_material import native_batch, binding
from torque_contract import derivative, contract_batch
from kz_production import stage_slice, layout
ENERGY=CONFIG['complete_energy']
R=Path(ENERGY['work']);OLD=WORK;MIN=Path(ENERGY['reference_run'])
SIZES=ENERGY['sizes']

def frozen():
    sys.path.insert(0,str(SOURCE))
    from local_vfunc import verify
    meta=verify(INTERFACE/'V_func_tmp')
    if meta['model_id']!=CONFIG['model_id'] or meta['source_parameters_sha256']!=CONFIG['parameters_sha256']:
        raise ValueError('Frozen V1 identity changed')
    return len(meta['files'])+1

def save(p,value):
    p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix(p.suffix+'.partial');tmp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');tmp.replace(p)
def wait(path):
    start=time.monotonic()
    while not path.exists():
        if (R/'FAILURE.json').exists():raise RuntimeError('Another worker failed; inspect FAILURE.json')
        if time.monotonic()-start>20*3600:raise TimeoutError(str(path))
        time.sleep(2)
def barrier(name,z):
    p=R/'barriers'/name;save(p/f'{z}.json',dict(done=True))
    for i in range(8):wait(p/f'{i}.json')
def points(L,z):return [(8*(2*x+y)+z,np.array([x/2,y/2,z/8])) for x,y in product(range(2),repeat=2)]
def folder(L,case):return R/f'L{L:02}/case_{case:02}'
def oldcase(L,case):
    source=ENERGY.get('reuse_cases',{}).get(str(L))
    return Path(source)/f'case_{case:02}' if source else R/f'coefficients/L{L:02}/case_{case:02}'
def stage_bloch(z):
    cache=Path(os.environ.get('TMPDIR','/tmp'))/f'V1_J12_density_{os.environ["PBS_JOBID"].split(".")[0]}_{z}';cache.mkdir(parents=True,exist_ok=True)
    src=np.load(OLD/'basis/bloch.npy',mmap_mode='r');slab=np.empty((8,8,NG,30),complex)
    for x,y in product(range(8),repeat=2):slab[x,y]=src[x,y,z]
    return cache,slab

def reference(z,slab,dv):
    out=R/'reference';out.mkdir(exist_ok=True)
    if (out/f'z{z}.npz').exists():return
    data=np.load(REF/'data.npz');i0=data['names'].tolist().index('Neel_anchor');mu=np.load(MIN/'results/BANDS.npz')['mu'][1,0]
    density=np.zeros((4,NG,4));N=np.zeros((4,10,10),complex);count=0.
    for ki,k in enumerate(data['K']):
        if abs(k[2]-z/8)>1e-10:continue
        a=np.load(REF/f'k_{ki+1:04}.npz');bits=a['fold_bits'];order=[next(i for i,b in enumerate(bits) if np.array_equal(b,c)) for c in product(range(2),repeat=2)]
        ix=np.concatenate([np.arange(i*60,(i+1)*60) for i in order]);h=a['H_B_eV'][i0][np.ix_(ix,ix)];e,u=np.linalg.eigh(h);f,_,_=fermi(e,mu);p=(u*f)@u.conj().T
        q=(np.array(list(product(range(2),repeat=2)))+k[:2])/2;psi=np.array([slab[tuple(np.rint(v*8).astype(int)%8)] for v in q])
        if ki==0:
            fields=np.load(MIN/'results/FIELDS.npz')['DFT'][0];nin=np.load(MIN/'results/FIELDS.npz')['NDFT'][0][[0,2,1,3]]
            h0,aa=operators(2,k,*static(OLD/'basis/pure_H0.au'));_,vu=hubbard(nin)
            reconstructed=h0+fold_project(psi,grid_to_cells(fields,2),2)+sum(aa[j].conj().T@vu[j]@aa[j] for j in range(4))
            from reference_projection import Projection,indices,spin_matrix
            old=Projection.__new__(Projection);old.ix=indices(data['grid']);old.bits=bits;old.block=30;old.w=np.load(REF/f'wr_{ki+1:04}.npy',mmap_mode='r')
            direct=spin_matrix(old.general(fields.T))[np.ix_(ix,ix)]
            local=fold_project(psi,grid_to_cells(fields,2),2)
            error=float(np.max(abs(local-direct)))
            h0error=float(np.max(abs(h0-a['H0_eV'][np.ix_(ix,ix)])))
            delta=reconstructed-h;shift=float(np.trace(delta).real/240)
            save(R/'reference/OPERATOR_PARITY.json',dict(local_projection_error_eV=error,H0_roundtrip_error_eV=h0error,DFT_H_uniform_gauge_shift_eV=shift,DFT_H_residual_after_uniform_alignment_eV=float(np.max(abs(delta-shift*np.eye(240))))))
            assert max(error,h0error)<1e-8,('Density map/native projection mismatch',error,h0error)
        density+=fold_density(psi,p,2)/128/dv;count+=f.sum()/128
        for j,site in enumerate([0,2,1,3]):
            aa=np.zeros((10,240),complex);aa[:5,::2]=data['A'][ki,site];aa[5:,1::2]=data['A'][ki,site];aa=aa[:,ix]
            N[j]+=aa@p@aa.conj().T/128
    np.savez(out/f'z{z}.npz',density=density,N=N,count=count)

def background(z,dv):
    barrier('reference',z)
    if z==0 and not (R/'reference/background.npz').exists():
        rows=[np.load(R/f'reference/z{i}.npz') for i in range(8)];act=sum(v['density'] for v in rows);n=sum(v['N'] for v in rows)
        full=grid_to_cells(read_density(MIN/'inputs/reference_actual/DENSITY.bin',2),2);nt=np.load(MIN/'results/FIELDS.npz')['NDFT'][0][[0,2,1,3]]
        assert abs(act[:,:,0].sum()*dv-84)<1e-7
        np.savez(R/'reference/background.npz',density=full-act,N=nt-n)
    wait(R/'reference/background.npz')

def wave_cache(L,z,cache,slab):
    paths={};checks=[];cells=np.array(list(product(range(L),repeat=2)))
    for ki,k in points(L,z):
        path=cache/f'L{L}_k{ki}.npy';paths[ki]=path
        if path.exists():continue
        psi=np.lib.format.open_memmap(path,mode='w+',dtype=complex,shape=(L*L,NG,30))
        for j,q in enumerate((cells+k[:2])/L):
            weights=ws_weights(np.r_[q,z/8])[:,:,z];ix=np.flatnonzero(abs(weights.reshape(-1))>1e-13)
            if len(ix)==1:psi[j]=slab.reshape(64,NG,30)[ix[0]]*weights.reshape(-1)[ix[0]]
            else:psi[j]=np.einsum('q,qgo->go',weights.reshape(-1)[ix],slab.reshape(64,NG,30)[ix],optimize=True)
        psi.flush();S,T,err=metric(psi);np.savez(cache/f'L{L}_k{ki}_metric.npz',T=T,gram_max=err)
        checks.append(dict(k=ki,gram_max=err))
    save(R/f'L{L:02}/metrics/z{z}.json',dict(points=checks))
    return paths

def spectrum(L,case,z,paths,cache,rvr):
    p=folder(L,case);p.mkdir(parents=True,exist_ok=True)
    definition=json.loads((ROOT/f'inputs/L{L:02}/CASE.json').read_text());sites=definition['angle_configurations'][case]['readout_sites_native_1based'];source=oldcase(L,case)
    batch=points(L,z)
    native_batch(source/'coeff/value.bin',np.array([k for _,k in batch]),[p/f'H_{ki:04}' for ki,_ in batch],p,60*L*L,rvr=rvr)
    for ki,k in batch:
        out=p/f'k_{ki:04}.npz'
        if out.exists():continue
        G=fold_unitary(L,k);T=np.load(cache/f'L{L}_k{ki}_metric.npz')['T'];err=float(np.load(cache/f'L{L}_k{ki}_metric.npz')['gram_max'])
        hpath=p/f'H_{ki:04}/H.k000001.h'
        if not hpath.exists():
            target=p/f'H_{ki:04}';native_batch(source/'coeff/value.bin',k[None],[target],p,60*L*L,rvr=rvr);hpath=target/'H.k000001.h'
        h=G@load_h(hpath,60*L*L)@G.conj().T
        checkpoint=source/f'k_{ki:04}/eigensystem.npz'
        if err<1e-10 and checkpoint.exists():
            b=np.load(checkpoint);e=b['e'];Q=G@b['v']
            assert np.max(abs(h@Q-Q*e))<1e-8
        else:
            e,u=np.linalg.eigh(T@h@T);Q=T@u
        np.savez(out,e=e,Q=Q,H=h)
        print('SPECTRUM',L,case,ki,'metric',err,flush=True)
    barrier(f'spectrum_L{L}_c{case}',z)
    if z==0:
        e=np.array([np.load(p/f'k_{i:04}.npz')['e'] for i in range(32)]);target=21*L*L
        mu=brentq(lambda m:np.mean([fermi(a,m)[0].sum() for a in e])-target,e.min()-20*SIGMA,e.max()+20*SIGMA,xtol=1e-12)
        save(p/'MU.json',dict(mu=mu,electrons=target))
    wait(p/'MU.json')

def density_stage(L,case,z,paths,dv,st):
    p=folder(L,case);out=p/f'density_z{z}.npz'
    if out.exists():return
    start=time.monotonic();mu=json.loads((p/'MU.json').read_text())['mu'];acc=np.zeros((L*L,NG,4));N=np.zeros((L*L,10,10),complex);bare=0.;F=0.;band=0.;count=0.
    for ki,k in points(L,z):
        b=np.load(p/f'k_{ki:04}.npz');e,Q=b['e'],b['Q'];f,_,free=fermi(e,mu);P=(Q*f)@Q.conj().T;h0,aa=operators(L,k,*st)
        psi=np.load(paths[ki],mmap_mode='r');acc+=fold_density(psi,P,L)/32/dv
        N+=np.array([a@P@a.conj().T for a in aa])/32;bare+=trace(P,h0)/32;F+=free/32;band+=float(e@f)/32;count+=float(f.sum())/32
        print('DENSITY_POINT',L,case,ki,flush=True)
    np.savez(out,density=acc,N=N,bare=bare,F=F,band=band,count=count,seconds=time.monotonic()-start)
    print('DENSITY_CASE',L,case,z,time.monotonic()-start,flush=True)

def assemble_density(L,case,dv):
    p=folder(L,case);rows=[np.load(p/f'density_z{i}.npz') for i in range(8)];bg=np.load(R/'reference/background.npz');idx=np.array([(x%2)*2+y%2 for x,y in product(range(L),repeat=2)])
    rho=sum(r['density'] for r in rows)+bg['density'][idx];N=sum(r['N'] for r in rows)+bg['N'][idx]
    charge=float(rho[:,:,0].sum()*dv);assert abs(charge-41*L*L)<1e-6,(L,case,charge)
    np.savez(p/'DENSITY.npz',N=N);write_density(p/'DENSITY.bin',cells_to_grid(rho,L),L)
    save(p/'ELECTRONIC.json',{key:float(sum(r[key] for r in rows)) for key in ['bare','F','band','count']})

def native(L):
    out=R/f'L{L:02}/native';out.mkdir(parents=True,exist_ok=True)
    if (out/'COMPLETE.json').exists():return
    input_file(out,L);(out/'input_list.txt').write_text('11\n'+''.join(f'case_{i:02} "{folder(L,i)}/DENSITY.bin" 0 0 1\n' for i in range(11)))
    with (out/'native.out').open('w') as log:subprocess.run([str(Path(ENERGY['native_binary'])),'grid.in'],cwd=out,stdout=log,stderr=subprocess.STDOUT,check=True)
    assert 'GGA_BRANCH_REEXPORT_COMPLETE' in (out/'native.out').read_text()
    for case in range(11):
        p=folder(L,case);ee=np.loadtxt(out/f'case_{case:02}_ENERGY.dat');n=np.load(p/'DENSITY.npz')['N'];eu,_=hubbard(n);b=json.loads((p/'ELECTRONIC.json').read_text());terms=dict(band_free=b['F'],minus_active_potential=-b['band']+b['bare'],EH=ee[0]*RY,Exc=ee[1]*RY,minus_core=-ee[2]*RY,ionic_density=ee[6]*RY,EU=eu)
        save(p/'ENERGY.json',dict(terms=terms,energy=sum(terms.values()),constant_reference_excluded=True))
    save(out/'COMPLETE.json',dict(done=True))

def potential_stage(L,case,z,paths,st):
    p=folder(L,case);mu=json.loads((p/'MU.json').read_text())['mu'];field=grid_to_cells(read_density(R/f'L{L:02}/native/case_{case:02}_TRUE.bin',L),L)*RY
    nt=np.load(p/'DENSITY.npz')['N'];_,vu=hubbard(nt);numerator=0.;denominator=0.
    for ki,k in points(L,z):
        dest=p/f'adjoint_{ki:04}.npz'
        if dest.exists():
            a=np.load(dest);numerator+=float(a['numerator']);denominator+=float(a['denominator']);continue
        b=np.load(p/f'k_{ki:04}.npz');e,Q,H=b['e'],b['Q'],b['H'];h0,aa=operators(L,k,*st);psi=np.load(paths[ki],mmap_mode='r')
        A=H-h0-fold_project(psi,field,L)-sum(a.conj().T@v@a for a,v in zip(aa,vu));ae=Q.conj().T@A@Q;f,fp,_=fermi(e,mu)
        num=float(np.dot(fp,np.diag(ae).real))/32;den=float(fp.sum())/32;np.savez(dest,Ae=ae,numerator=num,denominator=den);numerator+=num;denominator+=den
        print('NATIVE_PROJECTED',L,case,ki,flush=True)
    save(p/f'alpha_z{z}.json',dict(numerator=numerator,denominator=denominator))

def torque_stage(L,z,rvr):
    definition=json.loads((ROOT/f'inputs/L{L:02}/CASE.json').read_text());targets=[v['site_native_1based'] for v in definition['selected_shells'][0]['pairs']]
    for chunk in [list(range(0,4)),list(range(4,8)),list(range(8,11))]:
        entries=[];densities=[];ds=[]
        for case in chunk:
            p=folder(L,case);mu=json.loads((p/'MU.json').read_text())['mu'];al=[json.loads((p/f'alpha_z{i}.json').read_text()) for i in range(8)];alpha=sum(a['numerator'] for a in al)/sum(a['denominator'] for a in al)
            sites=[1] if case==0 or case>=7 else [1,*targets]
            derivs=[derivative(oldcase(L,case)/f'coeff/dsite_{s:04}.bin') for s in sites]
            for ki,k in points(L,z):
                dest=p/f'torque_k{ki:04}.npz'
                if dest.exists():continue
                b=np.load(p/f'k_{ki:04}.npz');e,Q=b['e'],b['Q'];ae=np.load(p/f'adjoint_{ki:04}.npz')['Ae'];ae.flat[::len(e)+1]-=alpha;f,fp,_=fermi(e,mu);Z=Q@(loewner(e,f,fp)*ae)@Q.conj().T;G=fold_unitary(L,k)
                densities.append(-G.conj().T@Z@G);entries.append((case,ki,k,dest,sites));ds.append(derivs)
        if entries:
            vals,im,stats=contract_batch(rvr,OLD/'basis/pure_H0.au',np.array([a[2] for a in entries]),np.array(densities),ds[0],workers=1,derivatives_by_point=ds)
            for row,v,imag in zip(entries,vals,im):np.savez(row[3],torque=v,sites=row[4],imaginary=imag)
            save(R/f'L{L:02}/torque_batch_{z}_{chunk[0]}.json',stats)

def collect(L):
    import csv
    d=json.loads((ROOT/f'inputs/L{L:02}/CASE.json').read_text());shell=d['selected_shells'][0];rec=[]
    for case in range(11):
        rows=[np.load(folder(L,case)/f'torque_k{i:04}.npz') for i in range(32)];rec.append(dict(sites=rows[0]['sites'].tolist(),torque=np.mean([r['torque'] for r in rows],axis=0).tolist()))
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
    definition=json.loads((ROOT/f'inputs/L{L:02}/CASE.json').read_text())
    targets=[p['site_native_1based'] for p in definition['selected_shells'][0]['pairs']]
    assert sorted(targets)==sorted([2,L,L+1,L*(L-1)+1])
    expected=[(None,0),(1,.25),(1,-.25),(1,.5),(1,-.5),(1,1),(1,-1),(2,.25),(2,-.25),(2,.5),(2,-.5)]
    for i,(site,angle) in enumerate(expected):
        case=definition['angle_configurations'][i]
        assert (case['perturb_site_native_1based'],case['angle_deg'])==(site,angle)
        if i%8!=z:continue
        out=oldcase(L,i);out.mkdir(parents=True,exist_ok=True)
        base=dict(np.load(ROOT/f'inputs/L{L:02}/structure.npz'))
        if site is not None:
            a=np.deg2rad(angle);v=base['spins'][site-1].copy()
            base['spins'][site-1]=[np.cos(a)*v[0]+np.sin(a)*v[2],v[1],-np.sin(a)*v[0]+np.cos(a)*v[2]]
        np.savez(out/'structure.npz',**base)
        sites=[1] if i==0 or i>=7 else [1,*targets]
        coefficients(out/'structure.npz',out/'coeff',sites,.001)
    barrier(f'coefficients_L{L}',z)

def main(z):
    compute_node();require_regression();binding();assert frozen()==13
    if not 0<=z<8 or not SIZES or len(set(SIZES))!=len(SIZES) or any(L not in (4,6,8) for L in SIZES):
        raise ValueError('Require kz worker 0..7 and distinct sizes from 4,6,8')
    R.mkdir(parents=True,exist_ok=True)
    for L in SIZES:
        if str(L) not in ENERGY.get('reuse_cases',{}):prepare_coefficients(L,z)
    meta=json.loads((INTERFACE/'V_func_tmp/V_FUNC.json').read_text());dv=abs(np.linalg.det(np.array(meta['primitive_cell_A'])/.529177210903))/NG
    cache,slab=stage_bloch(z);st=static(OLD/'basis/pure_H0.au');reference(z,slab,dv);background(z,dv)
    cached=CONFIG.get('rvr_slices',{}).get(str(z))
    rvr=Path(cached) if cached else cache/'RVR.bin'
    if not rvr.exists():stage_slice(OLD/'RVR/B_general',rvr,z)
    prefix,sectors,_,size=layout(rvr)
    original,_,_,full_size=layout(OLD/'RVR/B_general')
    if sectors!=(z/8,) or size!=full_size or prefix[:32]!=original[:32] or prefix[36:]!=original[36:]:
        raise ValueError('Cached RVR does not match the full source and kz')
    save(R/f'reused_RVR_z{z}.json',dict(path=str(rvr),bytes=size))
    for L in SIZES:
        paths=wave_cache(L,z,cache,slab)
        for case in range(11):
            spectrum(L,case,z,paths,cache,rvr);density_stage(L,case,z,paths,dv,st);barrier(f'density_L{L}_c{case}',z)
            if z==0:assemble_density(L,case,dv)
        barrier(f'all_density_L{L}',z)
        if z==0:native(L)
        wait(R/f'L{L:02}/native/COMPLETE.json')
        for case in range(11):potential_stage(L,case,z,paths,st)
        barrier(f'potentials_L{L}',z);torque_stage(L,z,rvr);barrier(f'torques_L{L}',z)
        if z==0:collect(L)
    if z==0:save(R/'COMPLETE.json',dict(V1_files_verified=frozen(),sizes=SIZES,J='J12 only'))

if __name__=='__main__':
    try:main(int(sys.argv[1]))
    except BaseException as e:
        if not (R/'FAILURE.json').exists():save(R/'FAILURE.json',dict(worker=int(sys.argv[1]),error=repr(e)))
        raise
