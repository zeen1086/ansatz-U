"""Eight PBS workers, one lossless local RVR slice and one writer per kz."""
from pathlib import Path
import hashlib
import json
import os
import shutil
import socket
import struct
import time

import numpy as np
from scipy.special import erfc

from common import CONFIG, ROOT, WORK, SIGMA, chemical_potential, compute_node, load_h, require_regression, save
from run_material import binding, case_root, coefficients, native_batch
from torque_contract import contract_batch, derivative

FULL_RVR = WORK/'RVR/B_general'
BASE = WORK/'basis/pure_H0.au'


def read(path):
    return json.loads(path.read_text())


def layout(source):
    with source.open('rb') as f:
        first=f.read(176)
        version,np_,nr,nw,nk,ndiff,_,_=struct.unpack_from('<8i',first,16)
        if first[:16]!=b'LANIO2_LRVR_V2  ' or version!=1 or nw!=30 or min(np_,nr,nk,ndiff)<=0:
            raise ValueError('Invalid full RVR header')
        prefix=bytearray(first+f.read(80+8*(np_+1)+8*(nr+ndiff)))
        sectors=struct.unpack(f'<{nk}d',f.read(8*nk))
        offset=f.tell()
    size=16*(ndiff+nr*nr*(2*np_+1))*nw*nw
    if source.stat().st_size!=offset+nk*size:raise ValueError('RVR payload byte count mismatch')
    return prefix,sectors,offset,size


def stage_slice(source, target, z):
    tick=time.monotonic()
    prefix,sectors,offset,size=layout(source)
    if not 0<=z<len(sectors) or abs(sectors[z]-z/8)>1e-12:
        raise ValueError(f'Sector {z} is absent from source grid')
    if shutil.disk_usage(target.parent).free<size+10*1024**3:
        raise RuntimeError('Insufficient node-local space for RVR slice')
    stamp=source.stat()
    struct.pack_into('<i',prefix,16+4*4,1)
    digest=hashlib.sha256()
    with source.open('rb') as src,target.open('xb') as dst:
        dst.write(prefix)
        dst.write(struct.pack('<d',sectors[z]))
        src.seek(offset+z*size)
        left=size
        while left:
            block=src.read(min(left,16*1024**2))
            if not block:raise EOFError('Truncated source RVR slice')
            dst.write(block)
            digest.update(block)
            left-=len(block)
        dst.flush()
        os.fsync(dst.fileno())
    if (stamp.st_size,stamp.st_mtime_ns)!=(source.stat().st_size,source.stat().st_mtime_ns):
        raise RuntimeError('RVR source changed during staging')
    local_prefix,k,local_offset,local_size=layout(target)
    if local_prefix!=prefix or k!=(sectors[z],) or local_size!=size:
        raise ValueError('Local RVR slice metadata mismatch')
    copied=hashlib.sha256()
    with target.open('rb') as f:
        f.seek(local_offset)
        for block in iter(lambda:f.read(16*1024**2),b''):
            copied.update(block)
    if copied.digest()!=digest.digest():raise ValueError('Local RVR payload digest mismatch')
    return dict(z=z,seconds=time.monotonic()-tick,bytes=size,sha256=digest.hexdigest(),host=socket.gethostname())


def data(task):
    L,m,h,i=task
    definition,case,p=case_root(L,i,m,h)
    k=np.array([[x/m,y/m,z/8,1/(m*m*8)] for x in range(m) for y in range(m) for z in range(8)])
    return definition,case,p,k


def validate_prediction(p,L,h,sites):
    v=read(p/'coeff/PREDICTION.json')
    if (v['status'],v['model_id'],v['L'],v['derivative_step_rad'],v['derivative_sites'])!=(
            'COMPLETE',CONFIG['model_id'],L,h,sites):
        raise ValueError(f'Frozen coefficient provenance mismatch: {p}')


def prepare(task):
    L,m,h,i=task
    _,case,p,_=data(task)
    if not (p/'coeff/PREDICTION.json').is_file():
        p.mkdir(parents=True,exist_ok=True)
        s=dict(np.load(ROOT/f'inputs/L{L:02}/structure.npz'))
        site=case['perturb_site_native_1based']
        if site is not None:
            x,y,z=s['spins'][site-1].copy()
            t=np.deg2rad(case['angle_deg'])
            s['spins'][site-1]=[np.cos(t)*x+np.sin(t)*z,y,-np.sin(t)*x+np.cos(t)*z]
        np.savez(p/'structure.npz',**s)
        coefficients(p/'structure.npz',p/'coeff',case['readout_sites_native_1based'],h)
    validate_prediction(p,L,h,case['readout_sites_native_1based'])


def spectrum_sector(task,z,rvr):
    definition,_,p,k=data(task)
    if (p/'SPECTRA.json').is_file():return
    indices=np.flatnonzero(k[:,2]==z/8)
    nd=definition['dimension']
    for start in range(0,len(indices),16):
        batch=indices[start:start+16]
        dirs=[p/f'k_{i:04}'/'value' for i in batch]
        native_batch(p/'coeff/value.bin',k[batch],dirs,p,nd,rvr=rvr)
        # A single BLAS call can use all 16 cores for an eigensystem.
        for i in batch:
            folder=p/f'k_{i:04}'
            if (folder/'eigensystem.npz').exists():continue
            e,v=np.linalg.eigh(load_h(folder/'value/H.k000001.h',nd))
            temp=folder/'eigensystem.partial.npz'
            np.savez(temp,e=e,v=v);temp.replace(folder/'eigensystem.npz')


def spectrum_finish(task):
    L,m,h,_=task
    definition,case,p,k=data(task)
    if (p/'SPECTRA.json').is_file():return
    energies=[]
    for i in range(len(k)):
        with np.load(p/f'k_{i:04}/eigensystem.npz') as f:energies.append(f['e'])
    e=np.asarray(energies)
    mu=chemical_potential(e,k[:,3],definition['electrons'])
    temp=p/'SPECTRA.partial.npz'
    np.savez(temp,k=k,e=e,mu=mu);temp.replace(p/'SPECTRA.npz')
    save(p/'SPECTRA.json',dict(status='COMPLETE',L=L,case=case,mesh=[m,m,8],derivative_step_rad=h,
         electrons=definition['electrons'],mu_eV=mu,model_id=CONFIG['model_id']))


def torque_inputs(tasks,z,pilot_root=None):
    entries=[]
    for task in tasks:
        L,m,h,_=task
        definition,case,p,k=data(task)
        sites=case['readout_sites_native_1based']
        validate_prediction(p,L,h,sites)
        sp=read(p/'SPECTRA.json')
        if sp['mesh']!=[m,m,8] or sp['derivative_step_rad']!=h or sp['model_id']!=CONFIG['model_id']:
            raise ValueError(f'Spectrum provenance mismatch: {p}')
        with np.load(p/'SPECTRA.npz') as s:
            if not np.array_equal(s['k'],k):raise ValueError(f'K mesh changed: {p}')
            mu=float(s['mu']);energies=s['e']
        ds=[derivative(p/f'coeff/dsite_{site:04}.bin') for site in sites]
        for i in np.flatnonzero(k[:,2]==z/8):
            target=(pilot_root or p)/f'k_{i:04}'/'TORQUE.npz'
            if target.exists():continue
            entries.append((k[i,:3],p/f'k_{i:04}/eigensystem.npz',target,ds,mu,energies[i],definition['dimension']))
    return entries


def torque_sector(tasks,z,rvr,pilot_root=None):
    entries=torque_inputs(tasks,z,pilot_root)
    for start in range(0,len(entries),16):
        block=entries[start:start+16]
        tick=time.monotonic()
        rho=np.empty((len(block),block[0][-1],block[0][-1]),complex)
        stats=[]
        for slot,(_,eig,_,_,mu,expected,nd) in enumerate(block):
            with np.load(eig) as f:e,v=f['e'],f['v']
            if not np.array_equal(e,expected) or v.shape!=(nd,nd):
                raise ValueError(f'Eigensystem mismatch: {eig}')
            x=(e-mu)/SIGMA
            occ=.5*erfc(x)
            rho[slot]=(v*occ)@v.conj().T
            stats.append(dict(count=occ.sum(), F=np.sum(e*occ-SIGMA*np.exp(-x*x)/(2*np.sqrt(np.pi)))))
            del v
        points=np.array([e[0] for e in block])
        ds=[e[3] for e in block]
        ans,imag,measure=contract_batch(rvr,BASE,points,rho,[],workers=len(block),derivatives_by_point=ds)
        del rho
        for slot,(_,_,p,*_) in enumerate(block):
            p.parent.mkdir(parents=True,exist_ok=True)
            tmp=p.with_name('TORQUE.partial.npz')
            np.savez(tmp,torque=ans[slot],imaginary=imag[slot],**stats[slot]);tmp.replace(p)
        print('TORQUE_KZ_LOCAL',json.dumps(dict(z=z,L=tasks[0][0],points=len(block),seconds=time.monotonic()-tick,
                                                source=str(rvr),**measure)),flush=True)


def torque_finish(task,out=None):
    L,m,h,i=task
    definition,case,p,k=data(task)
    out=out or p
    if (out/'TORQUE.json').exists():return read(out/'TORQUE.json')
    rows=[]
    for j in range(len(k)):
        with np.load(out/f'k_{j:04}/TORQUE.npz') as f:
            rows.append({key:f[key] for key in ('torque','count','F','imaginary')})
    weights=k[:,3]
    charge_error=float(abs(weights@np.array([v['count'] for v in rows])-definition['electrons']))
    values=weights@np.array([v['torque'] for v in rows])
    if charge_error>1e-8 or not np.isfinite(values).all():raise ValueError('Invalid integrated HF result')
    result=dict(status='COMPLETE',L=L,case=case,sites=case['readout_sites_native_1based'],mesh=[m,m,8],
        derivative_step_rad=h,torque_eV=values.tolist(),electron_error=charge_error,
        F_eV=float(weights@np.array([v['F'] for v in rows])),
        imaginary_gradient_eV=max(float(v['imaginary']) for v in rows),workflow='local_kz')
    save(out/'TORQUE.json',result)
    return result



if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase',choices=['stage','prepare','spectrum','spectrum-finish','torque','torque-finish'])
    parser.add_argument('--z',type=int,choices=range(8),required=True)
    parser.add_argument('--rvr',type=Path,required=True)
    parser.add_argument('--tasks',type=Path,help='JSON list of [L, mesh, inner_step_rad, case_index]')
    args=parser.parse_args();compute_node();require_regression();binding()
    if args.phase=='stage':
        print(json.dumps(stage_slice(FULL_RVR,args.rvr,args.z),indent=2))
    else:
        if args.tasks is None:parser.error('--tasks is required for calculation phases')
        tasks=read(args.tasks)
        if args.phase=='torque':torque_sector(tasks,args.z,args.rvr)
        else:
            for task in tasks:
                if args.phase=='prepare':prepare(task)
                elif args.phase=='spectrum':spectrum_sector(task,args.z,args.rvr)
                elif args.phase=='spectrum-finish':spectrum_finish(task)
                else:torque_finish(task)
