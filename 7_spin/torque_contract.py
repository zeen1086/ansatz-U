"""Exact adjoint of native local cal_Hsk: one RVR pass for all feedback derivatives."""
from concurrent.futures import ThreadPoolExecutor
from contextlib import nullcontext
from itertools import product
from pathlib import Path
from time import monotonic
import numpy as np


def read_array(f,dtype,shape):
    a=np.fromfile(f,dtype=dtype,count=int(np.prod(shape)))
    if a.size!=np.prod(shape) or not np.isfinite(a).all():raise ValueError('Truncated/nonfinite native array')
    return a.reshape(shape,order='F')


def derivative(path: Path):
    with path.open('rb') as f:
        if f.read(16)!=b'LANIO2_COEFF_V2 ':raise ValueError('Invalid coefficient magic')
        head=read_array(f,'<i4',(8,));model=f.read(64).decode();ueff,offset=read_array(f,'<f8',(2,))
        if head[0]!=1 or head[4]!=1 or offset!=0:raise ValueError('Feedback derivative coefficients required')
        n=int(head[3]);np_=int(head[1]);length=int(head[2])
        if n!=length**2:raise ValueError('Derivative supercell size mismatch')
        p=read_array(f,'<f8',(np_,4,n)).transpose(0,2,1)
        dn=read_array(f,'<c16',(10,10,n)).transpose(2,0,1)
        if f.read(1) or np.max(abs(dn-dn.conj().swapaxes(1,2)))>1e-10:raise ValueError('Invalid complete dN')
    return dict(model_id=model,L=length,p=p,N=dn,Ueff=float(ueff))


def projector_density(background: Path,length: int,k: np.ndarray,rho: np.ndarray) -> np.ndarray:
    with background.open('rb') as f:
        if f.read(16)!=b'LANIO2_AU_V1    ':raise ValueError('Invalid background magic')
        head=read_array(f,'<i4',(8,));nr=int(head[1])
        if head[0]!=16909060 or head[2]!=30 or head[3]!=5:raise ValueError('Unsupported projector layout')
        rr=read_array(f,'<i4',(3,nr)).T;f.seek(16*900*nr,1)
        pr=read_array(f,'<c16',(5,30,nr)).transpose(2,0,1)
    cells=np.array(list(product(range(length),repeat=2)));n=length**2
    q=np.column_stack(((cells+k[:2])/length,np.full(n,k[2])))
    pq=(np.exp(2j*np.pi*(q@rr.T))@pr.reshape(nr,-1)).reshape(length,length,5,30)
    pg=np.fft.ifftn(pq,axes=(0,1));out=[]
    for cell in cells:
        delta=cell-cells;wrapped=delta%length
        values=pg[wrapped[:,0],wrapped[:,1]]*np.exp(2j*np.pi*(delta@k[:2])/length)[:,None,None]
        a=np.zeros((10,60*n),complex);block=values.transpose(1,0,2).reshape(5,30*n)
        a[:5,::2]=block;a[5:,1::2]=block
        out.append((a@rho)@a.conj().T)
    return np.array(out)


def contract_batch(rvr: Path, background: Path, points: np.ndarray, densities: np.ndarray,
                   derivatives: list[dict], workers: int = 1, *,
                   derivatives_by_point: list[list[dict]] | None = None,
                   projectors_by_point: list[np.ndarray] | None = None):
    """Stream each orbital-pair block once for up to sixteen points at one kz."""
    points = np.asarray(points)
    if points.ndim != 2 or points.shape[1] != 3 or not 1 <= len(points) <= 16:
        raise ValueError('Expected 1..16 three-component k points')
    if not np.isfinite(points).all() or np.max(abs((points[:,2]-points[0,2]+.5)%1-.5)) > 1e-10:
        raise ValueError('A torque batch must have one common kz')
    series = derivatives_by_point if derivatives_by_point is not None else [derivatives]*len(points)
    if len(series) != len(points) or not all(series) or not 1 <= workers <= 16:
        raise ValueError('Expected feedback derivatives and 1..16 workers')
    length = series[0][0]['L']
    n = length**2
    cells = np.array(list(product(range(length), repeat=2)))
    if densities.shape != (len(points),60*n,60*n):
        raise ValueError('One complete density matrix is required per k point')
    for rho in densities:
        if not np.isfinite(rho).all() or np.max(abs(rho-rho.conj().T)) > 1e-10:
            raise ValueError('Complete Hermitian density matrix required')

    if projectors_by_point is not None:
        if len(projectors_by_point)!=len(points) or any(a.shape!=(n,10,60*n) or not np.isfinite(a).all() for a in projectors_by_point):
            raise ValueError('Invalid cached projectors for torque batch')

    read_seconds = 0.
    with rvr.open('rb') as f:
        if f.read(16)!=b'LANIO2_LRVR_V2  ':raise ValueError('Invalid RVR magic')
        head=read_array(f,'<i4',(8,));model=f.read(64).decode();f.read(64)
        np_,nr,nw,nk,ndiff=map(int,head[1:6]);nc=2*np_+1
        if head[0]!=1 or nw!=30 or min(np_,nr,nk,ndiff)<=0:raise ValueError('Invalid RVR dimensions')
        for point_derivatives in series:
            for d in point_derivatives:
                if d['model_id']!=model or d['L']!=length or d['p'].shape!=(np_,n,4):
                    raise ValueError('Derivative/RVR identity mismatch')
        read_array(f,'<f8',(3,3));read_array(f,'<f8',(1,))
        means=read_array(f,'<f8',(np_+1,))
        xy=read_array(f,'<i4',(2,nr)).T;difference=read_array(f,'<i4',(2,ndiff)).T
        kz=read_array(f,'<f8',(nk,));index=np.flatnonzero(abs((kz-points[0,2]+.5)%1-.5)<1e-10)
        if len(index)!=1:raise ValueError('Missing/duplicate projected kz')
        block_bytes=16*(ndiff+nr*nr*nc);start=f.tell()
        block_count=nw*nw
        if rvr.stat().st_size!=start+nk*block_count*block_bytes:raise ValueError('RVR byte count mismatch')
        f.seek(start+int(index[0])*block_count*block_bytes)

        raw = cells[:,None,:]+xy
        wrapped = raw%length
        shift = (raw-wrapped)//length
        site = wrapped[:,:,0]*length+wrapped[:,:,1]
        phases = np.exp(2j*np.pi*(shift@points[:,:2].T))
        pair_phases = [phase.conj()[:,:,None]*phase[:,None,:] for phase in phases.transpose(2,0,1)]
        raw = cells[:,None,:]+difference
        wrapped = raw%length
        shift = (raw-wrapped)//length
        target = wrapped[:,:,0]*length+wrapped[:,:,1]
        overlap_phases = np.exp(2j*np.pi*(shift@points[:,:2].T)).transpose(2,0,1)
        gamma = np.zeros((len(points),np_,n,4),complex)
        trace_overlap = np.zeros(len(points),complex)

        def update(slot, ii, jj, a, b, overlap, kernel):
            rho = densities[slot]
            uu=rho[jj,ii];dd=rho[jj+1,ii+1];du=rho[jj+1,ii];ud=rho[jj,ii+1]
            channels=np.stack([uu+dd,du+ud,1j*(ud-du),uu-dd],axis=1)*pair_phases[slot][:,None]
            vectors=channels.transpose(2,3,0,1).reshape(nr*nr,n*4,order='F')
            sensitivity=(kernel.T@vectors).reshape(nc,n,4,order='F')
            gamma[slot,:,:,0] += sensitivity[:np_,:,0]
            gamma[slot,:,:,1:] += sensitivity[np_+1:,:,1:]
            trace_overlap[slot] += np.sum((rho[b,a]+rho[b+1,a+1])*overlap_phases[slot]*overlap)

        # The next block is read only after all contractions using the current kernel finish.
        workers = min(workers, len(points))
        pool_context = ThreadPoolExecutor(max_workers=workers) if workers>1 else nullcontext()
        with pool_context as pool:
            for i in range(nw):
                ii=(site*60+2*i)[:,:,None]
                a=(np.arange(n)*60+2*i)[:,None]
                for j in range(nw):
                    tick=monotonic()
                    overlap=read_array(f,'<c16',(ndiff,));kernel=read_array(f,'<c16',(nr*nr,nc))
                    read_seconds+=monotonic()-tick
                    jj=(site*60+2*j)[:,None,:]
                    b=target*60+2*j
                    args=(ii,jj,a,b,overlap,kernel)
                    if pool is None:
                        for slot in range(len(points)):update(slot,*args)
                    else:
                        futures=[pool.submit(update,slot,*args) for slot in range(len(points))]
                        for future in futures:future.result()
        gamma[:,:,:,0] -= trace_overlap[:,None,None]*means[None,:np_,None]/n

    gradients = []
    for slot,(k, rho, sensitivity, point_derivatives) in enumerate(zip(points,densities,gamma,series)):
        if projectors_by_point is None:
            projected = projector_density(background,length,k,rho)
        else:
            aa=projectors_by_point[slot]
            projected=np.array([(a@rho)@a.conj().T for a in aa])
        gradients.append(np.array([np.sum(sensitivity*d['p'])-d['Ueff']*np.einsum('sij,sji->',projected,d['N'])
                                   for d in point_derivatives]))
    if any(not np.isfinite(row).all() for row in gradients):raise ValueError('Nonfinite HF gradient')
    imaginary=np.array([np.max(abs(row.imag)) for row in gradients])
    if np.max(imaginary)>1e-8:raise ValueError(f'Complex HF gradient: {np.max(imaginary)}')
    values = [-row.real for row in gradients]
    return (values if derivatives_by_point is not None else np.asarray(values)), imaginary, dict(
        rvr_blocks=block_count,rvr_bytes=block_count*block_bytes,rvr_read_and_validation_seconds=read_seconds)


def contract(rvr: Path, background: Path, k: np.ndarray, rho: np.ndarray, derivatives: list[dict]):
    """Compatibility interface for one k point."""
    values, imaginary, _ = contract_batch(rvr,background,k[None],rho[None],derivatives)
    return values[0], float(imaginary[0])
