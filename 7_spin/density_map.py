"""Common Wannier wave-function readout; preserve native cal_Hsk gauge and WS support."""
from pathlib import Path
from itertools import product
import numpy as np

GRID=np.array([40,40,36]);NG=int(GRID.prod())

def ws_weights(q: np.ndarray) -> np.ndarray:
    """Weights mapping the original 8^3 Bloch samples to the existing weighted WS interpolation."""
    rr=np.arange(-4,5);end=np.where(abs(rr)==4,.5,1.)
    weights=[]
    for x in q:
        samples=np.arange(8)/8
        weights.append(np.sum(end[:,None]*np.exp(2j*np.pi*rr[:,None]*(samples[None,:]-x)),axis=0)/8)
    return np.einsum('i,j,k->ijk',*weights)

def bloch_at(path: Path,q: np.ndarray) -> np.ndarray:
    src=np.load(path,mmap_mode='r');assert src.shape==(8,8,8,NG,30)
    weights=ws_weights(q);out=np.zeros((NG,30),complex)
    for idx in zip(*np.nonzero(abs(weights)>1e-13)):out+=weights[idx]*src[idx]
    return out

def cell_wave(bloch: np.ndarray,k: np.ndarray,L: int,cell: np.ndarray) -> np.ndarray:
    """One primitive-cell grid block, all supercell spatial orbitals (native x-major order)."""
    cells=np.array(list(product(range(L),repeat=2)));q=(cells+k[:2])/L;n=L*L
    mix=np.exp(2j*np.pi*((cell-cells[:,None,:])*q[None,:,:]).sum(axis=2))/n
    # mix[a,q] for site a and primitive Bloch q.
    return np.einsum('aq,qgo->gao',mix,bloch,optimize=True).reshape(bloch.shape[1],30*n)

def density(w: np.ndarray,rho: np.ndarray) -> np.ndarray:
    up=(w@rho[::2,::2])*w.conj();down=(w@rho[1::2,1::2])*w.conj()
    cross=np.sum((w@rho[::2,1::2])*w.conj(),axis=1)
    u=up.sum(axis=1).real;d=down.sum(axis=1).real
    return np.column_stack([u+d,2*cross.real,-2*cross.imag,u-d])

def project(w: np.ndarray,field: np.ndarray) -> np.ndarray:
    n=w.shape[1];out=np.zeros((2*n,2*n),complex)
    out[::2,::2]=w.conj().T@((field[:,0]+field[:,3])[:,None]*w)
    out[1::2,1::2]=w.conj().T@((field[:,0]-field[:,3])[:,None]*w)
    out[::2,1::2]=w.conj().T@((field[:,1]-1j*field[:,2])[:,None]*w)
    out[1::2,::2]=out[::2,1::2].conj().T
    return out

def selfcheck():
    rng=np.random.default_rng(1612);w=rng.normal(size=(19,6))+1j*rng.normal(size=(19,6));v=rng.normal(size=(12,4))+1j*rng.normal(size=(12,4));rho=v@v.conj().T;f=rng.normal(size=(19,4))
    lhs=np.einsum('ij,ji->',rho,project(w,f)).real;rhs=np.sum(density(w,rho)*f)
    assert abs(lhs-rhs)<1e-9,(lhs,rhs)
    for q in [[0,0,0],[.25,.375,.875]]:
        ww=ws_weights(np.array(q));ix=tuple(np.rint(np.array(q)*8).astype(int));assert abs(ww[ix]-1)<1e-14 and np.sum(abs(ww))-abs(ww[ix])<1e-12
    print('density/projector Pauli adjoint and WS interpolation selfchecks PASS')
if __name__=='__main__':selfcheck()

def fold_unitary(L: int,k: np.ndarray) -> np.ndarray:
    cells=np.array(list(product(range(L),repeat=2)));q=(cells+k[:2])/L
    return np.kron(np.exp(-2j*np.pi*(q@cells.T))/L,np.eye(60))

def fold_density(psi: np.ndarray,P: np.ndarray,L: int) -> np.ndarray:
    """Primitive-cell FFT of the density; O(L^4), avoiding an extra factor L^2."""
    n=L*L;g=psi.shape[1];acc=np.zeros((L,L,g,4),complex);cells=np.array(list(product(range(L),repeat=2)))
    for a in range(n):
        for b in range(n):
            block=P[60*a:60*(a+1),60*b:60*(b+1)]
            uu,dd,ud,du=block[::2,::2],block[1::2,1::2],block[::2,1::2],block[1::2,::2]
            mats=[uu+dd,ud+du,1j*(ud-du),uu-dd];d=tuple((cells[a]-cells[b])%L)
            batch=np.concatenate(mats,axis=1)
            tmp=(psi[a]@batch).reshape(g,4,30)
            acc[d]+=np.einsum('gco,go->gc',tmp,psi[b].conj())
    values=np.fft.ifftn(acc,axes=(0,1)).reshape(n,g,4)
    assert np.max(abs(values.imag))<1e-8,float(np.max(abs(values.imag)))
    return values.real

def fold_project(psi: np.ndarray,field: np.ndarray,L: int) -> np.ndarray:
    """Adjoint of fold_density for a local four-component real field."""
    n=L*L;g=psi.shape[1];ft=np.fft.fftn(field.reshape(L,L,g,4),axes=(0,1))/n
    out=np.zeros((60*n,60*n),complex);cells=np.array(list(product(range(L),repeat=2)))
    for a in range(n):
        for b in range(a,n):
            f=ft[tuple((cells[a]-cells[b])%L)]
            v=np.stack([f[:,0]+f[:,3],f[:,1]-1j*f[:,2],f[:,1]+1j*f[:,2],f[:,0]-f[:,3]],axis=1)
            temp=(v[:,:,None]*psi[b][:,None,:]).reshape(g,120)
            val=(psi[a].conj().T@temp).reshape(30,4,30)
            block=np.empty((60,60),complex)
            block[::2,::2]=val[:,0,:];block[::2,1::2]=val[:,1,:];block[1::2,::2]=val[:,2,:];block[1::2,1::2]=val[:,3,:]
            out[60*a:60*(a+1),60*b:60*(b+1)]=block
            out[60*b:60*(b+1),60*a:60*(a+1)]=block.conj().T
    return out

def cells_to_grid(a: np.ndarray,L: int) -> np.ndarray:
    trailing=a.shape[2:];out=np.empty((40*L,40*L,36,*trailing),a.dtype)
    for i,(x,y) in enumerate(product(range(L),repeat=2)):
        out[40*x:40*(x+1),40*y:40*(y+1)]=a[i].reshape((40,40,36,*trailing),order='F')
    return out.reshape((NG*L*L,*trailing),order='F')

def grid_to_cells(a: np.ndarray,L: int) -> np.ndarray:
    trailing=a.shape[1:];v=a.reshape((40*L,40*L,36,*trailing),order='F')
    return np.stack([v[40*x:40*(x+1),40*y:40*(y+1)].reshape((NG,*trailing),order='F') for x,y in product(range(L),repeat=2)])

def fft_selfcheck():
    rng=np.random.default_rng(17);L=2;n=4;psi=rng.normal(size=(n,7,30))+1j*rng.normal(size=(n,7,30));v=rng.normal(size=(60*n,3))+1j*rng.normal(size=(60*n,3));P=v@v.conj().T;f=rng.normal(size=(n,7,4));k=np.array([.2,.3,0]);G=fold_unitary(L,k)
    d=fold_density(psi,P,L);h=fold_project(psi,f,L)
    assert abs(np.einsum('ij,ji->',P,h).real-np.sum(d*f))<1e-8
    for i,c in enumerate(product(range(L),repeat=2)):
        w=cell_wave(psi,k,L,np.array(c))
        assert np.max(abs(d[i]-density(w,G.conj().T@P@G)))<1e-9
    print('primitive-cell FFT density/projector/gauge checks PASS')
