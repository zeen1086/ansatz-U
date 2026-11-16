"""Compact real-space functions; no training grid, DFT labels, or MPI needed."""
from itertools import product
import numpy as np

CUTOFF=6.0

def window(points):
    r=np.linalg.norm(points,axis=-1)
    return np.where(r<CUTOFF,np.cos(np.pi*r/(2*CUTOFF))**2,0.)


def local_evaluate(points,coefficients,geo):
    """Trilinear hats in actual Cartesian local coordinates; zero outside 6 Angstrom."""
    half=geo['half'];shape=2*half+1
    lookup=geo.get('lookup')
    if lookup is None:
        lookup=np.full(tuple(shape),-1,dtype=np.int32)
        lookup[tuple((geo['nodes']+half).T)]=np.arange(len(geo['nodes']))
        geo['lookup']=lookup
    out=np.zeros((len(points),coefficients.shape[1]))
    active=np.flatnonzero(np.linalg.norm(points,axis=1)<CUTOFF)
    pos=points[active]/geo['spacing_A'];low=np.floor(pos).astype(int);f=pos-low
    for dx in (0,1):
        for dy in (0,1):
            for dz in (0,1):
                shift=np.array([dx,dy,dz]);index=low+shift+half
                valid=np.all((index>=0)&(index<shape),axis=1)
                j=np.flatnonzero(valid);rows=lookup[tuple(index[j].T)];ok=rows>=0
                j=j[ok];rows=rows[ok]
                w=np.prod(np.where(shift,f[j],1-f[j]),axis=1)
                out[active[j]]+=coefficients[rows]*w[:,None]
    out*=window(points)[:,None]
    return out


def periodic_evaluate(points,center,cell,coefficients,geo):
    inv=np.linalg.inv(cell);frac=(points-center)@inv
    bounds=CUTOFF*np.linalg.norm(inv,axis=0)
    low=np.ceil(-frac.max(axis=0)-bounds).astype(int)
    high=np.floor(-frac.min(axis=0)+bounds).astype(int)
    result=np.zeros((len(points),coefficients.shape[1]))
    for shift in product(*(range(a,b+1) for a,b in zip(low,high))):
        relative=points-center+np.array(shift)@cell
        result+=local_evaluate(relative,coefficients,geo)
    return result
