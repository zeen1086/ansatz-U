"""Independent four-cell reference projection, preserving the original folded phases."""
from itertools import product
import numpy as np
BITS=np.array(list(product([0,1],repeat=2)))

def indices(grid):
    nx,ny,nz=map(int,grid);assert nx%2==ny%2==0
    base=np.arange(nx*ny*nz).reshape((nx,ny,nz),order='F')[:nx//2,:ny//2,:].ravel(order='F')
    return np.array([base+int(tx)*nx//2+int(ty)*ny//2*nx for tx,ty in BITS])


class Projection:
    def general(self,fields,chunk=8):
        n=len(fields);d=self.block;out=np.zeros((n,4*d,4*d),complex)
        ft=np.array([sum((-1)**int(b@t)*fields[:,ix] for t,ix in zip(BITS,self.ix)) for b in BITS])
        for j,l in product(range(4),repeat=2):
            bit=(self.bits[j]-self.bits[l])%2;b=next(i for i,v in enumerate(BITS) if np.array_equal(v,bit))
            left=self.w[:,j*d:(j+1)*d].conj().T;right=self.w[:,l*d:(l+1)*d]
            for start in range(0,n,chunk):
                c=ft[b,start:start+chunk];nr=len(c)
                batch=(c[:,:,None]*right[None]).transpose(1,0,2).reshape(len(right),nr*d)
                value=(left@batch).reshape(d,nr,d).transpose(1,0,2)
                out[start:start+nr,j*d:(j+1)*d,l*d:(l+1)*d]=value
        return out


def spin_matrix(c):
    # Local potential is v0*I+B.sigma, without the occupation Pauli factor 1/2.
    h=np.zeros((240,240),complex)
    h[0::2,0::2]=c[0]+c[3];h[1::2,1::2]=c[0]-c[3]
    h[0::2,1::2]=c[1]-1j*c[2];h[1::2,0::2]=c[1]+1j*c[2]
    return h
