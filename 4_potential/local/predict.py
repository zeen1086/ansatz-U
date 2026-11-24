"""Structure/unit-spin prediction from actual compact local kernels, without Hamiltonians."""
from pathlib import Path
import json
import numpy as np
from local_field import periodic_evaluate
from workflow import load_model
from symmetry import require


class Bundle:
    def __init__(self,path):
        path=Path(path)
        self.model=load_model(path/'model')
        self.cv=np.load(path/'C_v_local.npy',mmap_mode='r',allow_pickle=False)
        self.cb=np.load(path/'C_B_local.npy',mmap_mode='r',allow_pickle=False)
        with np.load(path/'local_grid.npz',allow_pickle=False) as z:self.geo={k:z[k].copy() for k in z.files}
        self.alpha=json.loads((path/'KERNELS.json').read_text())['scalar_bias_coefficient_per_Ni']
        require(self.cv.shape==(len(self.geo['nodes']),self.model.np+1),'scalar kernel shape')
        require(self.cb.shape==(len(self.geo['nodes']),self.model.np),'spin kernel shape')

    def raw_values(self,structure,points,prediction=None):
        """Scalar values before cell-mean gauge; subtract its mean over the full requested grid."""
        structure.validate();points=np.asarray(points,dtype=float)
        require(points.ndim==2 and points.shape[1]==3 and np.isfinite(points).all(),'Cartesian grid points')
        pred=prediction or self.model.predict(structure)
        p=self.model.physical_p(pred);out=np.zeros((len(points),4))
        for start in range(0,len(points),2048):
            part=points[start:start+2048];values=out[start:start+2048]
            for i,center in enumerate(structure.positions[structure.ni]):
                v=periodic_evaluate(part,center,structure.cell,self.cv,self.geo)
                b=periodic_evaluate(part,center,structure.cell,self.cb,self.geo)
                values[:,0]+=v[:,:-1]@p[i,:,0]+v[:,-1]*self.alpha
                values[:,1:]+=b@p[i,:,1:]
        return out,pred

    def grid_values(self,structure,points):
        values,pred=self.raw_values(structure,points)
        values[:,0]-=values[:,0].mean()
        return dict(v_B_eV=values,N=pred.N,VU_eV=pred.VU,p=self.model.physical_p(pred))
