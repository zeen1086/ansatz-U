"""Load the frozen predictor and validate its compact spatial kernels."""
from pathlib import Path
import json
import numpy as np
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
