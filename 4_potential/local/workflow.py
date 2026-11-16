"""Native spatial reconstruction and versioned direct-model storage."""
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import numpy as np

from features import FeatureSpec
from model import LocalModel, Reference
from symmetry import finite, require


@dataclass(frozen=True)
class LossConfig:
    v: float
    B: float
    N: float
    n: float
    m: float
    lambda_density: float
    multiply_N_by_U2: bool = False

    def __post_init__(self) -> None:
        require(all(np.isfinite(x) and x > 0 for x in (self.v, self.B, self.N, self.n, self.m)),
                'loss scales must be frozen positive values')
        require(self.lambda_density == 0, 'direct model excludes auxiliary density loss')


def save_model(model: LocalModel, directory: Path, training_record: dict | None = None) -> None:
    directory.mkdir(parents=True, exist_ok=False)
    meta = dict(schema='LANIO2_ACTIVE_LOCAL_GRID_FEEDBACK_V2', feature=asdict(model.spec),
                p_scales=model.p_scales.tolist(), p_representation=model.prep.tolist(),
                rho=model.rho, ueff=model.ueff, v_channels=getattr(model,'v_channels',model.np),
                head_coordinates='Qeff=Q+L(p0-pref); exact compensated feedback',
                feedback_coordinates={k:v.tolist() for k,v in model.feedback_coordinates.items()},
                feature_scales={k: v.tolist() for k, v in model.feature_scales.items()},
                feature_transforms={k:v.tolist() for k,v in model.feature_transforms.items()},
                p_feature_maps=None if model.p_feature_maps is None else {k:{p:a.tolist() for p,a in v.items()} for k,v in model.p_feature_maps.items()},
                reference={k: v.tolist() if isinstance(v, np.ndarray) else v
                           for k, v in asdict(model.reference).items()},
                orbital_order=['dz2', 'dxz', 'dyz', 'dx2-y2', 'dxy'],
                spin_order='up d1..5, down d1..5', scalar_gauge='cell_mean_zero',
                training_record=training_record, physical_release=False)
    (directory/'model.json').write_text(json.dumps(meta, indent=2)+'\n')
    np.savez(directory/'parameters.npz', **model.params)


def load_model(directory: Path, warm=False):
    meta = json.loads((directory/'model.json').read_text())
    schema = 'LANIO2_ACTIVE_LOCAL_GRID_FEEDBACK_V2'
    require(meta['schema'] == schema or (warm and meta['schema'] in ('LANIO2_GEOMETRIC_FEEDBACK_V1','LANIO2_SEPARATE_LOCAL_GRID_FEEDBACK_V1')), 'unsupported model schema')
    ref = Reference(**{k: v if k == 'source' else np.array(v)
                       for k, v in meta['reference'].items()})
    args = (FeatureSpec(**meta['feature']), np.array(meta['p_scales']), ref, meta['rho'])
    model = LocalModel(*args, p_representation=np.array(meta['p_representation']), ueff=meta['ueff'])
    model.feedback_coordinates = {k:np.array(v) for k,v in meta.get('feedback_coordinates', {k:np.eye(model.np).tolist() for k in ('se','vo')}).items()}
    model.v_channels=meta.get('v_channels',model.np)
    model.feature_scales = {k: np.array(v) for k, v in meta['feature_scales'].items()}
    model.feature_transforms={k:np.array(v) for k,v in meta.get('feature_transforms',{k:np.eye(len(v)).tolist() for k,v in model.feature_scales.items()}).items()}
    if meta.get('p_feature_maps') is not None:model.p_feature_maps={k:{p:np.array(a) for p,a in v.items()} for k,v in meta['p_feature_maps'].items()}
    with np.load(directory/'parameters.npz', allow_pickle=False) as data:
        require(set(data.files) == set(model.params), 'model parameter keys mismatch')
        for key, value in model.params.items():
            a = data[key]
            require(a.shape == value.shape and a.dtype == np.float64, f'parameter shape/type: {key}')
            finite(a, key); model.params[key] = a.copy()
    if warm and meta['schema']!='LANIO2_ACTIVE_LOCAL_GRID_FEEDBACK_V2':
        maps,_=model.feedback_maps()
        for sector,column,refname in (('se',0,'p_scalar'),('vo',1,'p_spin')):
            l=maps['L_'+sector]
            model.params['Q_'+sector]+=l@model.params['p_'+sector]
            model.params['b_Q_'+sector]+=l@(model.params['b_p_'+sector]-getattr(model.reference,refname)/model.scale[:,column])
    return model

