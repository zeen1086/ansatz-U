"""Native spatial reconstruction and versioned direct-model storage."""
from dataclasses import asdict
import json
from pathlib import Path
import numpy as np

from features import FeatureSpec
from model import LocalModel, Reference
from symmetry import REP, finite, require

SINGLE_RESPONSE_SCHEMA = 'LANIO2_ACTIVE_LOCAL_GRID_SINGLE_RESPONSE_V1'
COUPLED_SCHEMA = 'LANIO2_ACTIVE_LOCAL_GRID_FEEDBACK_V2'


def save_model(model: LocalModel, directory: Path, training_record: dict | None = None) -> None:
    directory.mkdir(parents=True, exist_ok=False)
    meta = dict(schema=SINGLE_RESPONSE_SCHEMA, feature=asdict(model.spec),
                p_scales=model.p_scales.tolist(), p_representation=model.prep.tolist(),
                ueff=model.ueff, v_channels=getattr(model,'v_channels',model.np),
                head_coordinates='N0=sigmoid(Q0); N=sigmoid(Q0+A(N0-Nref)); p=p0+K(N-Nref)',
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


def _import_coupled_weights(model: LocalModel, data, meta: dict) -> None:
    # Compile the normalized source maps into A=L K without changing the stored Q head.
    rho = meta['rho']
    require(np.isfinite(rho) and 0 < rho < 1, 'invalid source contraction threshold')
    projected_k = model.occupation_to_field_maps()
    projected_l = {}
    k_norm_sq = l_norm_sq = 0.
    for sector, multiplicity in (('se', 1), ('vo', 3)):
        key = 'L_' + sector
        raw_l = data[key]
        require(raw_l.shape == (15, model.np) and raw_l.dtype == np.float64, f'parameter shape/type: {key}')
        finite(raw_l, key)
        lmap = np.einsum('gab,bc,gdc->ad', REP[sector], raw_l, model.prep) / len(model.prep)
        projected_l[sector] = lmap
        coords = np.array(meta.get('feedback_coordinates', {}).get(sector, np.eye(model.np)))
        k_norm_sq += multiplicity * np.sum((coords @ projected_k['K_' + sector])**2)
        l_norm_sq += multiplicity * np.sum((lmap @ np.linalg.inv(coords))**2)
    norm_product = np.sqrt(k_norm_sq) * np.sqrt(l_norm_sq)
    require(np.isfinite(norm_product), 'source feedback norm overflow')
    scale = min(1., np.sqrt(4*rho/norm_product)) if norm_product else 1.
    for sector in ('se', 'vo'):
        model.params['K_' + sector] *= scale
        model.params['A_' + sector] = scale**2 * projected_l[sector] @ projected_k['K_' + sector]


def load_model(directory: Path):
    meta = json.loads((directory/'model.json').read_text())
    schema = meta['schema']
    require(schema in (SINGLE_RESPONSE_SCHEMA, COUPLED_SCHEMA), 'unsupported model schema')
    ref = Reference(**{k: v if k == 'source' else np.array(v)
                       for k, v in meta['reference'].items()})
    args = (FeatureSpec(**meta['feature']), np.array(meta['p_scales']), ref)
    model = LocalModel(*args, p_representation=np.array(meta['p_representation']), ueff=meta['ueff'])
    model.v_channels=meta.get('v_channels',model.np)
    model.feature_scales = {k: np.array(v) for k, v in meta['feature_scales'].items()}
    model.feature_transforms={k:np.array(v) for k,v in meta.get('feature_transforms',{k:np.eye(len(v)).tolist() for k,v in model.feature_scales.items()}).items()}
    if meta.get('p_feature_maps') is not None:model.p_feature_maps={k:{p:np.array(a) for p,a in v.items()} for k,v in meta['p_feature_maps'].items()}
    with np.load(directory/'parameters.npz', allow_pickle=False) as data:
        expected = set(model.params)
        if schema != SINGLE_RESPONSE_SCHEMA:
            expected.difference_update(('A_se', 'A_vo'))
            expected.update(('L_se', 'L_vo'))
        require(set(data.files) == expected, 'model parameter keys mismatch')
        for key, value in model.params.items():
            if key.startswith('A_') and schema != SINGLE_RESPONSE_SCHEMA:
                continue
            a = data[key]
            require(a.shape == value.shape and a.dtype == np.float64, f'parameter shape/type: {key}')
            finite(a, key); model.params[key] = a.copy()
        if schema != SINGLE_RESPONSE_SCHEMA:
            _import_coupled_weights(model, data, meta)
    return model

