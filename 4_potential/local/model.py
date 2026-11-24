"""One explicit local occupation response, followed by scalar and magnetic fields."""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from features import FeatureSpec, Structure, features
from symmetry import (REP, SECTORS, SPATIAL, finite, matrix_sigmoid, require,
                      sigmoid_pullback, to_channels, to_matrix)


@dataclass
class Reference:
    n_scalar: np.ndarray
    n_spin: np.ndarray
    p_scalar: np.ndarray             # physical eV expansion coefficients
    p_spin: np.ndarray
    source: str                     # registered provenance, or explicit synthetic

    @classmethod
    def from_anchor(cls, N: np.ndarray, p_eV: np.ndarray, axis: np.ndarray,
                    p_representation: np.ndarray, source: str) -> tuple[Reference, dict]:
        """Declared axial/D4h averaging of a mapped single-Ni anchor; keep its residual."""
        require(N.shape == (10, 10) and p_eV.ndim == 2 and p_eV.shape[1] == 4,
                "anchor matrix/coefficient shapes")
        finite(N, "anchor N"); finite(p_eV, "anchor p"); finite(axis, "anchor axis")
        require(axis.shape == (3,) and abs(np.linalg.norm(axis)-1) < 1e-10, "unit anchor axis")
        require(np.max(np.abs(N-N.conj().T)) < 1e-10, "Hermitian anchor required")
        require(p_representation.shape == (16, len(p_eV), len(p_eV)), "anchor representation shape")
        c = to_channels(N)
        nr = REP["se"].mean(axis=0); pr = p_representation.mean(axis=0)
        ref = cls(nr @ c["se"][:, 0], nr @ (c["vo"] @ axis),
                  pr @ p_eV[:, 0], pr @ (p_eV[:, 1:] @ axis), source)
        nfit, pfit = ref.values(axis[None], np.ones_like(p_eV))
        return ref, dict(N_discarded_Frobenius=float(np.linalg.norm(N-nfit[0])),
                         p_discarded_eV_norm=float(np.linalg.norm(p_eV-pfit[0])))

    def values(self, spins: np.ndarray, scale: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        n = {k: np.zeros((len(spins), a, c)) for k, (a, c) in SECTORS.items()}
        n["se"][:, :, 0] = self.n_scalar
        n["vo"] = self.n_spin[None, :, None] * spins[:, None, :]
        p = np.empty((len(spins), len(self.p_scalar), 4))
        p[:, :, 0] = self.p_scalar
        p[:, :, 1:] = self.p_spin[None, :, None] * spins[:, None, :]
        return to_matrix(n), p / scale


@dataclass
class Prediction:
    p: np.ndarray                   # normalized; use physical_p() at kernel boundary
    N: np.ndarray
    VU: np.ndarray                   # eV
    auxiliary: np.ndarray | None
    converged: bool
    iterations: int
    residual_p: float
    residual_N: float
    contraction_bound: float
    cache: dict

    def assert_converged(self) -> None:
        if not self.converged:
            raise RuntimeError(f"fixed point not converged: iterations={self.iterations}, "
                               f"r_p={self.residual_p:.3e}, r_N={self.residual_N:.3e}")


class LocalModel:
    def __init__(self, spec: FeatureSpec, p_scales: np.ndarray,
                 reference: Reference, rho: float | None = None, *, p_representation: np.ndarray,
                 ueff: float = 1., seed: int = 0):
        self.spec, self.reference = spec, reference
        self.ueff = float(ueff)
        require(rho is None or 0 < rho < 1, "invalid legacy contraction threshold")
        require(np.isfinite(ueff) and ueff >= 0, "invalid Ueff")
        p_scales = np.array(p_scales, dtype=float, copy=True)
        require(p_scales.ndim == 2 and p_scales.shape[1] == 2 and len(p_scales) > 0,
                "p_scales must be (ncoeff,2), scalar/B eV scales")
        require(np.isfinite(p_scales).all() and np.all(p_scales > 0), "positive fixed scales")
        self.p_scales = p_scales
        self.scale = p_scales[:, [0, 1, 1, 1]].copy()
        self.scale.flags.writeable = False; self.p_scales.flags.writeable = False
        self.np = len(self.scale)
        self.prep = np.array(p_representation, dtype=float, copy=True)
        require(self.prep.shape == (16, self.np, self.np), "spatial coefficient representation shape")
        finite(self.prep, "spatial coefficient representation")
        require(np.max(np.abs(self.prep @ self.prep.swapaxes(-1, -2) - np.eye(self.np))) < 1e-10,
                "spatial representation must be orthogonal")
        for a, action in enumerate(SPATIAL):
            for b, other in enumerate(SPATIAL):
                c = next(i for i, target in enumerate(SPATIAL) if np.array_equal(action @ other, target))
                require(np.max(np.abs(self.prep[a] @ self.prep[b] - self.prep[c])) < 1e-10,
                        "spatial representation must follow the registered D4h group")
        for scale in p_scales.T:
            require(np.max(np.abs(self.prep*scale[None, None, :] - self.prep*scale[None, :, None])) < 1e-10,
                    "fixed scales must commute with spatial symmetry")
        self.prep.flags.writeable = False
        for name, size in (("n_scalar", 15), ("n_spin", 15),
                           ("p_scalar", self.np), ("p_spin", self.np)):
            a = np.array(getattr(reference, name), dtype=float, copy=True)
            require(a.shape == (size,), f"reference {name} shape")
            finite(a, name)
            rep = self.prep if name.startswith("p_") else REP["se"]
            require(np.max(np.abs(np.einsum("gab,b->ga", rep, a) - a)) < 1e-10,
                    f"reference {name} must respect local D4h")
            a.flags.writeable = False
            setattr(reference, name, a)
        require(bool(reference.source), "reference provenance required")
        q = 10 * len(spec.lengths)
        dims = {"se": 1 + (len(spec.species) + 3)*q, "so": q,
                "vo": 1 + 3*q, "ve": 2*q}
        if spec.scalar_degree==4:dims['se']+=25
        if spec.vector_degree>=3:dims['vo']+=sum((d+1)*(d+2)//2 for d in range(2,spec.vector_degree+1))
        self.feature_scales = {k: np.ones(d) for k, d in dims.items()}
        self.feature_transforms = {k:np.eye(d) for k,d in dims.items()}
        self.p_feature_maps=None
        rng = np.random.default_rng(seed)
        self.params = {}
        for k, (orb, _) in SECTORS.items():
            self.params[f"Q_{k}"] = rng.normal(scale=.02, size=(orb, dims[k]))
        for k in ("se", "vo"):
            self.params[f"p_{k}"] = rng.normal(scale=.02, size=(self.np, dims[k]))
            self.params[f"b_p_{k}"] = np.zeros(self.np)
        for k in ("se", "vo"):
            self.params[f"b_Q_{k}"] = np.zeros(15)
            self.params[f"K_{k}"] = np.zeros((self.np, 15))
            self.params[f"A_{k}"] = np.zeros((15, 15))

    def zero_grads(self, active=None) -> dict[str, np.ndarray]:
        return {k: np.zeros_like(v) for k, v in self.params.items() if active is None or k in active}

    def prepare_features(self,raw):
        return {k:np.einsum('ngrc,rj->ngjc',v/self.feature_scales[k][None,None,:,None],self.feature_transforms[k]) for k,v in raw.items()}

    def field_features(self,z):
        if self.p_feature_maps is None:return z
        out={};character=np.r_[np.ones(8),-np.ones(8)]
        for sector in ('se','vo'):
            value=z[sector];maps=self.p_feature_maps[sector]
            plus=np.einsum('irc,rj->ijc',value.mean(axis=1),maps['plus'])
            minus=np.einsum('irc,rj->ijc',np.einsum('g,igrc->irc',character,value)/16,maps['minus'])
            out[sector]=np.zeros_like(value);n=plus.shape[1];m=minus.shape[1]
            out[sector][:,:,:n]=plus[:,None]
            out[sector][:,:,n:n+m]=character[None,:,None,None]*minus[:,None]
        return out

    def occupation_to_field_maps(self) -> dict[str, np.ndarray]:
        return {f"K_{k}": np.einsum("gab,bc,gdc->ad", self.prep, self.params[f"K_{k}"], REP[k])
                / len(self.prep) for k in ("se", "vo")}

    def occupation_response_maps(self) -> dict[str, np.ndarray]:
        return {f"A_{k}": np.einsum("gab,bc,gdc->ad", REP[k], self.params[f"A_{k}"], REP[k])
                / len(REP[k]) for k in ("se", "vo")}

    def _head(self, name: str, z: dict, spins: np.ndarray) -> dict:
        if name=='p':z=self.field_features(z)
        out = {}
        for k in (SECTORS if name == "Q" else ("se", "vo")):
            raw = np.einsum("ar,ngrc->ngac", self.params[f"{name}_{k}"], z[k])
            if k in ("se", "vo"):
                direction = np.ones((len(spins), 1)) if k == "se" else spins
                raw += self.params[f"b_{name}_{k}"][None, None, :, None] * direction[:, None, None, :]
            rep = REP[k] if name == "Q" else self.prep
            out[k] = np.einsum("gab,ngbc->nac", rep, raw) / len(rep)
        return out

    def _head_pullback(self, name: str, grad: dict, cache: dict, out: dict) -> None:
        phi=self.field_features(cache['phi']) if name=='p' else cache['phi']
        for k, g in grad.items():
            weight, bias = f"{name}_{k}", f"b_{name}_{k}"
            if weight not in out and bias not in out: continue
            rep = REP[k] if name == "Q" else self.prep
            gr = np.einsum("gab,nac->ngbc", rep, g) / len(rep)
            if weight in out:
                out[weight] += np.einsum("ngac,ngrc->ar", gr, phi[k])
            if bias in out:
                direction = np.ones((len(g), 1)) if k == "se" else cache["spins"]
                out[bias] += np.einsum("ngac,nc->a", gr, direction)

    @staticmethod
    def pack_p(x: dict) -> np.ndarray:
        return np.concatenate((x["se"], x["vo"]), axis=-1)

    @staticmethod
    def split_p(p: np.ndarray) -> dict:
        return {"se": p[:, :, :1], "vo": p[:, :, 1:]}

    @staticmethod
    def _K(maps: dict, n: np.ndarray) -> np.ndarray:
        nc = to_channels(n)
        return LocalModel.pack_p({k: np.einsum("ab,nbc->nac", maps[f"K_{k}"], nc[k])
                                  for k in ("se", "vo")})

    @staticmethod
    def _A(maps: dict, n: np.ndarray) -> np.ndarray:
        nc = to_channels(n)
        response = {k: np.zeros_like(v) for k, v in nc.items()}
        for k in ("se", "vo"):
            response[k] = np.einsum("ab,nbc->nac", maps[f"A_{k}"], nc[k])
        return to_matrix(response)

    def predict(self, structure: Structure, *, auxiliary: bool = False,
                tolerance: float = 1e-10, max_iterations: int = 100,
                prepared_features: dict[str, np.ndarray] | None = None,
                normalized: bool = False) -> Prediction:
        require(np.isfinite(tolerance) and tolerance > 0 and
                isinstance(max_iterations, int) and max_iterations >= 0, "invalid solver limits")
        require(not auxiliary, "direct model has no auxiliary density head")
        raw = features(structure, self.spec) if prepared_features is None else prepared_features
        phi = {}
        for k, (_, c) in SECTORS.items():
            require(raw[k].shape == (len(structure.ni), 16, self.params[f'Q_{k}'].shape[1], c),
                    'cached feature shape mismatch')
            scale = self.feature_scales[k]
            require(np.all(np.isfinite(scale)) and np.all(scale > 0), 'invalid fixed feature scale')
            phi[k] = raw[k] if normalized else np.einsum('ngrc,rj->ngjc',raw[k]/scale[None,None,:,None],self.feature_transforms[k])
        p0 = self.pack_p(self._head("p", phi, structure.spins))
        q0 = to_matrix(self._head("Q", phi, structure.spins))
        n0, initial_cache = matrix_sigmoid(q0)
        nr, _ = self.reference.values(structure.spins, self.scale)
        maps = self.occupation_to_field_maps()
        maps.update(self.occupation_response_maps())
        initial_delta = n0 - nr
        n, sigmoid_cache = matrix_sigmoid(q0 + self._A(maps, initial_delta))
        delta = n - nr
        p = p0 + self._K(maps, delta)
        finite(p, "predicted field coefficients")
        cache = dict(phi=phi, spins=structure.spins.copy(), maps=maps,
                     initial_sigmoid=initial_cache, initial_delta=initial_delta,
                     sigmoid=sigmoid_cache, delta=delta)
        return Prediction(p, n, self.ueff*(np.eye(10)/2 - n), None, True, 0, 0., 0., 0., cache)

    def physical_p(self, prediction: Prediction) -> np.ndarray:
        return prediction.p * self.scale

    def backward(self, pred: Prediction, gp: np.ndarray, gn: np.ndarray,
                 gaux: np.ndarray | None = None, *, active=None) -> dict[str, np.ndarray]:
        """Differentiate both sigmoid evaluations and the explicit A response."""
        require(gaux is None, "direct model has no auxiliary density loss")
        cache, out = pred.cache, self.zero_grads(active)
        field_grad = self.split_p(gp)
        occupation_grad = to_channels(gn)
        delta = to_channels(cache["delta"])
        for k, g in field_grad.items():
            key = f"K_{k}"
            occupation_grad[k] += np.einsum("ab,nac->nbc", cache["maps"][key], g)
            if key in out:
                g_map = np.einsum("nac,nbc->ab", g, delta[k])
                out[key] += np.einsum("gab,ad,gdc->bc", self.prep, g_map, REP[k]) / len(self.prep)
        gq = sigmoid_pullback(cache["sigmoid"], to_matrix(occupation_grad))
        response_grad = to_channels(gq)
        initial_delta = to_channels(cache["initial_delta"])
        initial_grad = {k: np.zeros_like(v) for k, v in response_grad.items()}
        for k in ("se", "vo"):
            key = f"A_{k}"
            initial_grad[k] = np.einsum("ab,nac->nbc", cache["maps"][key], response_grad[k])
            if key in out:
                g_map = np.einsum("nac,nbc->ab", response_grad[k], initial_delta[k])
                out[key] += np.einsum("gab,ad,gdc->bc", REP[k], g_map, REP[k]) / len(REP[k])
        gq += sigmoid_pullback(cache["initial_sigmoid"], to_matrix(initial_grad))
        self._head_pullback("p", field_grad, cache, out)
        self._head_pullback("Q", to_channels(gq), cache, out)
        return out
