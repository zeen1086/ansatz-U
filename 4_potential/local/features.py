"""Fixed smooth local features from structure and unit Ni spin directions."""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product
import numpy as np

from symmetry import SPATIAL, finite, require


@dataclass
class Structure:
    cell: np.ndarray                 # lattice vectors are rows, Angstrom
    positions: np.ndarray            # Cartesian, Angstrom; any cell images
    species: np.ndarray              # atomic numbers
    ni: np.ndarray                   # index into positions, ordered as N labels
    spins: np.ndarray                # (nNi,3), unit directions only

    def validate(self) -> None:
        require(self.cell.shape == (3, 3) and abs(np.linalg.det(self.cell)) > 1e-10,
                "invalid cell")
        require(self.positions.ndim == 2 and self.positions.shape[1] == 3,
                "positions must be (nat,3)")
        require(self.species.shape == (len(self.positions),), "species shape")
        require(self.ni.ndim == 1 and len(self.ni) > 0 and
                np.issubdtype(self.ni.dtype, np.integer), "Ni index dtype/shape")
        require(len(set(self.ni)) == len(self.ni) and
                np.all((self.ni >= 0) & (self.ni < len(self.positions))), "Ni indices")
        require(set(np.flatnonzero(self.species == 28)) == set(self.ni),
                "ni must enumerate every Ni atom exactly once")
        require(self.spins.shape == (len(self.ni), 3), "spin shape")
        for name in ("cell", "positions", "spins"):
            finite(getattr(self, name), name)
        require(np.max(np.abs(np.linalg.norm(self.spins, axis=1) - 1)) < 1e-10,
                "spins must already be unit directions, not DFT moments")


@dataclass(frozen=True)
class FeatureSpec:
    cutoff: float
    lengths: tuple[float, ...]
    species: tuple[int, ...] = (8, 28, 57)
    scalar_degree: int = 1
    vector_degree: int = 1

    def __post_init__(self) -> None:
        require(np.isfinite(self.cutoff) and self.cutoff > 0, "positive cutoff required")
        require(len(self.lengths) > 0 and np.all(np.isfinite(self.lengths)) and
                min(self.lengths) > 0, "positive fixed radial lengths required")
        require(len(set(self.species)) == len(self.species) and 28 in self.species,
                "unique species including Ni required")


def neighbors(s: Structure, center: int, cutoff: float) -> tuple[np.ndarray, ...]:
    """Enumerate all real displacements in support, including repeated images."""
    inv = np.linalg.inv(s.cell)
    delta = (s.positions - s.positions[center]) @ inv
    bounds = cutoff * np.linalg.norm(inv, axis=0)
    ds, atom = [], []
    for j, fractional in enumerate(delta):
        lo = np.ceil(-fractional - bounds).astype(int)
        hi = np.floor(-fractional + bounds).astype(int)
        for image in product(*(range(a, b + 1) for a, b in zip(lo, hi))):
            d = (fractional + image) @ s.cell
            r = np.linalg.norm(d)
            if 1e-12 < r < cutoff:
                ds.append(d); atom.append(j)
    return np.array(ds).reshape(-1, 3), np.array(atom, dtype=int)


def smooth_weights(d: np.ndarray, spec: FeatureSpec) -> np.ndarray:
    r2 = np.sum(d * d, axis=1)
    t = r2 / spec.cutoff ** 2
    inside = t < 1
    bump = np.zeros_like(t)
    bump[inside] = np.exp(-t[inside] / (1 - t[inside]))
    xyz = d / spec.cutoff
    x, y, z = xyz.T
    poly = np.array([np.ones(len(d)), x, y, z, x*x, y*y, z*z, x*y, x*z, y*z]).T
    radial = np.exp(-r2[:, None] / np.square(spec.lengths)) * bump[:, None]
    return (radial[:, :, None] * poly[:, None, :]).reshape(len(d), 10*len(spec.lengths))


def features(s: Structure, spec: FeatureSpec) -> dict[str, np.ndarray]:
    s.validate()
    require(set(s.species).issubset(spec.species), "unregistered structural species")
    all_spins = np.zeros((len(s.positions), 3)); all_spins[s.ni] = s.spins
    out = {k: [] for k in ("se", "so", "vo", "ve")}
    for site, center in enumerate(s.ni):
        d, atom = neighbors(s, center, spec.cutoff)
        e = s.spins[site]
        row = {k: [] for k in out}
        for action in SPATIAL:
            # g^{-1} E makes z's group axis a regular representation.
            weights = smooth_weights(d @ action, spec)
            geom = np.concatenate([weights[s.species[atom] == z].sum(axis=0)
                                   for z in spec.species])
            moment = weights.T @ all_spins[atom]
            dot = moment @ e
            next_m = np.roll(moment, 1, axis=0)
            cross = np.cross(moment, next_m)
            se = np.concatenate(([1.], geom, dot, dot*dot,
                                 np.sum(moment*moment, axis=1)))[:, None]
            row["se"].append(se)
            row["so"].append((cross @ e)[:, None])
            row["vo"].append(np.concatenate((e[None], moment, dot[:, None]*e,
                                               dot[:, None]*moment), axis=0))
            row["ve"].append(np.concatenate((np.cross(e, moment), cross), axis=0))
        for k in out:
            out[k].append(row[k])
    out={k: np.array(v) for k,v in out.items()}
    if spec.scalar_degree==4 or spec.vector_degree>=3:
        q=10*len(spec.lengths);start=1+len(spec.species)*q
        pick=np.array([0,4,5]);ni=spec.species.index(28)
        denominator=out['se'][:,:,1+ni*q+pick,0]
        a=out['se'][:,:,start+pick,0]/np.maximum(denominator,1e-14)
        def monomials(degrees):
            exponents=np.array([e for degree in degrees for e in product(range(degree+1),repeat=3) if sum(e)==degree])
            return np.prod(a[:,:,None,:]**exponents[None,None,:,:],axis=-1)
        if spec.scalar_degree==4:
            out['se']=np.concatenate([out['se'],monomials((3,4))[...,None]],axis=2)
        if spec.vector_degree>=3:
            vector=monomials(tuple(range(2,spec.vector_degree+1)))[...,None]*out['vo'][:,:,:1,:]
            out['vo']=np.concatenate([out['vo'],vector],axis=2)
    return out
