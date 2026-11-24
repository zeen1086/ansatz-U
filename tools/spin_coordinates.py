"""Fixed exponential coordinates for a two-site transverse energy Hessian."""
import numpy as np


CHART_KEYS = ('chart_base_spins', 'chart_sites', 'chart_basis', 'chart_eta')


def tangent_basis(spins: np.ndarray) -> np.ndarray:
    spins = np.asarray(spins, dtype=float)
    if (spins.ndim != 2 or spins.shape[1] != 3 or not np.isfinite(spins).all()
            or not np.allclose(np.linalg.norm(spins, axis=1), 1., rtol=0, atol=1e-10)):
        raise ValueError('Finite unit background spins required')
    axes = np.eye(3)[np.argmin(abs(spins), axis=1)]
    first = axes - np.sum(axes*spins, axis=1)[:, None]*spins
    first /= np.linalg.norm(first, axis=1)[:, None]
    return np.stack((first, np.cross(spins, first)), axis=-1)


def chart_spins(base: np.ndarray, sites: np.ndarray, basis: np.ndarray,
                eta: np.ndarray) -> np.ndarray:
    base = np.asarray(base, dtype=float)
    sites = np.asarray(sites)
    basis = np.asarray(basis, dtype=float)
    eta = np.asarray(eta, dtype=float)
    if (base.ndim != 2 or base.shape[1] != 3 or not np.isfinite(base).all()
            or not np.allclose(np.linalg.norm(base, axis=1), 1., rtol=0, atol=1e-10)):
        raise ValueError('Finite unit background spins required')
    if (sites.shape != (2,) or sites.dtype.kind not in 'iu' or len(set(sites)) != 2
            or np.any(sites < 1) or np.any(sites > len(base))):
        raise ValueError('Exactly two distinct 1-based native sites required')
    if (basis.shape != (2, 3, 2) or eta.shape != (2, 2)
            or not np.isfinite(basis).all() or not np.isfinite(eta).all()
            or not np.allclose(basis.transpose(0, 2, 1) @ basis, np.eye(2), rtol=0, atol=1e-10)
            or np.max(abs(np.einsum('si,sia->sa', base[sites-1], basis))) > 1e-10):
        raise ValueError('Finite orthonormal tangent bases and 2x2 coordinates required')
    r = np.linalg.norm(eta, axis=1)
    if np.any(r >= np.pi):
        raise ValueError('Exponential coordinates must remain inside the antipode')
    result = base.copy()
    result[sites-1] = (np.cos(r)[:, None]*base[sites-1]
                      + np.sinc(r/np.pi)[:, None]*np.einsum('sia,sa->si', basis, eta))
    return result


def derivative_files(sites: list[int], tangent: bool = False) -> list[str]:
    return ([f'dsite_{site:04}_a{alpha}.bin' for site in sites for alpha in (1, 2)]
            if tangent else [f'dsite_{site:04}.bin' for site in sites])


def curvature_cases(sites: list[int], steps: list[float]) -> list[dict]:
    sites_array = np.asarray(sites)
    h = np.asarray(steps, dtype=float)
    if (sites_array.shape != (2,) or sites_array.dtype.kind not in 'iu'
            or len(set(sites)) != 2 or min(sites) < 1):
        raise ValueError('Exactly two distinct 1-based native sites required')
    if (h.shape != (3,) or not np.isfinite(h).all() or h[0] <= 0 or h[-1] >= np.pi
            or not np.allclose(h, h[0]*np.array([1., 2., 4.]), rtol=1e-12, atol=0)):
        raise ValueError('Three positive outer steps h, 2h, 4h in radians required')
    cases = [dict(tangent_sites=list(sites), eta=np.zeros((2, 2)).tolist())]
    for step in h:
        for column in range(4):
            for sign in (1, -1):
                eta = np.zeros(4); eta[column] = sign*step
                cases.append(dict(tangent_sites=list(sites), eta=eta.reshape(2, 2).tolist()))
    return cases
