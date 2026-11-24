"""Real QE d orbitals; orthonormal full Hermitian spin-major coordinates."""
from __future__ import annotations

import numpy as np
from scipy.special import expit

SECTORS = {"se": (15, 1), "so": (10, 1), "vo": (15, 3), "ve": (10, 3)}


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError(message)


def finite(a: np.ndarray, name: str) -> None:
    require(np.isfinite(a).all(), f"nonfinite {name}")


def orbital_bases() -> tuple[np.ndarray, np.ndarray]:
    sym = [np.diag(np.eye(5)[i]) for i in range(5)]
    anti = []
    for i in range(5):
        for j in range(i + 1, 5):
            a = np.zeros((5, 5)); a[i, j] = 1 / np.sqrt(2)
            sym.append(a + a.T)
            anti.append(1j * (a - a.T))
    return np.array(sym), np.array(anti)


SYM, ANTI = orbital_bases()
PAULI = np.array([np.eye(2), [[0, 1], [1, 0]],
                  [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]])
BASIS = {
    k: np.array([[np.kron(PAULI[c], orb) / np.sqrt(2)
                  for c in ([0] if n == 1 else [1, 2, 3])]
                 for orb in (SYM if k in ("se", "vo") else ANTI)])
    for k, (_, n) in SECTORS.items()
}


def to_matrix(x: dict[str, np.ndarray]) -> np.ndarray:
    return sum(np.einsum("...ac,acij->...ij", x[k], BASIS[k]) for k in SECTORS)


def to_channels(a: np.ndarray) -> dict[str, np.ndarray]:
    return {k: np.einsum("acij,...ij->...ac", b.conj(), a).real
            for k, b in BASIS.items()}


def d4h() -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """16 local D4h point actions, without four-site shifts."""
    spatial, orbital = [], []
    for sw in (0, 1):
        for sx in (1, -1):
            for sy in (1, -1):
                for sz in (1, -1):
                    a = np.zeros((3, 3)); d = np.zeros((5, 5))
                    a[0, sw] = sx; a[1, 1 - sw] = sy; a[2, 2] = sz
                    d[0, 0] = 1; d[1, 1 + sw] = sx * sz
                    d[2, 2 - sw] = sy * sz; d[3, 3] = (-1) ** sw
                    d[4, 4] = sx * sy
                    spatial.append(a); orbital.append(d)
    rep = {}
    for k in SECTORS:
        b = SYM if k in ("se", "vo") else ANTI
        rep[k] = np.array([np.einsum("aij,bij->ab", b.conj(),
                                     np.einsum("ij,bjk,lk->bil", d, b, d)).real
                            for d in orbital])
    return np.array(spatial), rep


SPATIAL, REP = d4h()


def matrix_sigmoid(q: np.ndarray) -> tuple[np.ndarray, tuple[np.ndarray, ...]]:
    finite(q, "Q")
    require(np.max(np.abs(q - q.swapaxes(-1, -2).conj())) < 1e-10,
            "Q must be Hermitian")
    w, u = np.linalg.eigh(q)
    f = expit(w)
    # Finite precision can saturate: report this instead of clipping Q or N.
    require(np.all((f > 0) & (f < 1)), "matrix sigmoid saturated in float64")
    n = (u * f[..., None, :]) @ u.swapaxes(-1, -2).conj()
    actual = np.linalg.eigvalsh(n)
    require(np.all((actual > 0) & (actual < 1)), "matrix sigmoid lost strict bounds in float64")
    return n, (w, u, f)


def sigmoid_pullback(cache: tuple[np.ndarray, ...], grad: np.ndarray) -> np.ndarray:
    """Self-adjoint Frechet derivative, including degenerate eigenvalues."""
    w, u, f = cache
    dw = w[..., :, None] - w[..., None, :]
    df = f[..., :, None] - f[..., None, :]
    mid = expit((w[..., :, None] + w[..., None, :]) / 2)
    divided = mid * (1 - mid)
    np.divide(df, dw, out=divided, where=np.abs(dw) > 1e-7)
    local = u.swapaxes(-1, -2).conj() @ grad @ u
    out = u @ (divided * local) @ u.swapaxes(-1, -2).conj()
    return (out + out.swapaxes(-1, -2).conj()) / 2
