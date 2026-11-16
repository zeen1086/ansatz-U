"""Fixed paths and small numerical helpers for the reviewed V1 material run."""
from pathlib import Path
import json
import os
import subprocess
import numpy as np
from scipy.special import erfc

SOURCE = Path(__file__).resolve().parents[2]
CONFIG_PATH = Path(os.environ.get('ANSATZ_RUN_CONFIG', SOURCE/'examples/7_spin/production.json')).resolve()
CONFIG = json.loads(CONFIG_PATH.read_text())
ROOT = Path(CONFIG['input_root'])
INTERFACE = Path(CONFIG['interface'])
BUILD = Path(CONFIG['native_build'])
REF = Path(CONFIG['references'])
WORK = Path(CONFIG['work'])
PYTHON = CONFIG['python']
RY = 13.605693122994
SIGMA = .01*RY


def save(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix+'.partial')
    temp.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')
    temp.replace(path)


def compute_node() -> None:
    if 'ANSATZ_RUN_CONFIG' not in os.environ:
        raise RuntimeError('Set ANSATZ_RUN_CONFIG to an external production configuration')
    if 'PBS_JOBID' not in os.environ:
        raise RuntimeError('Material computations require the reviewed PBS compute job')


def run(*args: object) -> None:
    subprocess.run([str(x) for x in args], check=True)


def load_h(path: Path, dimension: int) -> np.ndarray:
    if path.stat().st_size != 16*dimension**2:
        raise ValueError(f'Truncated Hamiltonian: {path}')
    h = np.fromfile(path, '<c16').reshape(dimension, dimension, order='F')
    if not np.isfinite(h).all() or np.max(abs(h-h.conj().T)) > 1e-8:
        raise ValueError(f'Invalid Hamiltonian: {path}')
    return h


def chemical_potential(e: np.ndarray, weights: np.ndarray, electrons: int) -> float:
    if e.ndim != 2 or weights.shape != (len(e),) or abs(weights.sum()-1) > 1e-12:
        raise ValueError('Expected all-k spectra and normalized k weights')
    lo, hi = float(e.min()-20*SIGMA), float(e.max()+20*SIGMA)
    for _ in range(100):
        mu = (lo+hi)/2
        count = weights@np.sum(.5*erfc((e-mu)/SIGMA), axis=1)
        if count < electrons:
            lo = mu
        else:
            hi = mu
    if abs(count-electrons) > 1e-8:
        raise ValueError('Fixed electron count failed')
    return mu


def require_regression() -> None:
    report = json.loads((WORK/'regression/RESULT.json').read_text())
    if not report['passed'] or report['model_id'] != CONFIG['model_id']:
        raise RuntimeError('The real 2x2 Hamiltonian/band regression must pass first')


def require_calculation(length: int) -> None:
    require_regression()
