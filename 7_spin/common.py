"""Fixed paths and small numerical helpers for the reviewed V1 material run."""
from pathlib import Path
import json
import os
import subprocess
import numpy as np
from scipy.special import erfc

SOURCE = Path(__file__).resolve().parents[1]
CONFIG_PATH = Path(os.environ.get('ANSATZ_RUN_CONFIG', SOURCE/'examples/7_spin/production.json')).resolve()
CONFIG = json.loads(CONFIG_PATH.read_text())
ROOT = Path(CONFIG['input_root'])
INTERFACE = Path(CONFIG['interface'])
BUILD = Path(CONFIG['native_build'])
REF = Path(CONFIG['references'])
WORK = Path(CONFIG['work'])  # Shared RVR and basis inputs.
SPINS = Path(CONFIG.get('spin_work', CONFIG['work']))
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
    report = json.loads((SPINS/'regression/RESULT.json').read_text())
    if not report['passed'] or report['model_id'] != CONFIG['model_id']:
        raise RuntimeError('The real 2x2 Hamiltonian/band regression must pass first')


def require_calculation(length: int) -> None:
    require_regression()


def file_digest(path: Path) -> str:
    """Hash contents, rejecting a source modified while it is being read."""
    import hashlib
    before = path.stat()
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(8*1024*1024), b''):
            h.update(chunk)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns, before.st_ctime_ns) != (after.st_size, after.st_mtime_ns, after.st_ctime_ns):
        raise ValueError(f'Input changed during hashing: {path}')
    return h.hexdigest()


def shared_bundle(path: Path, identity: dict, build) -> Path:
    """Publish an immutable bundle atomically; validate it once per worker startup."""
    import fcntl
    import tempfile
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.with_name(path.name+'.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if path.exists():
            try:
                record = json.loads((path/'CACHE.json').read_text())
                if record['identity'] != identity:
                    raise ValueError('input identity mismatch')
                files = record['files']
                if not files or set(files) != {p.name for p in path.iterdir() if p.name != 'CACHE.json'}:
                    raise ValueError('incomplete cache file set')
                for name, info in files.items():
                    p = path/name
                    if p.stat().st_size != info['bytes'] or file_digest(p) != info['sha256']:
                        raise ValueError(f'corrupt payload: {name}')
            except (OSError, KeyError, TypeError, ValueError) as exc:
                raise ValueError(f'Invalid shared cache {path}: {exc}; use a new shared directory or archive this entry') from exc
            return path
        with tempfile.TemporaryDirectory(prefix=path.name+'.partial-', dir=path.parent) as tmp:
            staging = Path(tmp)/'bundle'
            staging.mkdir()
            build(staging)
            files = {}
            for p in staging.iterdir():
                if not p.is_file() or p.name == 'CACHE.json':
                    raise ValueError('Cache builder must produce only payload files')
                if p.suffix == '.npy':
                    a = np.load(p, mmap_mode='r', allow_pickle=False)
                    for row in a:
                        if not np.isfinite(row).all():
                            raise ValueError(f'Nonfinite cache payload: {p.name}')
                files[p.name] = dict(bytes=p.stat().st_size, sha256=file_digest(p))
            if not files:
                raise ValueError('Empty shared cache')
            save(staging/'CACHE.json', dict(identity=identity, files=files))
            staging.rename(path)
        return path
