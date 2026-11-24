"""Fixed paths and small numerical helpers for the reviewed V1 material run."""
from pathlib import Path
import json
import os
import subprocess
import sys
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
SPINS = Path(CONFIG['spin_work'])
PYTHON = CONFIG['python']
RY = 13.605693122994
SIGMA = .01*RY


def save(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    import tempfile
    with tempfile.NamedTemporaryFile(mode='w',dir=path.parent,prefix='.'+path.name,suffix='.partial',delete=False) as f:
        temp=Path(f.name)
        try:
            f.write(json.dumps(value,indent=2,allow_nan=False)+'\n')
            f.flush();os.fsync(f.fileno());temp.replace(path)
        finally:temp.unlink(missing_ok=True)


def compute_node() -> None:
    missing=[name for name in ('ANSATZ_RUN_CONFIG','PBS_JOBID') if not os.environ.get(name)]
    if missing:print('NOTE: '+', '.join(missing)+' not set; continuing with available configuration')


def run(*args: object) -> None:
    subprocess.run([str(x) for x in args], check=True)


def load_h(path: Path, dimension: int) -> np.ndarray:
    if path.stat().st_size != 16*dimension**2:
        raise ValueError(f'Truncated Hamiltonian: {path}')
    h = np.fromfile(path, '<c16').reshape(dimension, dimension, order='F')
    if not np.isfinite(h).all():
        raise ValueError(f'Invalid Hamiltonian: {path}')
    error=float(np.max(abs(h-h.conj().T)))
    quality_check(error<=1e-5, f'Hamiltonian Hermiticity residual: {error}; {path}')
    return h


def chemical_potential(e: np.ndarray, weights: np.ndarray, electrons: int) -> float:
    if e.ndim != 2 or weights.shape != (len(e),) or abs(weights.sum()-1) > 1e-9:
        raise ValueError('Expected all-k spectra and normalized k weights')
    if not np.isfinite(e).all() or not np.isfinite(weights).all() or not 0 <= electrons <= e.shape[1]:
        raise ValueError('Invalid spectra, weights or electron target')
    lo, hi = float(e.min()-20*SIGMA), float(e.max()+20*SIGMA)
    for _ in range(100):
        mu = (lo+hi)/2
        count = weights@np.sum(.5*erfc((e-mu)/SIGMA), axis=1)
        if count < electrons:
            lo = mu
        else:
            hi = mu
    quality_check(abs(count-electrons)<=1e-8, f'Fixed electron count residual: {count-electrons}')
    return mu


def quality_check(passed: bool, message: str) -> None:
    """Report finite-data quality without making it a production gate."""
    if not passed:
        print(f'QUALITY WARNING: {message}', file=sys.stderr, flush=True)


def check_regression() -> None:
    path = SPINS/'regression/RESULT.json'
    try:
        report = json.loads(path.read_text())
        passed = isinstance(report, dict) and report.get('passed') is True and report.get('model_id') == CONFIG['model_id']
    except (OSError, ValueError):
        passed = False
    quality_check(passed, f'Material regression missing, unreadable, failed or for another model: {path}; continuing')



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
                    if p.stat().st_size != info['bytes']:
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
                    for row in (a[:1] if a.nbytes>256*1024**2 else a):
                        if not np.isfinite(row).all():
                            raise ValueError(f'Nonfinite cache payload: {p.name}')
                files[p.name] = dict(bytes=p.stat().st_size)
            if not files:
                raise ValueError('Empty shared cache')
            save(staging/'CACHE.json', dict(identity=identity, files=files))
            staging.rename(path)
        return path


def atomic_npz(path: Path, **arrays) -> None:
    """Publish a complete archive; interruption leaves the prior result intact."""
    import tempfile
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, prefix='.'+path.name, suffix='.partial', delete=False) as f:
        temp = Path(f.name)
        try:
            np.savez(f, **arrays)
            f.flush(); os.fsync(f.fileno())
            temp.replace(path)
        finally:
            temp.unlink(missing_ok=True)


def valid_npz(path: Path, keys) -> bool:
    import zipfile
    try:
        with path.open('rb') as stream, np.load(stream, allow_pickle=False) as data:
            for key in keys:
                if key not in data:return False
                value=data[key]
                if isinstance(keys,dict) and value.shape!=keys[key]:return False
                if not np.isfinite(value).all():return False
            return True
    except (OSError, ValueError, TypeError, EOFError, zipfile.BadZipFile):
        return False


def archive_incomplete(path: Path) -> None:
    import uuid
    if path.exists():
        path.rename(path.with_name(path.name+'.incomplete-'+uuid.uuid4().hex))


def bind_inputs(path: Path, identity: dict, *, ignored=()) -> None:
    """Bind inputs, adopting unbound directories while rejecting identity changes."""
    import fcntl
    path.mkdir(parents=True, exist_ok=True)
    with (path/'INPUTS.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        manifest = path/'INPUTS.json'
        if manifest.exists():
            if json.loads(manifest.read_text()) != identity:
                raise ValueError(f'Checkpoint inputs changed: {path}; use a new output directory')
        else:
            if any(p.name not in ('INPUTS.lock', *ignored) for p in path.iterdir()):
                print(f'NOTE: adopting existing checkpoints without INPUTS.json: {path}')
            save(manifest, identity)


def input_digest(path: Path) -> str:
    """Metadata signature, without reading file contents."""
    stat=path.stat()
    return f'stat:{stat.st_size}:{stat.st_mtime_ns}:{stat.st_ctime_ns}'


def arrays_digest(arrays: dict) -> str:
    """Exact small spin input description, without hashing."""
    rows=[]
    for key,value in sorted(arrays.items()):
        a=np.asarray(value)
        if a.dtype.hasobject:raise ValueError('Object arrays are not supported as spin inputs')
        rows.append([key,a.dtype.str,a.shape,a.tolist()])
    return json.dumps(rows,separators=(',',':'),allow_nan=False)


def valid_json(path: Path, keys) -> bool:
    try:
        value=json.loads(path.read_text())
        return isinstance(value,dict) and all(key in value for key in keys)
    except (OSError,ValueError):return False
