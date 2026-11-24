#!/usr/bin/env python3
"""Rebuild one converged QE save into the current V1 v/B/N label package."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tools.fortran_io import read_array
from qe_labels import BOHR_ANGSTROM, NUMBER, read_save


def sha256(path: Path) -> str:
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def read_scf_log(path: Path, saved: dict) -> dict:
    """Bind a complete original SCF log to the save and restore its GGA branch."""
    text = path.read_text()
    if len(re.findall(r'Program\s+PWSCF\b', text)) != 1 or 'JOB DONE.' not in text:
        raise ValueError('Require one complete original PWSCF log, including JOB DONE.')
    if ('convergence has been achieved' not in text
            or 'convergence NOT achieved' in text or 'Error in routine' in text):
        raise ValueError('SCF log does not report successful convergence')
    atoms = re.findall(r'number of atoms/cell\s*=\s*(\d+)', text)
    if not atoms or any(int(n) != saved['nat'] for n in atoms):
        raise ValueError('SCF log atom count differs from the save')
    energies = re.findall(r'!\s+total energy\s*=\s*(' + NUMBER + r')\s+Ry', text)
    if not energies:
        raise ValueError('Missing final SCF energy in the log')
    energy = float(energies[-1].replace('D', 'e').replace('d', 'e'))
    if not np.isfinite(energy) or abs(energy - saved['total_energy_Ry']) > 1e-6:
        raise ValueError('SCF log final energy differs from the save')
    axes = re.findall(r'Fixed quantization axis for GGA:\s*(' + NUMBER + r')\s+('
                      + NUMBER + r')\s+(' + NUMBER + r')', text)
    axis = np.zeros(3)
    if axes:
        axis = saved['starting_axis']
        printed = np.array([[float(v.replace('D', 'e').replace('d', 'e')) for v in a] for a in axes])
        if not np.isfinite(printed).all() or np.max(abs(printed - axis)) > 5.1e-7:
            raise ValueError('Original GGA axis differs from the XML initial magnetic direction')
    return dict(lsign=bool(axes), ux=axis.tolist(), final_energy_Ry=energy,
                source='original SCF log; full-precision axis from input XML')


def read_fields(directory: Path, saved: dict, branch: dict) -> tuple[np.ndarray, np.ndarray]:
    """Read only the matched constraint-free dense-grid adapter output."""
    rows = (directory/'LABELS.meta').read_text().splitlines()
    if len(rows) != 4 or rows[0].strip() != 'V1_LABEL_FIELDS_RY':
        raise ValueError('Unsupported field export metadata')
    try:
        grid = np.array([int(v) for v in rows[1].split()])
        spin = np.array([float(v) for v in rows[2].split()])
        cell = np.array([float(v) for v in rows[3].split()]).reshape(3, 3) * BOHR_ANGSTROM
    except ValueError as exc:
        raise ValueError('Invalid field export metadata values') from exc
    if grid.shape != (3,) or not np.array_equal(grid, saved['grid_dense']):
        raise ValueError('Export grid differs from the SCF dense grid')
    if (spin.shape != (4,) or not np.isfinite(spin).all()
            or spin[0] != int(branch['lsign'])
            or not np.allclose(spin[1:], branch['ux'], rtol=0, atol=1e-12)):
        raise ValueError('Export GGA branch differs from the original SCF')
    if not np.isfinite(cell).all() or not np.allclose(cell, saved['cell'], rtol=0, atol=1e-8):
        raise ValueError('Export cell differs from the SCF save')
    with (directory/'V_LABELS.bin').open('rb') as stream:
        binary_grid = read_array(stream, '<i4', (3,))
        if not np.array_equal(binary_grid, grid):
            raise ValueError('Field binary grid differs from its metadata')
        fields = read_array(stream, '<f8', (int(np.prod(grid)), 4)).copy()
        if stream.read(1):
            raise ValueError('Unexpected trailing field records')
    moments = np.loadtxt(directory/'moments.dat', ndmin=2)
    if moments.shape != (saved['nat'], 5) or not np.isfinite(moments).all():
        raise ValueError('Incomplete or nonfinite atomic moments')
    if not np.array_equal(moments[:, 0], np.arange(1, saved['nat'] + 1)):
        raise ValueError('Atomic moments are not in QE atom order')
    moments = moments[saved['ni'], 1:4].copy()
    if np.any(np.linalg.norm(moments, axis=1) < 1e-8):
        raise ValueError('A Ni moment vanishes; its V1 spin direction is undefined')
    return fields, moments


def write_package(stage: Path, saved: dict, branch: dict, fields: np.ndarray,
                  moments: np.ndarray, provenance: dict) -> dict:
    scalar_mean = float(fields[:, 0].mean())
    arrays = dict(scalar=fields[:, 0] - scalar_mean, B=fields[:, 1:], N=saved['N'],
                  moments=moments, atom_ids=saved['ni'] + 1, grid=saved['grid_dense'])
    for name, value in arrays.items():
        np.save(stage/(name + '.npy'), value, allow_pickle=False)
    np.savez(stage/'structure.npz', cell=saved['cell'], positions=saved['positions'],
             species=saved['species'], ni=saved['ni'],
             spins=moments / np.linalg.norm(moments, axis=1)[:, None])
    manifest = dict(schema='LANIO2_V1_SCF_LABELS_V1', source=provenance,
                    units=dict(scalar='Ry', B='Ry', N='dimensionless', moments='mu_B',
                               cell='angstrom', positions='angstrom'),
                    grid=saved['grid_dense'].tolist(), grid_order='x fastest; ix + nx*(iy + ny*iz)',
                    field_convention='v*I + Bx*sigma_x + By*sigma_y + Bz*sigma_z',
                    scalar_components='local ionic + Hartree + XC; full-cell mean subtracted',
                    scalar_removed_mean_Ry=scalar_mean,
                    B_components='XC Pauli coefficients; no explicit constraint field or Hubbard term',
                    moments_source='get_locals on the converged density; not target directions',
                    gga=branch, occupation=saved['metadata'], atom_ids=(saved['ni'] + 1).tolist(),
                    files={name: sha256(stage/name) for name in
                           [*(key + '.npy' for key in arrays), 'structure.npz']})
    (stage/'LABELS.json').write_text(json.dumps(manifest, indent=2, allow_nan=False) + '\n')
    return manifest


def _quoted(value: str) -> str:
    if '\n' in value or '\r' in value:
        raise ValueError('Newlines are not allowed in QE paths')
    return "'" + value.replace("'", "''") + "'"


def extract(save: Path, scf_output: Path, executable: Path, output: Path) -> dict:
    save, scf_output, executable = (Path(p).resolve() for p in (save, scf_output, executable))
    output = Path(output).absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(f'Refusing to replace existing output: {output}')
    saved = read_save(save)
    if save.name != saved['prefix'] + '.save':
        raise ValueError('Save directory name differs from the XML prefix')
    branch = read_scf_log(scf_output, saved)
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise ValueError('Provide the executable built from 4_potential/export_labels.f90')
    sources = {name: save/name for name in ('data-file-schema.xml', 'occup.txt', 'charge-density.dat')}
    sources.update(scf_output=scf_output, exporter=executable)
    before = {name: sha256(path) for name, path in sources.items()}
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.' + output.name + '-', dir=output.parent) as temporary:
        stage = Path(temporary)
        source_lsign = '.true.' if branch['lsign'] else '.false.'
        axis = ', '.join(format(v, '.17g') for v in branch['ux'])
        (stage/'extract.in').write_text('&inputpp\n prefix=' + _quoted(saved['prefix'])
            + ',\n outdir=' + _quoted(str(save.parent) + '/')
            + ',\n source_lsign=' + source_lsign + ',\n source_ux=' + axis + '\n/\n')
        try:
            with (stage/'extract.out').open('w') as log:
                subprocess.run([str(executable), '-in', 'extract.in'], cwd=stage, stdout=log,
                               stderr=subprocess.STDOUT, check=True, timeout=900)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
            tail = '\n'.join((stage/'extract.out').read_text(errors='replace').splitlines()[-20:])
            raise RuntimeError(f'QE label extraction failed; no package published.\n{tail}') from exc
        if 'V1_LABEL_EXPORT_DONE' not in (stage/'extract.out').read_text():
            raise ValueError('QE label exporter did not complete')
        fields, moments = read_fields(stage, saved, branch)
        if before != {name: sha256(path) for name, path in sources.items()}:
            raise ValueError('SCF inputs or exporter changed during extraction')
        provenance = {name: dict(path=str(path), sha256=before[name]) for name, path in sources.items()}
        result = write_package(stage, saved, branch, fields, moments, provenance)
        for name in ('V_LABELS.bin', 'moments.dat', 'LABELS.meta'):
            (stage/name).unlink()
        if output.exists() or output.is_symlink():
            raise FileExistsError(f'Output appeared during extraction: {output}')
        stage.rename(output)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--save', type=Path, required=True, help='Converged prefix.save directory')
    parser.add_argument('--scf-output', type=Path, required=True, help='Complete original pw.x SCF log')
    parser.add_argument('--executable', type=Path, required=True, help='Built export_labels.x (one MPI rank)')
    parser.add_argument('--output', type=Path, required=True, help='New external per-configuration label directory')
    args = parser.parse_args()
    extract(args.save, args.scf_output, args.executable, args.output)
    print(f'Wrote V1 v/B/N labels to {args.output}; no SCF or training performed.')


if __name__ == '__main__':
    main()
