# Magnetic SCF preparation and V1 labels

`generate_configs.py` creates the magnetic SCF inputs described in the root
README. After running and converging each SCF, `extract_labels.py` converts
its saved density and complete complex Hubbard occupations into one V1 label
package. It does not run SCF, fit a model, or submit jobs.

## Build the label exporter

Use the same configured and built QE 7.6 tree as the SCF calculation, including
the project's magnetic source modifications. Its pseudopotentials must remain
available to QE. The exporter links the existing QE libraries; it does not
modify `0_QE`, the QE source, or the saved calculation. All build products go
to `BUILD`:

```sh
make -C /absolute/qe-7.6/PP/src \
  -f /absolute/ansatz-U/4_potential/Makefile.labels \
  BUILD=/absolute/external/label-build
```

The exporter requires one MPI rank and reads no wavefunctions. Launch the
Python command as one process in an environment where the matching QE
executable and libraries can run. It invokes the exporter directly and allows
up to 15 minutes for one label extraction.

## Extract one converged configuration

```sh
/path/to/existing/python /absolute/ansatz-U/4_potential/extract_labels.py \
  --save /absolute/cases/fm/tmp/lanio2.save \
  --scf-output /absolute/cases/fm/scf.out \
  --executable /absolute/external/label-build/export_labels.x \
  --output /absolute/external/labels/fm
```

Repeat for each successfully converged configuration, using its own save and
original SCF log. The input save must contain `data-file-schema.xml`,
`charge-density.dat`, and `occup.txt`. The supported contract is noncollinear,
no-SOC PBE La4Ni4O8 with four Ni-3d Dudarev U=1 eV shells and ortho-atomic
projectors, in the registered tetragonal 2×2×1 Ni layout. The save directory
name must match its XML prefix.

The complete original `pw.x` log is required: its convergence, atom count and
final energy are checked against the save. The original `Fixed quantization
axis for GGA` message selects the signed GGA branch; when active, the axis is
restored at full precision from the XML initial magnetic direction and checked
against the log. The read-back defaults of QE post-processing are not used to
choose this branch. A concatenation of multiple PWSCF runs is rejected.

The exporter reconstructs Hartree and XC fields from the saved density,
restores the original magnetic GGA branch, disables the explicit constraint
field, and adds the local ionic potential to the scalar channel. It exports
the dense FFT grid without padding. Actual atomic moments are integrated from
that density with QE `get_locals`. Hubbard N is read separately from native
`occup.txt`: XML `Hubbard_ns_nc` contains absolute values and cannot preserve
its complex phases. The XML magnitudes are still checked against `occup.txt`
to reject inconsistent saved occupations.

## Output contract

| File | Contents |
|---|---|
| `scalar.npy` | float64 `(ngrid,)`, Ry; ionic + Hartree + XC scalar, with its full-cell mean subtracted |
| `B.npy` | float64 `(ngrid,3)`, Ry; XC coefficients Bx, By, Bz in `v I + B·sigma`, without the explicit constraint or Hubbard term |
| `N.npy` | complex128 `(4,10,10)`; complete dimensionless occupation matrix for each Ni |
| `moments.npy` | float64 `(4,3)`; actual density-integrated moments in Bohr magnetons |
| `atom_ids.npy` | One-based QE atom indices corresponding to N and moments |
| `grid.npy` | Dense-grid dimensions `(nx,ny,nz)`; flat index `ix + nx*(iy + ny*iz)` |
| `structure.npz` | `cell, positions` in angstrom, atomic numbers `species`, zero-based `ni`, and unit `spins` from the actual moments |
| `LABELS.json` | Units, channel/order conventions, removed scalar mean, GGA branch, occupation checks, source paths/hashes and output hashes |
| `extract.in`, `extract.out` | The post-processing input and log |

N uses spin-major order: five up orbitals followed by five down orbitals, with
`dz2, dxz, dyz, dx2-y2, dxy` in each block. The conversion is
`N[s1,s2] = transpose(ns_nc(:,:,2*s2+s1,atom))` for zero-based spin indices.
The scalar/B arrays retain the raw grid origin; no training-specific shift or
conversion to eV is applied. B is not mean-subtracted. No model parameters,
sample weights or fitting protocol are created here.

An existing output directory is never replaced. The package is published only
after validating convergence, N Hermiticity/bounds, grid/cell/GGA consistency,
finite fields and nonzero Ni moments. Source save files must remain unchanged
during extraction. These checks establish label-format consistency, not SCF
grid convergence, target-angle accuracy or fitting accuracy.

Keep pseudopotentials, saves, logs and generated label packages outside the
repository. Only code and documentation are distributed here.

## Offline checks

```sh
/path/to/existing/python -B /absolute/ansatz-U/4_potential/test_labels.py
```

The tests cover synthetic complex occupations, ordering, scalar gauge, GGA
branch matching and transactional publication. They also run through the
repository's `tests/check.py`. A real QE build and numerical extraction must
be checked on the host with the matching QE libraries.
