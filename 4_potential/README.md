# Magnetic SCF labels and single-response fitting

`generate_configs.py` creates the magnetic SCF inputs described in the root
README. After running and converging each SCF, `extract_labels.py` converts
its saved density and complete complex Hubbard occupations into one V1 label
package. It does not run SCF, fit a model, or submit jobs.
`fit.py` then learns a new model from those packages using the current explicit
occupation correction and produces a bundle accepted by `local_vfunc.py export`.

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

## Fit the single-response model

Copy `fit.toml` outside the repository, list the label directories and their
positive sample weights, and run with an existing NumPy/SciPy environment:

```sh
/path/to/existing/python /absolute/ansatz-U/4_potential/fit.py \
  --config /absolute/external/fit.toml \
  --output /absolute/external/new-fit
```

Only the current model is trained:

```text
N0 = sigmoid(Q0)
N1 = sigmoid(Q0 + A (N0 - Nref))
p  = p0 + K (N1 - Nref)
VU = Ueff (I/2 - N1)
```

Each prediction evaluates two matrix sigmoids and one explicit correction.
The optimizer learns Q/p heads, their biases, A, K and spatial-kernel parameters.
It does not load an old model, require warm weights, optimize L, or solve a
feedback fixed point. Nref comes from the training occupations projected onto
the registered spatial symmetry and actual spin directions. Feature RMS scales
are fixed from the labels. Only the initial matrix logit clips reference
eigenvalues near 0/1; occupation labels and predictions are not clipped.

The training loader checks label hashes, units, full complex N, actual moment
directions, Ni order and common geometry/grid. Scalar/B labels are converted
from Ry to eV once for the compact full-grid objective. With normalized sample
weights w, the loss is

```text
sum_t w_t [ mean_grid(dv²) / sv²
          + mean_grid(|dB|²) / sB²
          + ||dN||_F² / (100*nNi*sN²) ]
  + regularization * mean(all_optimizer_parameters²)
```

The scales and optimizer limits are explicit in `fit.toml`. B uses vector MSE;
the report also gives component RMS (vector RMS divided by sqrt(3)). Complex N
uses its full Frobenius norm, with no extra off-diagonal multiplier. All listed
samples are used for training; the report makes no held-out accuracy claim.

Spatial kernels use a 6 Å support and the native cosine window. Each radial
length supplies four nodal basis functions, proportional to
`[1, z², x²+y², x²-y²]`, with normalized coordinates. The default six lengths
give 24 A1g/B1g channels. The fit learns linear mixtures within the same spatial
representation and a scalar bias; it stays within this finite analytic spatial
span. This is not an unrestricted free-grid kernel fit or a reproduction of
the historical frozen kernels. Material accuracy and the adequacy of this
span require checking the resulting v/B/N residuals. Native interpolation and
window conventions are identical in training and export. Ni centers must lie
on dense-grid nodes, as they do for the registered 80×80×36 grid.

The output directory contains:

| Output | Meaning |
|---|---|
| `RUN.json` | Bound configuration, source-code identities and label hashes |
| `history.jsonl` | Accepted optimizer iterations and losses |
| `checkpoint.npz` | Atomic parameter checkpoint with schema and content checks |
| `model_XXXXXX/` | Complete model/spatial-kernel bundle plus `structure.npz` and `FIT_REPORT.json` |
| `RESULT.json` | Final bundle path, per-state RMS errors, optimizer status and gradient convergence |

A valid bundle may be saved when the iteration budget is reached; this is not
reported as optimizer convergence. CLI exit code 2 indicates that the optimizer
did not report success. Export the bundle identified by `RESULT.json`:

```sh
/path/to/existing/python /absolute/ansatz-U/local_vfunc.py export \
  --source /absolute/external/new-fit/model_000200 \
  --code /absolute/ansatz-U/4_potential/local \
  --structure /absolute/external/new-fit/model_000200/structure.npz \
  --scalar-offset-eV 0 \
  --output /absolute/external/new-model
```

Use the actual bundle name; `000200` is only an example. Zero scalar offset in
this example keeps the fitted mean-zero gauge. Set a physically matched scalar
offset before combining with an absolute-energy background; fitting does not
determine that constant.

To resume, repeat the fit command with `--resume`. Labels, settings and code
must be unchanged, and coupled/feedback checkpoints are rejected. Each call
allows another configured iteration budget. The saved parameters are restored;
L-BFGS history is restarted. History beyond the last saved checkpoint is rolled
back after interruption. Run only one optimizer per output directory. New
labels or settings require a new directory. Existing model bundles are retained.

The implementation is in `fit.py`, `training_data.py`, and `local/kernels.py`;
the existing `local/model.py` and `local/workflow.py` remain the prediction,
gradient and model-storage implementation.

## Offline checks

```sh
/path/to/existing/python -B /absolute/ansatz-U/4_potential/test_labels.py
```

The repository's `tests/check.py` also runs label checks, dataset checks,
joint-loss finite differences, short synthetic fits and resume/export checks.
The kernel test compiles the existing native evaluator when gfortran is
available and compares its interpolation and window values with Python.
These checks do not establish material fitting accuracy. A real QE build and
numerical extraction must be checked on the host with the matching QE libraries.
