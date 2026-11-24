# Four-spin, three-potential fitting

The production chain uses only Fortran, shell, and BLAS/LAPACK:

```sh
make joint BUILD=/absolute/build
4_potential/fit/run_joint.sh /absolute/states.tsv \
  examples/4_potential/joint_fit.nml /absolute/new_run /absolute/build
# Save the repository path before entering the external evaluation directory:
repo="$PWD"
cd /absolute/new_run
/absolute/build/eval_joint "$repo/examples/4_potential/joint_eval.nml" \
  V_func_tmp "$repo/examples/4_potential/joint_texture.dat"
```

`RUN` and `BUILD` must be outside the source tree, and the run directory must not already exist. Existing outputs are not overwritten. There is no Python production dependency, automatic installation, DFT/Wannier rerun, or job submission.

## Inputs and readers

Each line of `states.tsv` has five TAB-separated fields: `name weight anchor unit archive_path`. Paths are absolute or relative to the list file; `unit=Ry/eV` and weights are nonnegative. `anchor=1` imposes an exact constraint; `weight=0,anchor=0` defines a held-out state whose residuals are reported only. See `examples/4_potential/joint_archives.tsv`.

The shell uses `unzip -p` to extract named members from existing DFT-label NPZ files without computing or converting numerical values:

| Member | dtype / shape / convention |
|---|---|
| `grid.npy` | int32/int64, `(3,)`; x varies fastest in the flattened real-space grid |
| `atom_ids.npy` | int32/int64, `(4,)`; Ni5,6,7,8 by default |
| `v_scalar_Ry.npy`, `B_material_Ry.npy` | float64, `(ngrid,)` and `(ngrid,3)` |
| Corresponding members for eV labels | `v_scalar_eV.npy`, `B_material_eV.npy` |
| `moments_muB.npy` | float64, `(4,3)`; direction normalization is explicit in configuration, and original amplitudes are reported |
| `N.npy` | complex128, `(4,10,10)`; `[up d1..5,down d1..5]` for each Ni |
| Optional `V_U_eV.npy` | Same shape as N; checks `Ueff*(I/2-N)` without adding it to scalar/B labels |

Fortran reads NPY v1/v2/v3, C/Fortran array order, and little-endian data natively. Object arrays, single precision, truncated files, dimension/atom-order mismatches, nonfinite values, and non-Hermitian N are rejected.

Normalized spins are fixed-length directions. DFT amplitudes are already represented in the labels and fitted coefficients: there is no additional factor of 0.372, 0.4, or S² and no independent longitudinal-amplitude degree of freedom.

You can also call `fit_joint fit.nml dataset.manifest [V_func_tmp [fit_report.txt]]` directly. If the last two arguments are omitted, it creates `V_func_tmp` and `fit_report.txt` in the calling directory; explicit external paths are also accepted. The manifest starts with the number of states, followed by two lines per state: `name weight anchor unit` and the explicit extraction directory. Each directory contains `grid/atom_ids/scalar/B/moments/N.npy` and optional `VU.npy`. The shell writes `dataset.manifest` directly, without generating `rows.manifest`; existing run results are unchanged.

## Model and fitting objective

- Scalar: the original scalar8 features—one constant, six `ei·ej-1` terms, and `alpha*(1-alpha)`.
- B: the original 12 base vectors and 28 nonlinear residual vectors, fitting all three components together.
- N: the original `scalar/chi/vector/cross` features with real-symmetric/imaginary-antisymmetric d-orbital matrix bases reconstruct the full complex 10×10 matrix. All four Ni sites depend on all four input spins, allowing occupation eigenvalues to change.
- The original 64 square-cell crystal operations, Ni ordering, and QE real-d-orbital convention are retained. The features preserve global spin SO(3)/SU(2) and time-reversal properties.
- The grid requires `nx=ny` divisible by 4 and even `nz`, with Ni at the original quarter-cell positions. The four-input nonrigid N model is not automatically extended to arbitrary supercells.

The objective consists of three least-squares blocks for the raw scalar field, B field, and N matrices, with explicit state weights. It uses the old Python feature families; **it does not reproduce the old joint projected-operator loss or its frozen training results**.

SVD solves the problem directly, with hard anchors eliminated through a null-space construction. Reports include the effective relative cutoff, ranks, singular values, and per-state residuals. `fit_report.txt` collects fitting parameters, state weights and anchors, ranks/singular values, residuals, and occupation checks. Unidentified directions use the declared SVD minimum-norm solution; no regularizer or weight adjustment is added silently.

`remove_scalar_mean` is off by default. When enabled, the removed mean for each state is recorded explicitly; results cannot be compared directly with another potential zero. The original `fit_new` diagonal/full-component iteration and restart interfaces remain available. The new SVD entry point uses neither iteration counts nor learning rates.

The new `V_func_tmp` (formerly `model.joint`) retains its schema: 16-byte `LANIO2_JOINT_V1 `, twelve int32 values, two float64 values, two int32 flags, then column-major float64 arrays of shapes `(8,ngrid)`, `(40,ngrid)`, `(8,15)`, `(4,10)`, `(12,15)`, and `(6,10)`. Scalar/B coefficients are in eV; N coefficients are dimensionless. If training occupations are out of bounds, the candidate and report are retained with gate=0, and the evaluator rejects the candidate. N is not clipped.

Readers still check magic, version, byte order, dimensions, and total byte count. Existing `model.joint` contents need no conversion; only the default filename changed.

## Evaluation and angular derivatives

Input consists of four lines `ex ey ez axis_x axis_y axis_z` in DFT Ni order. Both directions and axes must be unit vectors. `derivative_site=1..4` selects the perturbed input. Analytic derivatives use `dei/dtheta=axis×ei` and include scalar/B terms and N on every affected Ni site.

`ueff` is read from the evaluation input: `VU=ueff*(I/2-N)` and `dVU=-ueff*dN`. Changing it does not solve a new DFT density.

Without a fourth evaluation argument, the following files and `eval_report.txt` are written in the calling directory. Binary layouts are preserved without splitting records:

| Former name | Current name | Units |
|---|---|---|
| `fields_Ry.bin` | `V_fit` | Ry; scalar/Bx/By/Bz |
| `dfields_Ry_per_rad.bin` | `dV_fit` | Ry/rad; same components |
| `N_VU_eV.bin` | `N_VU_fit` | Dimensionless N; VU in eV |
| `dN_dVU_eV_per_rad.bin` | `dN_dVU_fit` | dN in 1/rad; dVU in eV/rad |

The optional fourth argument is an `output_prefix`, such as `/absolute/new_run/probe`, producing four files such as `probe.V_fit` plus `probe.txt`. Field files store the grid in the first record and `(ngrid,4)` float64 values in the second, with x varying fastest. N files begin with atom order and Ueff; derivative files begin with atom order and perturbed Ni. N/VU and their derivatives each retain two `(10,10,4)` complex128 records, including all orbital/spin off-diagonal elements.

This entry point covers DFT labels → three-potential fitting → potentials and analytic derivatives. The new `V_func_tmp` is **incompatible** with the old five-record file of the same name and the local `.au` format; it cannot restart old `fit_new` calculations.

This scalar8/B40 fitting and grid-evaluation interface remains available, but its dedicated `joint` RVR and `ansatzU` backends were retired on 2026-11-16. Current [V1 RVR](../../6_RVR/README.md) uses frozen local kernels and the one-way occupation-to-field model. The old backend remains in the [pre-cleanup revision](https://github.com/zeen1086/ansatz-U/tree/3e1307ac721e91a66eb3851dae2df6bc4c9c99e1/6_RVR/joint).
