# V1 Hamiltonian and Energy Usage Guide

## Directory Layout and Data Placement

All Python and Fortran sources and build files reside directly in `7_spin/`. Runtime data is split into:

```text
7_spin/
├── *.py, *.f90, Makefile.native, README.md
├── shared/                  # Spin-independent for a fixed model/geometry/basis
│   ├── model/V_func_tmp/
│   ├── RVR/                 # B_general, PROJECTION.json, etc.
│   ├── basis/               # pure_H0.au, bloch.npy, etc.
│   ├── wannier/references/
│   ├── complete_reference/, pseudopotentials/
│   └── complete_energy/     # Fixed complement density, persistent wavefunction/operator cache
└── spins/                   # Configuration inputs, results, and task state
    ├── inputs/Lxx/          # structure.npz, CASE.json
    ├── H_case/              # Single-configuration H/dH
    ├── regression/          # H/band checks for specified reference spin configurations
    ├── J/Lxx/.../case_XX/   # Coefficients, eigenstates, chemical potential, band free energy/torques
    └── complete_energy/     # Configuration density, complete energy, torques, J, and checkpoints
```

Both data directories are ignored by Git. Link them to external storage or specify external absolute paths in the configuration; do not push large calculation data to GitHub. Use a new set of directories when the model, geometry, basis, or reference density changes. Matching filenames alone do not establish compatibility.

The complete-energy workflow persistently stores wavefunctions, overlap matrices and transforms, H0, and projectors in `shared/complete_energy/cache/` for reuse across spins. Structure files contain spin directions, so the entire files belong in `spins/`.

## Single-Configuration Hamiltonian

Prepare the model, `B_general`, and `pure_H0.au` following the [RVR usage guide](../6_RVR/README.md), then build and run from the repository root:

```sh
make local-rvr BUILD=/absolute/build
/absolute/python local_vfunc.py cal_Hsk \
  --model /absolute/model/V_func_tmp \
  --structure /absolute/L04_structure.npz \
  --rvr /absolute/RVR/B_general \
  --background /absolute/basis/pure_H0.au \
  --kpoints /absolute/kpoints.dat \
  --binary /absolute/build/local/cal_Hsk \
  --output /absolute/repo/7_spin/spins/H_case
```

`cal_Hsk.f90` assembles `H=H0+Hv+HB+HU`. Default outputs are per-k H and `HAMILTONIAN.json`. Add `--derivative-sites 1 2 --step 0.001` for angular derivatives, `--components` for components, or `--keep-intermediates` for debugging files. Use a new output directory.

## Batch Energy, Torque, and J

Fill in an external configuration as described below, then run locally or in a PBS compute job. Distinguish the readout definitions:

| Entry point | Energy and output |
|---|---|
| `run_material.py`, `kz_production.py` | Gaussian band free energy `Fband`; `F_eV` in `TORQUE.json`, with torque `−dFband/dθ` |
| `complete_energy.py` | Complete trial-density free energy with fixed complement density; `energy` in `ENERGY.json`, with torque `−dF/dθ` |
| `summarize_energy.py` | Summarize J12 and compare single-site curvature from complete energy and torque; this curvature is not J12 |

Both free energies use a fixed electron count, a common chemical potential across all k points, and Gaussian smearing `σ=0.01 Ry`:

```text
f = erfc((ε−μ)/σ)/2
Fband = Σk,n wk [ε f − σ exp(−((ε−μ)/σ)²)/(2√π)]
F = Fband − Σk wk Tr[P(H−H0)]
    + EH[ρtotal] + Exc[ρtotal] − Exc_core
    + ∫ n_total vion + Ueff/2 Tr[Ntotal−Ntotal²]
```

`P` is the active-space density matrix; `ρtotal/Ntotal` include the fixed reference complement. The input-potential subtraction includes v/B/HU; the Hubbard functional is added only once. Complete energy omits configuration-independent reference constants, so compare configuration differences only at the same size, reference data, and energy gauge. Energy is in eV per supercell, torque in eV/rad, and J in meV. Do not mix energies or J values from the two readouts.

Complete energy uses the Wannier overlap metric. The band-free-energy entry diagonalizes H directly and requires an orthonormal input basis. At nonorthogonal interpolation points, identical smearing parameters alone do not make the band terms of the two entries equivalent.

Native assembly represents the fixed scalar zero as `cI`. The complete-energy entry replaces it with `cS` in the folded basis, solves using `S^(-1/2)`, and saves the corrected H for torque evaluation. Use consistent scalar-gauge conventions when comparing results. Existing matching model coefficients may be supplied through `reuse_cases`.

## Two-Site Transverse Curvature (K4)

To select K4 instead of J12, add this object inside the external configuration's
`complete_energy` section (site numbers are 1-based **native Ni** indices):

```json
"curvature": {
  "sites": [1, 2],
  "inner_step_rad": 0.001,
  "outer_steps_rad": [0.004363323129985824, 0.008726646259971648, 0.017453292519943295]
}
```

Use fresh spin and shared complete-energy work directories. Supply the existing
`inputs/Lxx/structure.npz` geometry and background spins; K4 generates its own
25 cases and does not read `CASE.json`. The existing supported sizes remain
4, 6 and 8; `sites` must fit every configured size. Exactly two distinct sites
and three positive outer steps `h, 2h, 4h` are supported. Omitted steps use the
values above. The worker and summarizer entry points are unchanged. The band/J
drivers do not implement K4; their supported interface uses CASE definitions.

For each selected spin `e0`, the code fixes orthonormal tangent columns `T` at
the background and uses the sphere exponential chart
`e(eta) = cos(|eta|) e0 + sinc(|eta|/pi) T eta`.
Both outer probes and inner four-point derivatives use this **same** chart,
including at displaced points. The native derivative filenames are
`dsite_XXXX_a1.bin` and `dsite_XXXX_a2.bin`. For standalone evaluation, a chart
structure contains `chart_base_spins`, `chart_sites`, `chart_basis`, `chart_eta`
and the matching displaced `spins`; use `local_vfunc.py evaluate --tangent`
with `--derivative-sites` in the same order as `chart_sites`. Chart structures
are rejected by the global y-axis derivative mode. Coefficient checkpoints
bind the chart arrays, derivative mode, step and coordinate source.

The convention is `tau = -dF/deta`, `K = d²F/deta² = -d(tau)/deta`, with no
Néel-sign conversion or shell averaging. `F` is the same complete trial-density
Gaussian free energy as above, at fixed electron number. Every probe updates
the model fields, model occupations, eigenstates and active density. Other
spin directions and the reference complement density/occupation stay fixed.
This is not a constrained self-consistent DFT Hessian or a spin-relaxed Hessian.

`Lxx/RESULT.json` contains the raw, unsymmetrized 4x4 matrix
`K_eV_per_rad2`, with row/column order
`[(site1,1),(site1,2),(site2,1),(site2,2)]`, including both on-site blocks.
It also records `F0_eV`, `tau0_eV_per_rad`, the background and tangent bases,
all 25 energies/torques, per-step matrices, Richardson step indicators,
reciprocity residuals, and energy-versus-torque **diagonal** curvature errors.
The summarizer writes `K4_eV_per_rad2.csv` with the same order.
Coverage refers only to these four coordinates; no uncomputed sites are zeros.
Results are periodic-supercell derivatives, without image-multiplicity division.
Nonzero background torque must be retained when interpreting local curvature;
positive curvature alone does not establish a stationary stable state.

The analytic checks need no model payload, native binary or production job:

```sh
/absolute/existing/python -B tests/test_spin_curvature.py
```

They compare all four blocks with an independent Heisenberg Hessian on collinear
and nonstationary noncollinear backgrounds, including energy mixed differences,
step dependence, basis covariance, native ordering, and synthetic driver I/O.
These checks do not validate actual-material accuracy, k/size convergence, or
the QE/native production path. Those require separate data and benchmarks.

## Inputs and Build

Use an existing Python 3.11+ environment with NumPy and SciPy. Native local programs require Fortran and BLAS/LAPACK.
Complete energy additionally requires a compiled QE 7.6 environment, pseudopotentials, and reference densities compatible with the reference data.

```sh
make local-rvr BUILD=/absolute/build
# Example Intel/OpenMP production build:
make local-rvr BUILD=/absolute/intel-build FC=ifort \
  FFLAGS='-O2 -check bounds -free -extend-source -assume byterecl -qopenmp' \
  LIBS='-mkl=parallel'
```

See the [local RVR interface](../6_RVR/README.md) for V1 model export and compatibility. New exports package the single-response model, loader and predictor with the frozen parameters and the feature, symmetry and field kernels from `--code`. Packages must declare `N_single_response_to_v_B` and include the complete runtime. Export the original coupled V1 weights with the current command to create a compatible package. Use fresh spin and complete-energy work paths when model or runtime inputs change.

Standalone `local_vfunc.py evaluate` generates only the production-required `value.bin/dsite_XXXX.bin` and `PREDICTION.json` by default. `prediction.npz`, which duplicates these arrays, is generated only with `--keep-intermediates`. Batch drivers retain coefficients, indices, and completion markers needed for restart; they do not apply the temporary-file cleanup used after the single-configuration `cal_Hsk` command.

Copy the [configuration template](../examples/7_spin/production.json) to an external run directory and replace all absolute paths:

```sh
export ANSATZ_RUN_CONFIG=/absolute/run/production.json
```

Configuration paths and data conventions:

- `python`: an existing interpreter on the execution host; `native_build`: build root containing `local/cal_Hsk`.
- `interface`: directory containing `V_func_tmp/V_FUNC.json`; replace the all-zero placeholder `model_id` with the actual ID from `V_FUNC.json`.
- `input_root`: contains `inputs/Lxx/structure.npz` and `inputs/Lxx/CASE.json`.
  CASE retains the production pairing, angle, and site definitions; see the local interface for supercell and Ni ordering.
  Each entry in `angle_configurations` accepts an optional `global_spin_rotation_deg` rotation vector.
  Its direction is the Cartesian rotation axis and its magnitude is the active right-hand angle in degrees.
  For example, `[0, 180, 0]` rotates the entire case about (+y) after the local site perturbation.
  Supply the original spins in `structure.npz` and specify the global rotation once in the case.
  The complete-energy reader applies it to each frozen reference occupation as
  `(U ⊗ I5) N_bg (U ⊗ I5)†` and to the three magnetization components of the same reference density.
  Omitting the rotation preserves the original frame, including for locally tilted cases.
- `references`: common-Wannier `data.npz`, `k_0001.npz ...`, and `wr_0001.npy ...`.
- `work`: shared data directory `shared/`, containing `RVR/B_general`, `RVR/PROJECTION.json`,
  `basis/pure_H0.au`, and `basis/bloch.npy`.
- `rvr_slices`: optional `{"0":"/node-local/z0/B_general", ...}`. Cache paths must be visible on the corresponding worker
  nodes and match the full RVR model, basis, dimensions, and kz. If omitted, slices are extracted from the full source.
- `spin_work`: configuration-result root `spins/`, containing `J/` and `regression/`; `regression/RESULT.json` is an optional diagnostic report, not a prerequisite for calculation.

The `complete_energy` configuration additionally requires:

- `work`: spin-dependent complete-energy output directory `spins/complete_energy/`; use a new directory whenever the model, inputs, angular step, or grid changes.
- `shared_work`: shared reference-density, wavefunction, overlap-transform, and operator-cache directory `shared/complete_energy/`; use a new directory when the model or reference inputs change.
- `reference_run`: original 2×2 reference package containing `results/BANDS.npz`, `results/FIELDS.npz`,
  and `inputs/reference_actual/DENSITY.bin`. Native density-tiling regression also reads
  `results/native/reference_actual_ENERGY.dat` and `reference_actual_TRUE.bin`.
- `pseudo_dir`: `La.upf`, `Ni.upf`, and `O.upf` from the same source.
- `native_binary`: native QE adapter built below.
- `sizes`: calculation sizes; only 4, 6, and 8 are accepted. No automatic addition of sizes or sampling.
- `reuse_cases`: optional `{"4":"/existing/L04/mesh2_h0.001", ...}`, from which
  `case_00 ... case_10/coeff` are read. These coefficients must have a `CHECKPOINT.json` generated by this version,
  with structure/spin, model, derivative sites, step, and file checks matching the current calculation. External directories are read-only.
  Existing eigenstates may also be reused if complete and passing residual checks; otherwise they are rediagonalized.
  If omitted, complete V1 coefficients for these configurations are generated in the new output directory.

The LaNiO2 package uses a primitive-cell grid of 40×40×36, 30 spatial orbitals, full 8³ WS data,
21 active electrons/Ni and 41 total valence electrons/Ni. Its frozen model and QE reference use Ueff=1 eV.
The complete-energy driver reads `Ueff_eV` from the model's `V_FUNC.json` and checks the Hamiltonian
coefficient headers against it. A calculation at another Ueff requires matching model and reference data.
Build the native QE adapter from an existing QE `PP/src`, placing outputs in a separate directory:

```sh
mkdir -p /absolute/qe_adapter_build
make -C /absolute/qe/PP/src \
  -f /absolute/repo/7_spin/Makefile.native \
  TASK_DIR=/absolute/qe_adapter_build native
```

This make recipe uses the production Intel compiler flags `-heap-arrays -check bounds`. Use the same
QE `make.inc` and static libraries; do not mix module files from other versions.

## Batch Band-Free-Energy Readout

Run locally or in an allocated PBS compute job. `run_material.py` exposes the following stages:
assemble H by kz group, determine a common chemical potential across all k points, then compute HF torques.

```sh
/absolute/python 7_spin/run_material.py spectrum --L 4 --case 1 --mesh 2
/absolute/python 7_spin/run_material.py torque --L 4 --case 1 --mesh 2
/absolute/python 7_spin/run_material.py collect --L 4 --mesh 2
```

`collect` requires all configurations specified in CASE; it does not label a single-configuration output as complete J by default.
Missing, failed, or mismatched material-regression reports produce only a `QUALITY WARNING`; calculation continues. `V1_K_WORKERS` controls threads for batched native assembly within one kz,
while `K_WORKERS` controls parallel batch contraction. Configure mathematical-library threads for the actual node.

Local-slice calculations across configurations use explicit stages. The task file contains
`[L, mesh, step_rad, case]` lists as shown in the [batch example](../examples/7_spin/batch_tasks.json):

```sh
/absolute/python 7_spin/kz_production.py stage --z 0 --rvr /local/z0.bin
/absolute/python 7_spin/kz_production.py prepare --z 0 --rvr /local/z0.bin --tasks /run/tasks.json
/absolute/python 7_spin/kz_production.py spectrum --z 0 --rvr /local/z0.bin --tasks /run/tasks.json
```

Only one task runs `prepare`; run `spectrum` once for each z=0..7.
After all spectra finish, one task runs `spectrum-finish` to determine the common chemical potential. Then run `torque` for each z=0..7,
followed by one `torque-finish` task after all torques finish.
These stages do not submit or cancel jobs. Missing PBS_JOBID or ANSATZ_RUN_CONFIG produces only a NOTE; without the latter, the repository configuration template is read. Its actual data paths must still be accessible.

## Complete-Energy Readout

See the energy definition above. Complete energy uses overlap transforms at points outside the original sampling and computes density response at fixed electron count.

Each of eight workers is assigned one kz. The existing setup of 16 cores per node and node-local slices is recommended:

```sh
# Optional: set one ID in a common launcher and pass it to all eight workers.
export ANSATZ_RUN_ID="$(/absolute/python -c 'import uuid; print(uuid.uuid4().hex)')"
# Set z=0..7 across compute nodes, using the same ANSATZ_RUN_CONFIG and ANSATZ_RUN_ID.
/absolute/python 7_spin/complete_energy.py 0
# After all tasks finish:
/absolute/python 7_spin/summarize_energy.py
```

The calculation uses a fixed 2×2×8 k mesh, Gaussian smearing of 0.01 Ry, an inner five-point derivative step of 0.001 rad,
and outer angles of ±0.25°/±0.5°/±1°. It computes only nearest-neighbor J12 and reverse-pairing checks.
Outputs retain full densities, spectra, native energies, per-k torques, and directional J values.
Completion markers do not establish physical convergence with size or finer sampling.

Main results: band free energy is stored in each configuration's `TORQUE.json/F_eV`; complete free energy is stored in
`Lxx/case_XX/ENERGY.json`, and J12 in `Lxx/RESULT.json`.
Summarization generates `J12_by_size.csv`, per-size `energy_components.csv`,
and `ENERGY_TORQUE_CHECK.json`. The latter compares single-site energy curvature with the torque derivative; it is not J12.

### Restart After Interruption

1. Confirm that all previous workers have stopped. Keep inputs, configuration, model, and executables unchanged.
2. Rerun the appropriate stage commands for band free energy. For complete energy, if `ANSATZ_RUN_ID` is unset, an ID is generated and reused under a lock at `run_ids/<key>`, keyed by `PBS_JOBID` (parent-process PID locally). Eight workers from one launcher share an ID; workers spanning different PBS jobs or parent processes must explicitly receive the same ID.
3. Complete, input-bound checkpoints are reused directly. Corrupt/missing spectra, densities, or per-k torques are recomputed. Unfinished coefficient, Hamiltonian, or QE directories are retained as `*.incomplete-*` and regenerated.
4. Inspect `attempts/<ANSATZ_RUN_ID>/FAILURE.json` and job logs. Failure markers print only a NOTE; waiting has no hard timeout. Start markers with the same ID may be overwritten, while flock prevents concurrent active workers. Explicitly use a new ID to isolate old barriers; reusing an ID does not automatically clear them.

Result files are written to temporary files and atomically replaced. QE stages use fresh temporary directories and publish after output validation. If a native stage finished but its energy JSON did not, restart writes the missing summary directly. The band-free-energy workflow removes temporary data such as eigenstates only after the whole configuration's `TORQUE.json` is published.

`INPUTS.json` binds structure/spin, CASE, model, basis, grid/step, and relevant code. Coefficients and native outputs also record input and file checks. Existing directories without `INPUTS.json` generally print a NOTE and register the current identity. A complete-energy directory containing existing configuration results must have `INPUTS.json` to resume. An existing manifest with a different identity still raises an error; use the corresponding inputs or a new directory. File identities use sizes/timestamps, while structure/spin use actual array descriptions; file-content digests are neither recorded nor checked. Old hash identities that differ from new records are not silently overwritten.

Retained `*.incomplete-*` directories and hidden temporary directories left by forcibly terminated processes are not reusable checkpoints. They may be removed after confirming that no writers are active; normal restart does not require deleting them first. Keep inputs and caches unchanged during execution.

Set `spin_work` explicitly for spin-dependent outputs. Complete energy also requires explicit `shared_work` and `work` directories for shared caches and current-run results.

## Small Checks and Native Regression

```sh
/absolute/python 7_spin/complete_readout.py
/absolute/python -B check.py
```

`complete_readout.py` includes finite-difference checks of energy, the Hubbard term, and fixed-electron-number response.
`density_map.py` retains an independent reference-projection algorithm used by complete-energy initialization and offline tests.
After configuring reference density and the QE adapter and allocating a PBS compute job, run the native density-tiling check:

```sh
/absolute/python 7_spin/check_large_native.py
```

Results are written to `complete_energy.shared_work/native_check/`, which must not already exist.
These checks do not replace material convergence or independent DFT accuracy validation.

## Shared Cache for Complete Energy

First rebuild with `make local-rvr BUILD=/absolute/build` and set
`complete_energy.shared_work` in the configuration (recommended: `7_spin/shared/complete_energy`). Running
`complete_energy.py 0` through `complete_energy.py 7` automatically prepares and reuses the cache; no separate cache command is required.

```text
shared/complete_energy/
├── INPUTS.json                 # Input and calculation-code identity record
├── reference/                  # Fixed reference complement
└── cache/
    ├── Lxx/kXXXX/
    │   ├── wave.npy            # Interpolated wavefunctions
    │   ├── overlap.npy         # S
    │   ├── transform.npy       # S^(-1/2)
    │   ├── metric.json
    │   └── CACHE.json
    ├── Lxx/operators_<k>/
    │   ├── h0.npy              # H0 in the folded basis
    │   ├── projectors.npy      # Projectors in the folded basis
    │   ├── native_projectors.npy
    │   ├── native.bin          # H0/projectors for native H assembly
    │   └── CACHE.json
    └── z*/RVR.bin              # Generated if no external slice is specified
```

Entries are built by L and k on first use. Subsequent spin spectra reuse overlap transforms; density and response readouts reuse
H0/projectors. Native H assembly and torque contraction also read the cache for the corresponding basis. H, model coefficients,
eigenstates, chemical potential, density, energy, and torque are still recomputed per spin and stored in `spins/`.
Single-configuration `local_vfunc.py cal_Hsk` and standalone band-free-energy entries retain their original direct assembly.

Input identities use file size, mtime/ctime, and existing model/basis IDs; `cache_identity()` no longer reads entire inputs to compute hashes. Shared-cache reuse checks identities, file sets, and byte counts without retaining content-digest fields, empty fields, or digest lists. Model-package `files` records file lengths; `verify` checks existence and length. New model/basis exports use random 64-character IDs, with unchanged native binary-header length.

New `.npy` caches of at most 256 MiB are checked row by row; larger arrays receive only a first-row check. Keep inputs and caches unchanged during execution. These checks cannot detect all same-length content corruption.

Caches use per-entry file locks, write to temporary directories, and publish atomically as a whole. Interrupted
`*.partial-*` directories are not read; retries regenerate unpublished entries. Remove leftover temporary directories only after confirming
that no writers remain. Missing files, byte-count mismatches, or identity mismatches in published caches raise errors
rather than being silently overwritten. Archive damaged entries before rebuilding with the same inputs; use new
`shared_work` and configuration-output directories when inputs or source code change. Existing reference/cache directories without `INPUTS.json` print a NOTE and register the current identity; an existing mismatched identity still raises an error.

Caches are not automatically deleted after calculation. Wavefunctions alone for 32 k points occupy approximately
`32 × L² × 57600 × 30 × 16` bytes: about 14.2/31.9/56.6 GB for L=4/6/8,
plus overlap and operator matrices. Available space is checked before building wavefunctions; concurrent tasks still require
space for the aggregate total. When only spin changes, reuse the same cache without copying it for each configuration.

## Runtime Checks

Material regression, finite electron-count deviations, and small positive overlap eigenvalues produce quality notes without blocking execution. Electron-count deviations are written to results or logs; completion does not mean these diagnostics passed.

Missing or damaged files, dimension/format errors, model/cache identity mismatches, nonfinite values, and non-positive-definite overlap matrices remain functional errors. Fix the inputs or regenerate the affected results; invalid caches cannot serve as valid restart checkpoints.

Hamiltonian Hermiticity residuals above 1e-5 and reference-projection residuals above 1e-6 produce only diagnostic notes. Adjoint/direct derivative residuals above 1e-6 still raise errors.

Re-exporting models or rebuilding projections generates new IDs. Update the configuration and use the corresponding new RVR; do not mix files from different exports. Generated records use declared IDs and file metadata for input identity.
