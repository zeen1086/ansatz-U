# ansatz+U

Fortran tools for LaNiO₂ potential fitting, Wannier/RVR processing, Hamiltonians, bands, and torques. The directory layout follows the stages of [ansatz](https://github.com/zeen1086/ansatz).

| Directory | Contents |
|---|---|
| `0_QE` | QE 7.6 source modifications, installation, and build scripts |
| `1_scf`–`3_band_no_mag` | Structure, SCF, and band input templates |
| `4_potential` | Potential splitting, fitting, and grid conversion; `local/` contains frozen V1 inference sources |
| `5_wannier` | NSCF/Wannier templates and wavefunction/projector format conversion |
| `6_RVR` | Current V1 local spatial-kernel projection, merging, and background/projector readers |
| `7_spin` | V1 Hamiltonians, batch calculations, and complete-energy readout |
| `common`, `examples`, `tests` | Shared numerical modules, input templates, and offline tests |

## Current V1 workflow

At fixed geometry, the frozen iteration-30 V1 weights predict the full complex occupation once, `N_dir = sigmoid(Q0)`. The scalar and magnetic coefficients are `p = p0 + K(N_dir - N_ref)`; the Hamiltonian uses `VU = Ueff (I/2 - N_dir)`. Projecting its occupied states gives a separate output occupation. See the [native RVR guide](6_RVR/README.md) for commands.

```text
local_vfunc.py export / export-background
  → prepare-projection
  → local/rvr prepare/project/combine → B_general
  → local_vfunc.py cal_Hsk → local/cal_Hsk
  → 7_spin batch energy, torque, and J12 readout
```

Use an existing Python 3.11+ environment with NumPy/SciPy, gfortran, and BLAS/LAPACK. On macOS, the default library is Accelerate. To build only the current V1 tools:

```sh
make local-rvr BUILD=/absolute/external/build
```

The `--code` argument of `local_vfunc.py export` can point to this repository's `4_potential/local/`. Keep model parameters, Wannier data, RVR data, and calculation outputs at external paths. `export-background` generates pure H0; assembly uses `H0 + Hv + HB + HU`. The `.au` binary layout is unchanged: legacy KB/N0 placeholder fields are checked but do not enter the current calculation.

New calculations retain essential results by default. The unified `local/rvr prepare/project/combine` commands check and regenerate missing shards, then remove caches/shards after a successful merge. A completed single-configuration `local_vfunc.py cal_Hsk` run keeps matrices and summary JSON. Add `--keep-intermediates` to the corresponding command for debugging. Production coefficients and checkpoints needed for restart are retained; see the [local interface](6_RVR/README.md).

[7_spin](7_spin/README.md) provides same-kz batch assembly, RVR adjoint contraction, node-local slices, and complete energy/torque/J12 with a frozen complementary-space density. Use the [external JSON configuration template](examples/7_spin/production.json). Complete-energy and band-free-energy readouts have different definitions; implementation regression tests are not independent DFT validation.

## Retired interfaces and other stages

The old three-component RVR, scalar8/B40 `joint` RVR, `speedup`, HR rearrangement, and rigid-N0/tensor `ansatzU` entry points were removed on 2026-11-16. Their commands and templates are no longer supported on this branch. Sources remain available in the [pre-cleanup revision 3e1307a](https://github.com/zeen1086/ansatz-U/tree/3e1307ac721e91a66eb3851dae2df6bc4c9c99e1). Their formats and background definitions differ; renaming an old model does not make it compatible with current V1.

Build fitting tools separately: `make joint BUILD=/absolute/build` builds `fit_joint/eval_joint`; see [fit](4_potential/fit/README.md). Their parameters are not V1 parameters. `make legacy` builds only `fit_new`; `make all` builds V1 and these fitting tools. See [7_spin](7_spin/README.md) for single-configuration and batch usage.

## Shared ansatz workflow

Requirements: Python 3.11+, NumPy, Fortran, and BLAS/LAPACK, plus the appropriate MPI/Wannier programs for those stages. Use an existing Python interpreter and keep generated inputs and builds outside the repository:

```sh
export ANSATZ_BUILD_DIR=/absolute/external/ansatz-build
/path/to/existing/python configure.py --output-dir /absolute/external/case
/path/to/existing/python build.py 4_potential/V_split.f90 4_potential/make_R_vec.f90
```

Regenerate inputs and rebuild after changing `parameter`. `run.py` checks build signatures and preserves the calling directory. See `0_QE/README.md` for QE installation. Set pseudopotential names and numerical grids in the templates to match the actual inputs. Structure, projection, and lattice templates come from the original LaNiO₂ inputs and do not modify the frozen model.

Set `ANSATZ_PARAMETER=/absolute/external/parameter` to select external configuration. Wannier scripts preserve the calling directory: enter the external run directory before invoking them. Default repository paths are explicit placeholders and do not read neighboring historical directories.

`V_split` reads `polar2` and an explicit mapping file. Its first line is `n_spins n_channels n_cells`, followed by one line of 1-based spin indices per split cell. Cells follow x/y/z order, with z varying fastest; see `examples/4_potential/channel_map.dat`. Output retains the original `V1...` files and adds `spinN/V_HYC_dp` and `spin_info` for direct fitting input.

For labels containing periodic endpoints, the historical five-record fitting model supports `make_R_vec base_V_HYC_dp fitted_extended_model output_model output_R_vec` to remove endpoints while preserving all five records. Restart fitting with a model whose grid matches the labels. Current V1 uses its own export package; see `6_RVR/README.md` for spatial-kernel/Wannier grid, translation, and weight matching.

`check.py` runs offline tests. See `VALIDATION.md` for numerical summaries and `PROVENANCE.md` for source provenance. Supply models, wavefunctions, caches, and run outputs through external paths.
