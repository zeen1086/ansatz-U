# ansatz+U

One supported workflow for the frozen LaNiO₂ V1 model: export existing weights, prepare Wannier/RVR inputs, and calculate Hamiltonians, energies, and torques. The repository does not provide a complete V1 training command. Model weights, reference densities, and Wannier data must be supplied separately.

| Directory | Contents |
|---|---|
| `0_QE` | QE 7.6 source modifications, installation, and build scripts |
| `1_scf`–`3_band_no_mag` | Structure, SCF, and band input templates |
| `4_potential` | Potential extraction helpers; `local/` contains the frozen V1 inference sources |
| `5_wannier` | NSCF/Wannier templates and wavefunction/projector format conversion |
| `6_RVR` | Current V1 local spatial-kernel projection, merging, and background/projector readers |
| `7_spin` | V1 Hamiltonians, batch calculations, and complete-energy readout |
| `common`, `examples`, `tests` | Shared numerical modules, input templates, and offline tests |

## Current V1 workflow

At fixed geometry, the frozen iteration-30 V1 weights use one explicit occupation response: `N0 = sigmoid(Q0)`, then `N1 = sigmoid(Q0 + A(N0 - N_ref))`, with the normalized maps compiled into `A = L K`. The scalar and magnetic coefficients are `p = p0 + K(N1 - N_ref)`; the Hamiltonian uses `VU = Ueff (I/2 - N1)`. There is no fixed-point iteration or refit. Projecting its occupied states gives a separate output occupation. See the [native RVR guide](6_RVR/README.md) for export compatibility and commands.

```text
local_vfunc.py export / export-background
  → prepare-projection
  → local/rvr prepare/project/combine → B_general
  → local_vfunc.py cal_Hsk → local/cal_Hsk
  → 7_spin energy, torque, and J12 or two-site K4 readout
```

Use an existing Python 3.11+ environment with NumPy/SciPy, gfortran, and BLAS/LAPACK. On macOS, the default library is Accelerate. The default build produces the two current executables, `local/rvr` and `local/cal_Hsk`:

```sh
make BUILD=/absolute/external/build
```

The `--code` argument of `local_vfunc.py export` can point to this repository's `4_potential/local/`. Keep model parameters, Wannier data, RVR data, and calculation outputs at external paths. `export-background` generates pure H0; assembly uses `H0 + Hv + HB + HU`. The `.au` binary layout includes reserved KB/N0 fields, which are checked but do not enter the calculation.

New calculations retain essential results by default. The unified `local/rvr prepare/project/combine` commands check and regenerate missing shards, then remove caches/shards after a successful merge. A completed single-configuration `local_vfunc.py cal_Hsk` run keeps matrices and summary JSON. Add `--keep-intermediates` to the corresponding command for debugging. Production coefficients and checkpoints needed for restart are retained; see the [local interface](6_RVR/README.md).

[7_spin](7_spin/README.md) provides batch assembly and energy/torque readout from the same model. Use the [external JSON configuration template](examples/7_spin/production.json). The band-free-energy and complete-energy readouts are distinct physical observables; complete energy includes a frozen complementary-space density and requires a matching QE 7.6 adapter and reference data. J12 and two-site K4 are supported complete-energy analyses. Implementation regression tests are not independent DFT validation.

## Choose your starting point

- **Existing exported model, RVR, and background:** start with the single-configuration command in [6_RVR](6_RVR/README.md#4-reuse-across-supercells-and-spin-configurations), then use [7_spin](7_spin/README.md) for energies and torques.
- **Existing frozen weights and Wannier data:** follow [6_RVR](6_RVR/README.md) from export through projection.
- **Preparing QE and Wannier inputs:** use the preparation tools below. These stages do not train the frozen V1 model.

## QE and Wannier preparation

Requirements: Python 3.11+, NumPy, Fortran, and BLAS/LAPACK, plus the appropriate MPI/Wannier programs for those stages. Use an existing Python interpreter and keep generated inputs and builds outside the repository:

```sh
export ANSATZ_BUILD_DIR=/absolute/external/ansatz-build
/path/to/existing/python configure.py --output-dir /absolute/external/case
/path/to/existing/python build.py 4_potential/V_split.f90
```

Regenerate inputs and rebuild after changing `parameter`. `run.py` checks build signatures and preserves the calling directory. See `0_QE/README.md` for QE installation. Set pseudopotential names and numerical grids in the templates to match the actual inputs. Structure, projection, and lattice templates come from the original LaNiO₂ inputs and do not modify the frozen model.

Set `ANSATZ_PARAMETER=/absolute/external/parameter` to select external preparation configuration. Wannier scripts preserve the calling directory: enter the external run directory before invoking them. Replace the explicit placeholder paths with your input paths.

`V_split` reads `polar2` and an explicit mapping file. Its first line is `n_spins n_channels n_cells`, followed by one line of 1-based spin indices per split cell. Cells follow x/y/z order, with z varying fastest; see `examples/4_potential/channel_map.dat`. Output retains the original `V1...` files and adds extracted `spinN/V_HYC_dp` and `spin_info`.

V1 projection takes model-aligned raw `wan_1 ... wan_30` and `weighted_translations.dat`; see [6_RVR](6_RVR/README.md) for grid, translation, and weight matching. Exported models must declare `N_single_response_to_v_B` and include their complete runtime. Export the original coupled V1 weights with the current command to create a compatible package.

`check.py` runs offline tests. See `VALIDATION.md` for numerical summaries and `PROVENANCE.md` for source provenance. Supply models, wavefunctions, caches, and run outputs through external paths.
