# Run Hamiltonian, energy and torque calculations

Run from the repository root with an existing Python 3.11+ environment containing NumPy/SciPy. These commands use the registered LaNiO₂ V1 inputs. Supply the model, RVR and background prepared in [6_RVR](../6_RVR/README.md). Keep data and outputs outside the repository.

## 1. Build and configure

```sh
export ANSATZ_PYTHON=/absolute/existing/python
make local-rvr BUILD=/absolute/build
cp examples/7_spin/production.json /absolute/run/production.json
export ANSATZ_RUN_CONFIG=/absolute/run/production.json
```

Create the external run directory first. Edit the copied configuration:

| Setting | Required value/data |
|---|---|
| `python`, `native_build` | Existing interpreter and build directory containing `local/cal_Hsk` |
| `interface`, `model_id` | Directory containing `V_func_tmp`; ID from its `V_FUNC.json` |
| `work` | Shared directory containing `RVR/B_general`, `RVR/PROJECTION.json`, `basis/pure_H0.au`, `basis/bloch.npy` |
| `references` | Matching Wannier reference package: `data.npz`, `k_*.npz`, `wr_*.npy` |
| `input_root` | Directory containing `inputs/Lxx/structure.npz` and `CASE.json` |
| `spin_work` | Band-energy and spin-result directory |
| `complete_energy.work`, `complete_energy.shared_work` | Separate complete-energy result and shared-cache directories |
| `complete_energy.reference_run` | Matching reference package with `results/BANDS.npz`, `results/FIELDS.npz`, `inputs/reference_actual/DENSITY.bin` |
| `complete_energy.pseudo_dir`, `complete_energy.native_binary` | Matching `La.upf`, `Ni.upf`, `O.upf`, and the QE adapter built below |
| `complete_energy.sizes` | Requested sizes from `4, 6, 8` |

Use matching geometry, basis, Ueff and scalar-energy gauge throughout. Band free energy assumes an orthonormal basis; complete trial-density free energy includes the overlap metric and fixed reference complement. Their energy and J outputs are different observables.

## 2. Single Hamiltonian

```sh
"$ANSATZ_PYTHON" local_vfunc.py cal_Hsk \
  --model /absolute/model/V_func_tmp --structure /absolute/L04_structure.npz \
  --rvr /absolute/RVR/B_general --background /absolute/basis/pure_H0.au \
  --kpoints /absolute/kpoints.dat --binary /absolute/build/local/cal_Hsk \
  --output /absolute/new/H_case
```

The output directory must be new. Outputs are per-k Hamiltonians and `HAMILTONIAN.json`. Add `--derivative-sites 1 2 --step 0.001` for derivatives.

## 3. Band free energy and torque

```sh
"$ANSATZ_PYTHON" 7_spin/run_material.py spectrum --L 4 --case 1 --mesh 2
"$ANSATZ_PYTHON" 7_spin/run_material.py torque --L 4 --case 1 --mesh 2
# Run every case defined in CASE.json before collecting:
"$ANSATZ_PYTHON" 7_spin/run_material.py collect --L 4 --mesh 2
```

Results are under `spin_work/J/`; each case's `TORQUE.json` contains `F_eV` and `torque_eV`. Set `V1_K_WORKERS`, `K_WORKERS` and BLAS threads to fit the allocated resources.

For separate kz workers, use a [task list](../examples/7_spin/batch_tasks.json). Each command below is a stage, not one unattended sequence:

```sh
"$ANSATZ_PYTHON" 7_spin/kz_production.py stage --z 0 --rvr /local/z0.bin
"$ANSATZ_PYTHON" 7_spin/kz_production.py prepare --z 0 --rvr /local/z0.bin --tasks /absolute/tasks.json
"$ANSATZ_PYTHON" 7_spin/kz_production.py spectrum --z 0 --rvr /local/z0.bin --tasks /absolute/tasks.json
"$ANSATZ_PYTHON" 7_spin/kz_production.py spectrum-finish --z 0 --rvr /local/z0.bin --tasks /absolute/tasks.json
"$ANSATZ_PYTHON" 7_spin/kz_production.py torque --z 0 --rvr /local/z0.bin --tasks /absolute/tasks.json
"$ANSATZ_PYTHON" 7_spin/kz_production.py torque-finish --z 0 --rvr /local/z0.bin --tasks /absolute/tasks.json
```

Run `stage`, `spectrum` and `torque` for each z=0..7, with a separate slice path per worker. Run `prepare` once first. Run `spectrum-finish` once after all spectra finish, then the torque workers, then `torque-finish` once. These commands do not submit jobs.

## 4. Complete energy, torque and J12/K4

Build the QE adapter against the existing matching QE 7.6 libraries:

```sh
mkdir -p /absolute/qe_adapter_build
make -C /absolute/qe/PP/src -f /absolute/repo/7_spin/Makefile.native \
  TASK_DIR=/absolute/qe_adapter_build native
```

This recipe uses Intel compiler flags; use the matching QE `make.inc` and libraries. Set `complete_energy.native_binary` to the resulting `native_large.x`.

Start eight concurrent workers, indexed 0..7, with the same configuration and run ID. Generate the ID once in the shared launcher and pass it to every worker:

```sh
export ANSATZ_RUN_ID="$("$ANSATZ_PYTHON" -c 'import uuid; print(uuid.uuid4().hex)')"
# Worker 0; launch the other seven with arguments 1 through 7:
"$ANSATZ_PYTHON" 7_spin/complete_energy.py 0
# After all eight workers finish:
"$ANSATZ_PYTHON" 7_spin/complete_readout.py summary
```

The default J12 calculation uses the supplied CASE definitions. For a two-site 4×4 curvature matrix K4, add `"curvature": {"sites": [1, 2]}` inside `complete_energy`; sites are 1-based native Ni indices. K4 uses `structure.npz`, generates its 25 probes, and does not need CASE.json. Optional `inner_step_rad` and `outer_steps_rad` set derivative steps; the latter must contain three positive values h, 2h, 4h. Use fresh result/cache directories when changing this setup.

Outputs: `Lxx/case_XX/ENERGY.json` for complete energies, `Lxx/RESULT.json` for J12/K4, and summary CSV files including `J12_by_size.csv` or `Lxx/K4_eV_per_rad2.csv`. K4 is a periodic-supercell curvature in the two sites' tangent coordinates; the accompanying single-site energy/torque consistency check is not J12.

## 5. Restart and checks

Stop all previous workers before restarting. Keep inputs, code, model and executables unchanged, then rerun the stages. Complete-energy workers must share a run ID; use a new ID to isolate old synchronization markers. Inspect `attempts/<ID>/FAILURE.json` and logs if workers wait. Checkpoints and caches are reused only when their implemented identity and payload checks pass. Changed inputs or source require new shared-cache and result directories; do not edit identity records to force reuse.

```sh
"$ANSATZ_PYTHON" 7_spin/complete_readout.py selfcheck
"$ANSATZ_PYTHON" -B tests/check.py
```

With the actual reference package and QE adapter available, run the density-tiling check:

```sh
"$ANSATZ_PYTHON" 7_spin/native_input.py check
```

It additionally reads `reference_run/results/native/reference_actual_ENERGY.dat` and `reference_actual_TRUE.bin`, and writes `complete_energy.shared_work/native_check/`, which must not already exist. Numerical checks and completion markers do not establish material accuracy or k/grid/size convergence.
