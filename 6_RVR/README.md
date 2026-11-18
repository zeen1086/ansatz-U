# V1 RVR Usage Guide

Run the following commands from the repository root. Replace example absolute paths with real paths; `ANSATZ_PYTHON` must point to an existing Python environment with NumPy and SciPy. Keep models, builds, and calculation outputs outside the repository.

## 1. Build and Prepare the Model

```sh
export ANSATZ_PYTHON=/absolute/existing/python
make local-rvr BUILD=/absolute/build
```

This produces `/absolute/build/local/rvr` and `/absolute/build/local/cal_Hsk`.

Skip the following two steps if an exported model and pure H0 background already exist:

```sh
"$ANSATZ_PYTHON" local_vfunc.py export \
  --source /absolute/frozen_V1/result/final \
  --code /absolute/frozen_V1/model_code \
  --structure /absolute/L02_structure.npz \
  --scalar-offset-eV=-7.187822007469299 \
  --output /absolute/model/V_func_tmp

"$ANSATZ_PYTHON" local_vfunc.py export-background \
  --source /absolute/DFT_only/primitive.bin \
  --output /absolute/basis/pure_H0.au
```

The scalar constant above applies to the current frozen V1. Export with the original four-Ni geometry. The structure NPZ contains `cell, positions, species, ni, spins`: geometry in Å, atomic numbers, zero-based Ni indices, and unit spin directions. The background must be a compatible DFT-only `primitive.bin` in the same basis as the projected Wannier functions.

## 2. Prepare Projection Inputs

Provide `wan_1 ... wan_30`: little-endian complex128 stream files in Fortran order with shape `(n_translation, n_grid)`, on the same grid as the model. The first line of `weighted_translations.dat` is the translation count; each subsequent line is `Rx Ry Rz weight`.

```sh
"$ANSATZ_PYTHON" local_vfunc.py prepare-projection \
  --model /absolute/model/V_func_tmp \
  --wan-dir /absolute/wan \
  --translations /absolute/weighted_translations.dat \
  --background /absolute/basis/pure_H0.au \
  --grid-weight 1 \
  --kz 0 0.125 0.25 0.375 0.5 0.625 0.75 0.875 \
  --output /absolute/RVR
```

`--grid-weight 1` applies to the existing discretely normalized inputs; other sources require their actual integration weight. The program does not automatically truncate Wannier support, add boundary weights, or renormalize. Use a new output directory.

## 3. Compute, Resume, and Read Outputs

```sh
/absolute/build/local/rvr prepare /absolute/RVR/projection.nml
/absolute/build/local/rvr project /absolute/RVR/projection.nml
/absolute/build/local/rvr combine /absolute/RVR/projection.nml
```

`prepare` creates the shared spatial-field cache; `project` computes orbital/kz shards; `combine` merges them into reusable `B_general`. Shard tasks can be selected using one-based orbital and kz indices:

```sh
# All kz values for orbital 1
/absolute/build/local/rvr project /absolute/RVR/projection.nml 1
# kz index 3 for orbital 1
/absolute/build/local/rvr project /absolute/RVR/projection.nml 1 3
# Retain intermediate files for debugging
/absolute/build/local/rvr combine /absolute/RVR/projection.nml --keep-intermediates
```

Different orbital/kz tasks may run concurrently, with only one writer per file. Combine after all projections finish. Do not use `mpirun` to launch multiple processes writing the same shard.

After interruption, rerun the same commands. Spatial caches are reused after header and mean checks; shards are reused after header/byte-count checks and first/last-block sampling. Missing shards are computed, and unfinished `.partial` files are regenerated from scratch. Format/byte-count errors, nonfinite first/last blocks, or identity mismatches raise errors; move the affected files elsewhere before retrying. If the final file was published but cleanup was interrupted, rerunning `combine` resumes cleanup.

Keep the model, geometry, Wannier functions, translations/weights, and preparation records unchanged during restart. If inputs change, prepare the projection again in a new output directory.

By default, caches and shards are removed only after successful combination; failures retain intermediate files. The final directory contains:

```text
B_general          # Reusable matrix elements
PROJECTION.json    # Input provenance and preparation record
projection.nml     # Projection configuration
kz.dat             # Projected kz values
```

`PROJECTION.json` is not a completion marker. Successful combination prints `B_GENERAL_COMPLETE`.

## 4. Reuse Across Supercells and Spin Configurations

For fixed primitive-cell geometry, different L×L×1 supercells and spin textures share the same `B_general`, model package, and `pure_H0.au`, while computing their own Hamiltonians. kx/ky may change; kz must be in the projected list. Changing the model, geometry, or Wannier basis requires reprojection.

```sh
"$ANSATZ_PYTHON" local_vfunc.py cal_Hsk \
  --model /absolute/model/V_func_tmp \
  --structure /absolute/L04_structure.npz \
  --rvr /absolute/RVR/B_general \
  --background /absolute/basis/pure_H0.au \
  --kpoints /absolute/kpoints.dat \
  --binary /absolute/build/local/cal_Hsk \
  --derivative-sites 1 2 --step 0.001 \
  --output /absolute/H_case
```

The k-point file starts with the point count, followed by `kx ky kz weight` on each line; weights sum to 1. Derivative sites use native one-based Ni numbering `x*L+y+1` and denote angular derivatives about the global y axis. Omit `--derivative-sites` to compute only H.

Default outputs are per-k `value.k000001.h`, requested `dsite_XXXX.k000001.h`, and `HAMILTONIAN.json`. Matrices are little-endian complex128 in Fortran order, in eV, with dimension `60L²`; index and k-point information is stored under `native_indices` in the JSON.

Add `--components` for H0/Hv/HB/HU components, or `--keep-intermediates` to retain coefficients, native inputs, and diagnostic files. Intermediate files are removed only after all results pass validation; failures retain them. Single-configuration `cal_Hsk` does not automatically resume an existing output directory.

See the [production usage guide](../7_spin/README.md) for batch energy, torque, and production checkpoints.

Completed `B_general` files receive only header, dimension, geometry, and byte-count checks, without another full scan. First/last-block shard sampling cannot detect all interior-block corruption. Spatial caches are no longer read back in full after writing.
