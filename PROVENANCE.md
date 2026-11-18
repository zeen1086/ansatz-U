# Provenance

2026-11-16: `6_RVR` was consolidated around the current V1 workflow. The legacy three-component, scalar8/B40, and rigid-N0 entries below refer to the [pre-cleanup version 3e1307a](https://github.com/zeen1086/ansatz-U/tree/3e1307ac721e91a66eb3851dae2df6bc4c9c99e1), not the current executable interfaces. See [6_RVR](6_RVR/README.md) for current files and dependencies. `background_projectors` and `project_at` were moved unchanged from `local_ansatz` to `frozen_io`; the `.au` format and checks remain compatible.

Common tools originate from [zeen1086/ansatz](https://github.com/zeen1086/ansatz), commit `05fa36cbf6fb7740bd55293466e79fa13af545c7`. QE modifications, Wannier, RVR, and build tools were retained. Material-specific fitting, symmetry projection, and spin-channel branches were removed; input templates use the original LaNiO₂ structure and projections. The QE license is retained at `0_QE/License`.

Source mapping for the LaNiO₂ branch:

| Original working-tree file | Historical import location |
|---|---|
| `fit/fit_new.f90` | `4_potential/fit/fit_new.f90`, `fit_mask.f90` |
| `1_make_RVR/cal_RVR_split.f90` | `6_RVR/cal_RVR_split_lanio2.f90` |
| `1_make_RVR/combine.f90` | `6_RVR/combine_lanio2.f90` |
| `9_super_wan/rearrange_HR_new.f90` | `6_RVR/rearrange_HR_lanio2.f90` |
| `9_super_wan/H_sk/cal_Hsk_new.f90` | `7_spin/cal_Hsk_lanio2.f90` |
| `9_super_wan_AFM/H_sk/cal_Hsk_new.f90` | `7_spin/AFM/cal_Hsk_lanio2.f90` |
| `9_super_wan/combine.f90` | `7_spin/combine_lanio2.f90` |
| `9_super_wan/band.f90` | `7_spin/band_lanio2.f90` |
| `9_super_wan/calc_tot_ene.f90` | `7_spin/calc_tot_ene_lanio2.f90` |

The new `frozen_io`, `tensor_rvr`, `local_ansatz`, `numerics`, and `run_ansatzU` port and extend the existing LaNiO₂ local Python reference. Full differences, original source identities, and numerical records are stored in external archives.

The four-spin, three-potential `joint_features/joint_model` interface implements scalar8, B12+28, complex N, and 64 crystal operations from the existing `covariant_basis.py` and `b_features.py`. `joint_npy` reads existing DFT labels natively; `joint_svd/fit_joint` implements a new weighted SVD objective; `eval_joint` outputs the three potentials and analytical angular derivatives. `run_joint.sh` only unpacks, builds, and launches; the production chain does not invoke Python. This objective follows the user's choice and does not inherit numerical validation or training results from the old operator loss.

The default three-potential parameter name changed from `model.joint` to `V_func_tmp`; evaluation files use `V_fit/dV_fit/N_VU_fit/dN_dVU_fit`. The new format is incompatible with the old five-record model of the same name.

Subsequent three-potential stages retain the original names, using `joint/` subdirectories to distinguish the new format:

| Original stage/shared function | New entry or extension | Purpose |
|---|---|---|
| `cal_RVR_split` | `6_RVR/joint/cal_RVR_split.f90`, `joint_rvr.f90` | Retain inverse-translation Wannier indexing and grid contraction; read the new V_func_tmp and project scalar8/B40/fixed scalar-reference kernels |
| `combine` | `6_RVR/joint/combine.f90` | Check row bindings and Hermiticity; generate the new B_general while retaining the original background/projectors and N coefficients |
| `local_ansatz::assemble` | Shared `background_projectors` | Extract the original H0+fixed Hv and p_R interpolation/periodic phases without changing numerical definitions; the original local interface calls the same function |
| `cal_Hsk` | `7_spin/joint/cal_Hsk.f90` | Assemble scalar-potential differences, B, and all four Ni HU terms; read out three torque contributions, bands/EF, and Φ from the same H |
| `run_ansatzU` | Explicit `joint_rvr_path` added | Select the original four-spin 2×2-cell backend; reject duplicate background/B inputs or implicit large-cell N extension |

The primitive30 read mapping follows the `+ikR` output of `5_wannier/make_wan/make_wan_fun_mpi.f90` and the original RVR inverse-translation convention directly, without re-Wannierization or an added amplitude factor. Production numerics in the new stages use Fortran/BLAS/LAPACK; the original fitting features and weights are unchanged.
