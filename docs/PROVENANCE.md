# Source provenance

The frozen LaNiO₂ V1 evaluator is implemented in `4_potential/local/`. `local_vfunc.py` exports its parameters and runtime; `6_RVR` projects its fields, and `7_spin` assembles Hamiltonians and evaluates energies and torques. The original coupled V1 weights are compiled into a single explicit response without refitting.

QE modifications and Wannier preparation tools originate from [zeen1086/ansatz](https://github.com/zeen1086/ansatz), commit `05fa36cbf6fb7740bd55293466e79fa13af545c7`. Their original copyright notices and the [QE license](../0_QE/License) are retained.

RVR uses the inverse-translation Wannier convention and the `+ikR` output of `5_wannier/make_wan/make_wan_fun_mpi.f90`. Shared binary readers are in `6_RVR/frozen_io.f90` and `common/npy_io.f90`. Model weights, reference data, and calculation results are external inputs; see [README](../README.md) and [validation](VALIDATION.md).

The [bundled pseudopotentials](../1_scf/pseudo/README.md) match the recorded SHA256 values for all 26 V1 training states. `1_scf/generate_configs.py` is a new input generator, checked against recovered SCF target angles from `lanio2_hubbard_g2_20261031_r03/cases/s0_fmz`, `lanio2_hubbard_controlled_training_20261103_r26/cases/rel78_*`, and the `Neel_t14` input archived by `local_feedback_real_20261112_r01`. It currently covers those six states only. The missing original targets are not inferred from converged magnetic moments, which can differ from the imposed directions.
