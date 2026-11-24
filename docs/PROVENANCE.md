# Source provenance

The frozen LaNiO₂ V1 evaluator is implemented in `4_potential/local/`. `local_vfunc.py` exports its parameters and runtime; `6_RVR` projects its fields, and `7_spin` assembles Hamiltonians and evaluates energies and torques. The original coupled V1 weights are compiled into a single explicit response without refitting.

QE modifications and Wannier preparation tools originate from [zeen1086/ansatz](https://github.com/zeen1086/ansatz), commit `05fa36cbf6fb7740bd55293466e79fa13af545c7`. Their original copyright notices and the [QE license](../0_QE/License) are retained.

RVR uses the inverse-translation Wannier convention and the `+ikR` output of `5_wannier/make_wan/make_wan_fun_mpi.f90`. Shared binary readers are in `6_RVR/frozen_io.f90` and `common/npy_io.f90`. Model weights, reference data, and calculation results are external inputs; see [README](../README.md) and [validation](VALIDATION.md).
