# Source provenance

The frozen LaNiO₂ V1 evaluator is implemented in `4_potential/local/`. `local_vfunc.py` exports its parameters and runtime; `6_RVR` projects its fields, and `7_spin` assembles Hamiltonians and evaluates energies and torques. The original coupled V1 weights are compiled into a single explicit response without refitting.

New training in `4_potential/fit.py` directly uses that explicit-response model and its gradients, with trainable A and K and no old feedback trainer or warm model. References and normalization come from the validated label packages. Spatial kernels use the configured finite D4h radial-polynomial basis, with the existing native interpolation/window contract. This new-fit configuration is not the historical frozen-weight training protocol.

QE modifications and Wannier preparation tools originate from [zeen1086/ansatz](https://github.com/zeen1086/ansatz), commit `05fa36cbf6fb7740bd55293466e79fa13af545c7`. Their original copyright notices and the [QE license](../0_QE/License) are retained.

RVR uses the inverse-translation Wannier convention and the `+ikR` output of `5_wannier/make_wan/make_wan_fun_mpi.f90`. Shared binary readers are in `6_RVR/frozen_io.f90` and `common/npy_io.f90`. Model weights, reference data, and calculation results are external inputs; see [README](../README.md) and [validation](VALIDATION.md).

`4_potential/generate_configs.py` generates fresh configurations from a preceding SCF input. The full-angle uniform-polar sampling and integer upper-hemisphere sampling follow the rules recorded in `hubbard_training_candidates_20261102_r14`; a user-selected seed makes new samples repeatable. FM, relative-angle and Néel definitions were checked against archived SCF targets. The generator does not depend on historical random-angle files or converged magnetic moments. Pseudopotentials are supplied externally.
