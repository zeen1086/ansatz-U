# V1 Validation Summary

Verified on 2026-11-24. The evaluator uses `N0 = sigmoid(Q0)`, `N1 = sigmoid(Q0 + A(N0-N_ref))`, and `p = p0 + K(N1-N_ref)`. Original coupled V1 weights are normalized and compiled into `A=L K` without refitting.

| Check | Result |
|---|---|
| Repository offline suite | 75 passed; 1 MPI check skipped because the compiler and mpirun are unavailable |
| DFT input preparation | Three UPF hashes match all 26 archived training-state records; six recovered target-angle sets match archived SCF inputs; primitive-to-supercell generation, preserved parameters, isolated scratch, and invalid-input rejection passed |
| Root layout and moved preparation commands | Offline suite passed after relocation; external-directory configure/build/run, argument forwarding, working-directory preservation, and stale-build rejection passed |
| Clean external builds through default `make` and `tools/build.py` | Passed; Python build fingerprints verified |
| Single-response algebra, full complex N, symmetry, occupation bounds, and parameter derivatives | Passed |
| Actual iteration-30 source conversion | 14 parameter arrays identical to the reference conversion |
| Native synthetic RVR/H/dH, L=1,2,3, independent real-grid integration | 120 component comparisons; maximum error 1.353e-16 eV |
| Cached versus direct native assembly | 120 component comparisons; maximum error 2.269e-16 eV |
| Batch versus single native assembly | 40 component comparisons; identical results |
| Actual V1 export → packaged runtime → synthetic RVR → native H/dH | Maximum error 9.155e-16 eV; all 15 compact/debug matrices identical |
| Model format, packaged runtime, explicit output paths, cache identity, and restart checks | Passed |

The shared NPY parsing and byte-count checks, current numerical kernels, and compiled V1 weight values are preserved. Tests cover implementation consistency and reject unsupported model packages. They do not establish independent DFT accuracy, material-grid projection accuracy, or size/k-grid convergence. No new DFT or production calculation was run.

The magnetic SCF generator currently covers six recovered target configurations. Original targets for the other 20 training states have not been recovered. Generated inputs use atomic initialization and inherited base numerical settings; these checks do not establish historical SCF-branch or label reproduction.

See [7_spin](../7_spin/README.md#small-checks-and-native-regression) for test commands and [6_RVR](../6_RVR/README.md) for the required model export format. Model, basis, runtime, or input changes require matching data and fresh work directories.
