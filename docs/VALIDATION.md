# V1 Validation Summary

Verified on 2026-11-24. The evaluator uses `N0 = sigmoid(Q0)`, `N1 = sigmoid(Q0 + A(N0-N_ref))`, and `p = p0 + K(N1-N_ref)`. Original coupled V1 weights are normalized and compiled into `A=L K` without refitting.

| Check | Result |
|---|---|
| Repository offline suite | 98 passed; 1 MPI check skipped because the compiler and mpirun are unavailable |
| DFT input preparation | 26-state coverage, original sampling ranges, integer upper-hemisphere angles, seed repeatability, subset consistency, external pseudopotential paths, primitive-to-supercell generation, preserved parameters, isolated scratch, and invalid-input rejection passed |
| SCF-to-V1 label package | Nine offline tests passed for complex N/order, source consistency, PBE, GGA branch, grid/cell, scalar gauge and atomic publication; real saved-state N/XML parsing passed. Native exporter syntax/interface and build dry-run checks passed; linking and running against actual QE libraries was not performed |
| Single-response fit and spatial kernels | Twelve checks passed: label datasets/reference, complete joint-loss gradients, short synthetic fit, checkpoint/history recovery, direct bundle export, D4h fields and actual compiled native interpolation/window parity |
| Root layout and moved preparation commands | Offline suite passed after relocation; external-directory configure/build/run, argument forwarding, working-directory preservation, and stale-build rejection passed |
| Clean external builds through default `make` and `tools/build.py` | Passed; Python build fingerprints verified |
| Single-response algebra, full complex N, symmetry, occupation bounds, and parameter derivatives | Passed |
| Actual iteration-30 source conversion | 14 parameter arrays identical to the reference conversion |
| Native synthetic RVR/H/dH, L=1,2,3, independent real-grid integration | 120 component comparisons; maximum error 1.353e-16 eV |
| Cached versus direct native assembly | 120 component comparisons; maximum error 2.269e-16 eV |
| Batch versus single native assembly | 40 component comparisons; identical results |
| Actual V1 export → packaged runtime → synthetic RVR → native H/dH | Maximum error 9.155e-16 eV; all 15 compact/debug matrices identical |
| Model format, packaged runtime, explicit output paths, cache identity, and restart checks | Passed |
| Merged bundle loader and existing exported packages | Four-module exports and existing five-module layouts give identical synthetic p/N; missing loaders and mixed runtime modules are rejected; both layouts retain runtime identity checks |

The shared NPY parsing and byte-count checks, current numerical kernels, and compiled V1 weight values are preserved. Tests cover implementation consistency and reject unsupported model packages. They do not establish independent DFT accuracy, material-grid projection accuracy, or size/k-grid convergence. No new DFT or production calculation was run.

The magnetic SCF generator produces 26 configurations using the original sampling rules and a selected seed. Generated inputs use atomic initialization and inherited base numerical settings; these checks do not establish SCF convergence or material accuracy.

The label extractor reads complete complex N from native `occup.txt` and checks its magnitudes against the saved XML. It rebuilds dense-grid local fields with the original SCF magnetic GGA branch and no explicit constraint field. The offline adapter test uses synthetic output and does not establish a material-level extraction regression; that requires the matching QE 7.6 build.

The new fitter uses only the current single explicit response and starts from the label packages. Short synthetic optimization verifies that loss decreases and its bundle is accepted by the existing exporter, including unchanged predictions after export. Spatial kernels are learned within the configured finite A1g/B1g radial-polynomial span. No real-material fit or reproduction of the historical frozen kernels/accuracy was performed. Optimizer termination, gradient convergence and material accuracy remain separate results.

See [7_spin](../7_spin/README.md#small-checks-and-native-regression) for test commands and [6_RVR](../6_RVR/README.md) for the required model export format. Model, basis, runtime, or input changes require matching data and fresh work directories.
