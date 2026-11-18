# Validation Summary

The local, tensor, and `joint` RVR tables below are historical validation records for the [pre-cleanup code](https://github.com/zeen1086/ansatz-U/tree/3e1307ac721e91a66eb3851dae2df6bc4c9c99e1). These legacy entry points were retired on 2026-11-16; “current entry” and `make all` in the tables refer to that version. See [6_RVR](6_RVR/README.md) for the current V1 cleanup regression and [production](7_spin/README.md) for production readout tests.

| Check | Result |
|---|---|
| Clean build of LaNiO₂ sources, diagonal/full-component fitting and restart, original-stage combination | Passed |
| Common ansatz offline tests, Makefile and Python build entry points | 16 passed, 1 skipped because no MPI compiler was available; both build entry points passed |
| 12 regression cases for nspin column bounds and RVR lengths/markers | Passed |
| H for frozen 2×2 inputs against the existing Python reference | Maximum difference 3.20e-14 eV |
| Eigenvalues, EF, and B/U torques against an independent reference | Maximum differences 1.07e-13 eV, 1.78e-14 eV, 5.01e-15 eV/rad |
| Φ against the existing Python reference | Maximum difference 5.26e-11 meV/rad² |
| Three-dimensional spin directions, N eigenvalue preservation, HB/HU finite differences, 27-translation tensor kernel | Passed |
| Nondefault Ueff/sigma/electron count/3D axes/readout sites, two tensor amplitudes, and Ry/eV | Independent-reference H difference ≤7.97e-15 eV, torque difference ≤1.18e-14 eV/rad |
| Central differences of B/U derivatives and fixed-electron-number free-energy derivatives | Matrix difference ≤1.50e-11 eV/rad, total torque difference ≤3.08e-8 eV/rad |
| Original RVR combination → two H_sk sets → background combination | nspin=1/2 passed; unsupported extra channels rejected |
| Current entry against the existing Python local implementation (U=1, sigma=.05) | H difference 1.08e-14 eV, Φ difference 8.38e-10 meV/rad² |

Code thresholds: H/eigenvalues/EF 1e-10 eV, torque 1e-9 eV/rad, electron count 1e-8, Φ 1e-6 meV/rad²; finite-difference matrices 2e-8 eV/rad and free-energy torque 2e-7 eV/rad. The additional comparisons above used existing small synthetic models.

Full-size grids, reprojection, angular-step and spatial-range convergence, and independent material validation remain incomplete. B/U are derivative contributions of the same H; see the README for the spatial and orbital scope of the +U entry. Test and run outputs are stored outside the repository.

## Native Four-Spin, Three-Potential SVD Interface

The following implementation checks cover fitting/evaluation with the new `V_func_tmp` (formerly `model.joint`) and are recorded separately from the local `.au` path above.

| Check | Result |
|---|---|
| Clean external `make all` build; shell unpacking → Fortran input → SVD fitting → model readback → three-potential/derivative evaluation | Passed; no Python invoked in the production chain |
| 32-state synthetic example with 2 hard anchors, 1 zero-weight anomalous holdout, mixed Ry/eV and C/Fortran array order | Passed; scalar8, B40, N_vector12, N_cross6, and N_chi4 all reached their respective full ranks |
| Scalar/B coefficients against independent NumPy constrained SVD | Maximum differences 7.13e-15 / 2.56e-15 eV |
| Scalar/B, full N, and HU against original Python feature evaluation | Maximum differences 8.68e-18 eV, 4.45e-16, 4.45e-16 eV |
| Central differences of analytical scalar/B and N derivatives | Native checks ≤1.76e-12 eV/rad, 2.05e-11 /rad; independent Python checks ≤5.34e-13, 6.41e-12 |
| 3D directions, global SO(3)/SU(2), time reversal, 64 crystal operations, nonrigid N eigenvalue changes, and neighbor dN | Passed |
| Reading existing real 93° DFT labels on an 80×80×36 grid: full complex N, units, atom order, and VU relation | Passed; the full real dataset was not trained |
| 16 boundary/parameter cases | Passed; included rejection of single precision/truncation/order/non-Hermitian/HU mismatches, anchor conflicts, out-of-range occupations, overwrite protection, U=0/1.7, perturbed Ni selection, and normalization switch |
| Explicit scalar-mean removal switch and model metadata | Passed; disabled by default, with the removed amount recorded in the report |
| Existing offline regression | 16 passed, 1 skipped because no MPI compiler was available |

Thresholds for the new interface are 1e-10 in the respective units for coefficient/potential/N readback and independent evaluation, and 1e-8 for central differences of angular derivatives. The default SVD relative cutoff is 1e-10; the effective cutoff after applying the floating-point floor is reported. The 5 meV physical target was not used in place of code-consistency thresholds.
The initial SVD validation covered only DFT-label reading, three-potential fitting, and evaluation/derivative implementation. The new-feature RVR/Hamiltonian/torque connection, which was not covered then, is described below. Full real 26-state fitting, equivalence to the old projected-operator loss, and nonrigid N extension to large supercells remain unverified.

2026-11-08 filename regression: a clean out-of-tree `make all` build passed. An independently built baseline from pre-rename commit `2fe6859` reused the 32-state small-grid synthetic inputs. The new `V_func_tmp` and old `model.joint`, and the four renamed evaluation files and their predecessors, were byte-for-byte identical, with zero numerical difference. `fit_report.txt` was identical apart from the explicit input directory. Default/explicit paths, output prefixes, format-header/version/truncation rejection, and overwrite protection for models/run directories/evaluation outputs passed. The new shell script does not generate `rows.manifest`; historical files are unchanged. No real retraining, projection, or HPC calculation was performed.

## Three-Potential RVR → Full H → Three Torque Contributions/Φ

The following uses a synthetic periodic four-Ni cell outside the source tree to validate the native production chain. No real large-grid projection or HPC job was run.

| Check | Result |
|---|---|
| Clean `make all`, including new cal_RVR_split/combine and ansatzU; gfortran double precision/BLAS/LAPACK | Passed; the final build also enabled `-fcheck=all` |
| 120 spatial orbitals, 256 grid points, complex Wannier functions, 3 translations, 3 unequally weighted 3D K points; complete 120-row Fortran projection/combination/readback | Passed; maximum difference of 49 spatial kernels against independent grid contraction 4.45e-16 eV |
| `hbg+[Hv-Hv_ref]+HB+HU`, independent Bloch background and original four-Ni projectors | Maximum full-H difference 7.46e-16 eV; each component checked separately, with no scalar or HU double counting |
| Eigenvalues, EF, FD occupations, and fixed electron count | Maximum differences 7.33e-15 eV, 1.25e-16 eV, 5.89e-15; target electron-count error 0 |
| v/B/U angular derivatives against independent central differences | Maximum differences 2.63e-11 / 6.57e-12 / 1.49e-11 eV/rad |
| Three torque contributions and total torque from the same full H/occupations | Maximum differences 1.05e-16 / 8.33e-17 eV/rad; nonzero scalar and U contributions |
| Scalar/B/U and total Φ | Maximum differences 7.26e-12 / 6.81e-12 meV/rad² |
| Fixed-electron-number free-energy angular differences against total torque | Maximum difference 8.34e-10 eV/rad |
| Nonrigid full complex N and neighboring-Ni derivatives; U=0/1.7, sigma=.083, 19.3 electrons/Ni, another perturbed Ni | Passed; maximum imaginary part of N 1.09e-3, neighbor dN 1.46e-4/rad; HU/dHU exactly zero at U=0 |
| Global SU(2) rotation of full H and each component; invariance of three torque contributions and EF | Maximum matrix difference 4.48e-16 eV; torque 1.39e-16 eV/rad; EF 5.21e-18 eV |
| Original primitive30-to-four-cell-120 translation/phase/support mapping | Native projection and independent expansion byte-for-byte identical; Bloch-folding difference 7.81e-18; truncation of existing support rejected |
| Move aside original V_func_tmp, frozen model, reference fields, and Wannier directories; read only portable B_general | 53 first-state outputs byte-for-byte identical, without reprojection or grid reads |
| Model/format/background/overwrite protection, mixed-row rejection, quadrature weights and cal_TF support actually used | Passed; independent scaling difference for weight 1.7 ≤8.33e-16 eV |
| Original local U=0, U=.37 3D texture, U=1 Table 3, nondefault electron count/FD, and original three-potential grid evaluation | U=0/Table 3/nondefault numerics and original three-potential evaluation byte-for-byte identical; final-build 3D U=.37 maximum low-order difference 4.85e-19, below the original 1e-10 threshold |
| Original repository offline tests | 16 passed; 1 skipped because no MPI compiler was available |

The prescribed code thresholds were retained: matrices/eigenvalues/EF 1e-10 eV, torque 1e-9 eV/rad, electron count 1e-8, Φ 1e-6 meV/rad²; derivative finite differences 2e-8 eV/rad and fixed-number free-energy torque 2e-7 eV/rad. They were not relaxed to physical RMS thresholds.

These conclusions cover the synthetic 2×2 implementation and legacy-entry regression above. Binding the reference scalar field/common gauge of the real frozen model, full DFT fitting/projection and real Table 3, spatial/angular-step convergence, and large-cell N extension remain incomplete. Physical validation of the old local model does not carry over.
