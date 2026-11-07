# QE 7.6 modifications

Source overlays for the official [QE 7.6 tag](https://gitlab.com/QEF/q-e/-/tree/qe-7.6),
based on the magnetic and nonmagnetic sources used on ISSP.

- `common/`: direction constraints, potential export and real-space wavefunction export.
- `magnetic/`: magnetic build entry point.
- `nonmagnetic/`: fixed scalar-potential reader, Hamiltonian update and build entry point.

Start with two clean QE 7.6 trees, including the pinned external dependencies
(`git clone --branch qe-7.6 --recursive https://gitlab.com/QEF/q-e.git` or the complete release).
`install.py` verifies the original files against `upstream.json` before copying common and
variant sources. It rejects modified or incompatible source trees.

```bash
python3 install.py magnetic /path/to/qe-7.6-magnetic
python3 install.py nonmagnetic /path/to/qe-7.6-nonmagnetic
cd /path/to/qe-7.6-magnetic
bash HYC/build_mpi_mkl.sh
# Repeat the build in the nonmagnetic tree.
```

Load Intel MPI/MKL first. The helper builds `pw.x`, `pp_exam.x`, `bands.x`,
`pw2wannier90.x`, `wfck2r.x` and the other PP tools, using sequential MKL.
The nonmagnetic tree also provides `bin/pw_no_mag.x`. Source `HYC/mpi_env.sh` before MPI runs.
Use the magnetic build for magnetic SCF/export and the nonmagnetic build for the fixed-scalar stages.

`pp_exam.x < print_V.in` writes `V_HYC_dp`: a three-int grid record followed by float64
`(V0,Bx,By,Bz)` in Ry, x fastest. The generated input enables `write_v_qe` for fitting.
The nonmagnetic executable requires a matching `V_HYC_dp` in its working directory;
it reads the scalar channel once and reapplies zero local magnetic field at each update.

The added `wfck2r.x` option `write_hyc=true` writes sequential complex128 records
`WFN/wfn_k*_b*_vec_dp`, with physical FFT-grid points (padding removed) and spinor components.
The original direct-access and Octave modes remain available with `write_hyc=false`.
Use the workflow's generated `wfck2r.in` for Wannier conversion.

For direction constraints, set `constrained_magnetization='atomic vector direction'`
and select species with `constrain_species` in `&SYSTEM`. Direction is constrained;
moment magnitude remains free. The original copyright notices and [GPL license](License) are retained.
