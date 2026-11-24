# ansatz+U

Run the commands below from the repository root. Replace `/absolute/...` with
your paths. Keep pseudopotentials, calculation data, models and outputs outside
the repository.

## 1. Set up and build

Use an existing Python 3.11+ environment with NumPy/SciPy, gfortran and
BLAS/LAPACK (Accelerate on macOS). For DFT, build QE 7.6 using [0_QE](0_QE/README.md).

```sh
export ANSATZ_PYTHON=/absolute/existing/python
make BUILD=/absolute/external/build
```

This builds `local/rvr` and `local/cal_Hsk` under the specified build directory.

## 2. Prepare DFT inputs

Set the structure, pseudopotential paths and numerical settings in the input
templates and `parameter`, then generate the inputs:

```sh
"$ANSATZ_PYTHON" tools/configure.py --output-dir /absolute/external/case
```

For preparation executables, build with `tools/build.py SOURCE` and run with
`tools/build.py run [--mpi N] SOURCE [arguments]`, using the same Python interpreter.

## 3. Run the workflow

| Step | Commands and required inputs |
|---|---|
| Generate magnetic configurations → run SCF → extract labels → fit → export | Follow [4_potential](4_potential/README.md), starting from the base `scf.in` |
| Project the model into the Wannier basis | Follow [6_RVR](6_RVR/README.md); supply the exported model, matching Wannier functions and DFT-only background |
| Calculate Hamiltonians, energies and torques | Follow [7_spin](7_spin/README.md); supply the model, `B_general`, background and spin configurations |

## Check the installation

```sh
"$ANSATZ_PYTHON" -B tests/check.py
```

Test coverage and numerical limits: [validation](docs/VALIDATION.md).
