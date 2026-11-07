#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
PYTHON=${PYTHON:-python3}
cd "$ROOT"
"$PYTHON" build.py 7_spin/cal_Hsk.f90 7_spin/combine.f90 7_spin/cal_torque_hf.f90
THREADS=$("$PYTHON" parameters.py threads)
export OMP_NUM_THREADS="$THREADS" MKL_NUM_THREADS=1 OMP_DYNAMIC=FALSE MKL_DYNAMIC=FALSE
cd "$SCRIPT_DIR"
mkdir -p data tot_HR
"$PYTHON" "$ROOT/run.py" 7_spin/cal_Hsk.f90
"$PYTHON" "$ROOT/run.py" 7_spin/combine.f90
"$PYTHON" "$ROOT/run.py" 7_spin/cal_torque_hf.f90
