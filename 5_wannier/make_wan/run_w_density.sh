#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT=$(cd "$SCRIPT_DIR/../.." && pwd)
PYTHON=${PYTHON:-python3}
cd "$ROOT"
"$PYTHON" build.py 5_wannier/make_wan/make_w_density.f90
cd "$SCRIPT_DIR"
"$PYTHON" "$ROOT/run.py" 5_wannier/make_wan/make_w_density.f90 | tee make_w_density.log
