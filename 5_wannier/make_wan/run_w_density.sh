#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT=$(cd "$SCRIPT_DIR/../.." && pwd)
PYTHON=${PYTHON:-python3}
"$PYTHON" "$ROOT/build.py" "$ROOT/5_wannier/make_wan/make_w_density.f90"
"$PYTHON" "$ROOT/run.py" 5_wannier/make_wan/make_w_density.f90 | tee make_w_density.log
