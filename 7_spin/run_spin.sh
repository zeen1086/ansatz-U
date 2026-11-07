#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
: "${ANSATZ_BUILD_DIR:?Set an external build directory}"
if (($# != 1)); then
    echo 'Usage: run_spin.sh /absolute/path/run.nml' >&2
    exit 1
fi
make -C "$ROOT" BUILD="$ANSATZ_BUILD_DIR" "$ANSATZ_BUILD_DIR/ansatzU"
exec "$ANSATZ_BUILD_DIR/ansatzU" "$1"
