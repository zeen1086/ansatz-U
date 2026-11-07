#!/usr/bin/env bash
# Run inside an existing allocation; no jobs are submitted.
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT=$(cd "$SCRIPT_DIR/../.." && pwd)
PYTHON=${PYTHON:-python3}
"$PYTHON" "$ROOT/build.py" --fc "${MPIFC:-${FC:-mpifort}}" "$ROOT/5_wannier/make_wan/make_wan_fun_mpi.f90"
WAN_MIN=$("$PYTHON" parameters.py wan_min)
WAN_MAX=$("$PYTHON" parameters.py wan_max)
MPI_RANKS=${MPI_RANKS:-2}
if ((MPI_RANKS < 2)); then
    echo 'The master/worker Wannier generator needs at least two MPI ranks.' >&2
    exit 1
fi
for ((band=WAN_MIN; band<=WAN_MAX; band++)); do
    "$PYTHON" "$ROOT/run.py" --mpi "$MPI_RANKS" 5_wannier/make_wan/make_wan_fun_mpi.f90 "$band"
done
