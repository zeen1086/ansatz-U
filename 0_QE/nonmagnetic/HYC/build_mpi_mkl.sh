#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
: "${MKLROOT:?Load Intel oneAPI before building}"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 MKL_DYNAMIC=FALSE
./configure --enable-parallel --disable-openmp --with-scalapack=intel \
 MPIF90=mpiifort F90=ifort CC=icc \
 BLAS_LIBS="-L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl" \
 LAPACK_LIBS=" " \
 SCALAPACK_LIBS="-L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64"
# QE may autodetect MPI ScaLAPACK; the sequential BLAS kernel remains mandatory.
if grep -E '^(BLAS_LIBS|LAPACK_LIBS|SCALAPACK_LIBS).*mkl_(intel|gnu)_thread' make.inc; then
 echo 'Refusing threaded MKL' >&2; exit 1
fi
make -j"${BUILD_JOBS:-2}" pw
make -C PP -j"${BUILD_JOBS:-2}" all
ln -sfn ../PW/src/pw.x bin/pw_no_mag.x
