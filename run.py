#!/usr/bin/env python3
"""Run a verified Fortran build from the current calculation directory."""
from pathlib import Path
import argparse
import os

from build import verify
from parameters import ROOT

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mpi", type=int, help="number of MPI ranks")
    parser.add_argument("source", type=Path)
    parser.add_argument("arguments", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    source = args.source if args.source.is_absolute() else ROOT / args.source
    binary = verify(source)
    command = [str(binary)] + args.arguments
    if args.mpi is not None:
        if args.mpi < 1:
            parser.error("--mpi must be positive")
        if source.stem in ("make_wan_fun_mpi", "make_dis_mpi") and args.mpi < 2:
            parser.error("this master/worker program needs at least two MPI ranks")
        command = ["mpirun", "-np", str(args.mpi)] + command
    os.execvp(command[0], command)
