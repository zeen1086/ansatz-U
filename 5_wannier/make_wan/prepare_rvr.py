#!/usr/bin/env python3
"""Map sequential Wannier records onto R_vec and write raw RVR input streams."""
from pathlib import Path
import argparse
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from fortran_io import new_outputs, read_array
from parameters import load


def prepare(source: Path, output: Path, p: dict) -> None:
    with (Path(p["nscf_dir"])/"V_HYC_dp").open("rb") as stream:
        grid = read_array(stream,"<i4",(3,))
    if np.any(grid < 1):
        raise ValueError("invalid Wannier grid")
    with (Path(p["func_dir"])/"R_vec").open("rb") as stream:
        ni=int(read_array(stream,"<i4",(1,))[0])
        vectors=read_array(stream,"<i4",(ni,3))
    if np.any(vectors < 0) or np.any(vectors >= grid):
        raise ValueError("R_vec indices outside the Wannier grid")
    nd=int(np.prod([2*p["convert_size_"+axis]+1 for axis in "xyz"]))
    new_outputs(output,["wan_data"]+[f"cal_TF_{row}" for row in range(1,p["nwan"]+1)])
    (output/"wan_data").mkdir()
    indices=vectors[:,0]+grid[0]*(vectors[:,1]+grid[1]*vectors[:,2])
    for row in range(1,p["nwan"]+1):
        with (source/f"wan_{row}").open("rb") as stream:
            values=read_array(stream,"<c16",(nd,int(np.prod(grid)),p["nspin"]))
            if stream.read(1):
                raise ValueError("unexpected trailing Wannier records")
        if float(np.sum(np.abs(values)**2)) <= 0:
            raise ValueError(f"zero Wannier norm in row {row}")
        (output/"wan_data"/f"wan_{row}").write_bytes(values[:,indices,:].tobytes(order="F"))
        # The original converter keeps every direction; support truncation is separate.
        np.savetxt(output/f"cal_TF_{row}",np.ones((1,nd),dtype=int),fmt="%8d")


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path,help="directory for new wan_data and cal_TF; existing outputs are rejected")
    args=parser.parse_args()
    p=load()
    prepare(Path(p["nscf_dir"])/"make_wan",args.output,p)
