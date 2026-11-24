#!/usr/bin/env python3
"""Convert standard Wannier90 U/disentanglement and HR files to the HYC records."""
from pathlib import Path
import argparse
import sys
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tools.fortran_io import new_outputs, write_record
from tools.parameters import load


def read_u(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with path.open() as stream:
        next(stream)
        nk, columns, rows = map(int, next(stream).split())
        if min(nk, rows, columns) < 1:
            raise ValueError("invalid Wannier U matrix dimensions")
        result = np.empty((rows, columns, nk), dtype="<c16", order="F")
        kpoints = np.empty((nk, 3))
        for k in range(nk):
            line = next(stream)
            while not line.strip():
                line = next(stream)
            kpoints[k] = [float(x) for x in line.split()]
            values = np.array([[float(x) for x in next(stream).split()] for _ in range(rows*columns)])
            if values.shape != (rows*columns, 2):
                raise ValueError("invalid complex Wannier matrix row")
            result[:,:,k] = (values[:,0]+1j*values[:,1]).reshape((rows, columns), order="F")
        if any(line.strip() for line in stream):
            raise ValueError("unexpected trailing Wannier U data")
    if not np.isfinite(result).all() or not np.isfinite(kpoints).all():
        raise ValueError("non-finite Wannier U matrix")
    return result, kpoints


def read_hr(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with path.open() as stream:
        next(stream)
        nw, nr = int(next(stream)), int(next(stream))
        if nw < 1 or nr < 1:
            raise ValueError("invalid HR dimensions")
        degeneracy = []
        while len(degeneracy) < nr:
            degeneracy += [int(x) for x in next(stream).split()]
        if len(degeneracy) != nr or min(degeneracy) < 1:
            raise ValueError("invalid HR degeneracy")
        vectors, matrices = [], []
        for _ in range(nr):
            matrix = np.zeros((nw,nw), dtype="<c16", order="F")
            seen, vector = set(), None
            for _ in range(nw*nw):
                row = next(stream).split()
                if len(row) != 7:
                    raise ValueError("invalid HR row")
                v = tuple(map(int,row[:3]))
                i,j = int(row[3])-1, int(row[4])-1
                if not 0 <= i < nw or not 0 <= j < nw or (i,j) in seen:
                    raise ValueError("invalid or duplicate HR matrix index")
                if vector is not None and vector != v:
                    raise ValueError("mixed translations in HR block")
                vector = v
                seen.add((i,j))
                matrix[i,j] = complex(float(row[5]),float(row[6]))
            vectors.append(vector)
            matrices.append(matrix)
        if len(set(vectors)) != nr or any(line.strip() for line in stream):
            raise ValueError("duplicate translations or trailing HR data")
    matrices = np.stack(matrices, axis=2)
    if not np.isfinite(matrices).all():
        raise ValueError("non-finite HR elements")
    return np.array(vectors), np.array(degeneracy), matrices


def prepare(source: Path, output: Path, p: dict) -> None:
    seed = source / p["prefix"]
    u, kpoints = read_u(seed.with_name(seed.name+"_u.mat"))
    dis, kd = read_u(seed.with_name(seed.name+"_u_dis.mat"))
    nk, nw, nb = p["nk1"]*p["nk2"]*p["nk3"], p["nwan"], p["num_bands"]
    if u.shape != (nw,nw,nk) or dis.shape != (nb,nw,nk) or not np.allclose(kpoints,kd,atol=1e-10,rtol=0):
        raise ValueError("Wannier matrices disagree with parameter or k-point ordering")
    # This workflow uses the full excluded-band complement as the outer window.
    win = (source / (p["prefix"]+".win")).read_text().lower()
    import re
    active = "\n".join(line.split("!")[0].split("#")[0] for line in win.splitlines())
    if not re.search(r"use_ws_distance\s*=\s*(?:\.false\.|false|f)(?:\s|$)",active) or re.search(r"\bdis_win_(min|max)\s*=",active) or re.search(r"use_ws_distance\s*=\s*(?:\.true\.|true|t)\b",active):
        raise ValueError("require no outer dis_win bounds and use_ws_distance=false")
    vectors, deg, hr = read_hr(seed.with_name(seed.name+"_hr.dat"))
    if hr.shape[:2] != (nw,nw):
        raise ValueError("HR nwan disagrees with parameter")
    names = ["u_mat_HYC", "u_dis_mat_HYC", "window_TF_HYC", "k_list", "ir_vec_HYC", "H_deg_HYC"]
    names += [f"hr_HYC_{n+1}" for n in range(len(deg))]
    new_outputs(output,names)
    for name, matrix in (("u_mat_HYC",u),("u_dis_mat_HYC",dis)):
        with (output/name).open("wb") as stream:
            write_record(stream,np.array(matrix.shape,dtype="<i4"))
            write_record(stream,matrix)
    with (output/"window_TF_HYC").open("wb") as stream:
        write_record(stream,np.ones((nb,nk),dtype="<i4",order="F"))
    np.savetxt(output/"k_list",kpoints,fmt="%12.8f%12.8f%12.8f")
    np.savetxt(output/"ir_vec_HYC",vectors,fmt="%8d%8d%8d")
    np.savetxt(output/"H_deg_HYC",np.column_stack((np.arange(1,len(deg)+1),deg)),fmt="%8d%8d")
    for n in range(len(deg)):
        with (output/f"hr_HYC_{n+1}").open("wb") as stream:
            write_record(stream,hr[:,:,n])


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path,help="directory for new converted records; existing outputs are rejected")
    args=parser.parse_args()
    prepare(Path(load()["nscf_dir"]),args.output,load())
