#!/usr/bin/env python3
"""Compile selected Fortran programs, or verify their saved build fingerprints."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import shlex
import shutil
import subprocess
import sys

from parameters import ROOT, fortran_module, load

DEFAULT_LIBS = "-framework Accelerate" if sys.platform == "darwin" else "-llapack -lblas"

FIT_MODULES = ["wannier_weight_mod", "cri3_block_mask_mod", "cri3_c3_basis_mod"]


def build_root() -> Path:
    return (ROOT / os.environ.get("ANSATZ_BUILD_DIR", ".build")).resolve()


def locations(source: Path) -> tuple[Path, Path]:
    relative = source.resolve().relative_to(ROOT)
    work = build_root() / relative.with_suffix("")
    binary = (source.with_suffix(".x") if build_root() == ROOT / ".build"
              else build_root() / "bin" / relative.with_suffix(".x"))
    return work, binary


def dependencies(source: Path) -> list[Path]:
    if source.name == "fit_new.f90":
        return [ROOT / "4_potential/fit/fit_mask.f90", source]
    if source.name == "run_ansatzU.f90":
        return [ROOT / name for name in ("common/numerics.f90", "6_RVR/frozen_io.f90",
                "6_RVR/tensor_rvr.f90", "7_spin/local_ansatz.f90")] + [source]
    result = []
    if source.parent.name == "fit" and source.stem not in FIT_MODULES:
        result += [source.parent / (name + ".f90") for name in FIT_MODULES]
    if source.stem in ("cal_RVR_zgemm", "test_rvr_kernel"):
        result.append(ROOT / "6_RVR/speedup/rvr_kernel_mod.f90")
    if source.stem in ("cal_Hsk", "cal_Hsk_trace", "cal_torque_hf", "test_spin_channels"):
        result.append(ROOT / "7_spin/spin_channels_mod.f90")
    return result + [source]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def state(source: Path, fc: str, flags: list[str], libs: list[str], values: dict | None = None) -> dict:
    compiler = shlex.split(fc)
    resolved = shutil.which(compiler[0])
    if resolved is None:
        raise FileNotFoundError(f"Fortran compiler not found: {compiler[0]}")
    compiler[0] = str(Path(resolved).resolve())
    return {"parameter": hashlib.sha256(fortran_module(values if values is not None else load()).encode()).hexdigest(),
            "compiler": compiler, "compiler_mtime_ns": Path(compiler[0]).stat().st_mtime_ns,
            "flags": flags, "libs": libs,
            "sources": {str(p.relative_to(ROOT)): digest(p)
                        for p in dependencies(source) + [ROOT / "build.py", ROOT / "parameters.py"]}}


def verify(source: Path) -> Path:
    source = source.resolve()
    work, binary = locations(source)
    try:
        saved = json.loads((work / "build.json").read_text())
        expected = state(source, shlex.join(saved["compiler"]), saved["flags"], saved["libs"])
        if expected != saved["state"] or digest(binary) != saved["binary_sha256"]:
            raise ValueError("changed parameters, sources, compiler or binary")
    except (OSError, ValueError, KeyError) as exc:
        raise RuntimeError(f"Rebuild {source.relative_to(ROOT)} before running: {exc}") from exc
    return binary


def build(sources: list[Path], fc: str, flags: list[str], libs: list[str], check: bool = False) -> None:
    values = load()
    for source in sources:
        source = source.resolve()
        relative = source.relative_to(ROOT)
        if relative.parts[0] == "0_QE":
            raise ValueError("build these source replacements inside their matching QE tree")
        work, binary = locations(source)
        current = state(source, fc, flags, libs, values)
        if not check:
            try:
                saved = json.loads((work / "build.json").read_text())
                if current == saved["state"] and digest(binary) == saved["binary_sha256"]:
                    print(f"Up to date: {relative}")
                    continue
            except (OSError, ValueError, KeyError):
                pass
        work.mkdir(parents=True, exist_ok=True)
        module = work / "parameter_mod.f90"
        module.write_text(fortran_module(values))
        cmd = current["compiler"] + flags
        if check:
            cmd += ["-fsyntax-only"]
        cmd += [str(module)] + [str(p) for p in dependencies(source)]
        if not check:
            binary.parent.mkdir(parents=True, exist_ok=True)
            cmd += ["-o", str(binary) + ".new"] + libs
        print(shlex.join(cmd), flush=True)
        subprocess.run(cmd, cwd=work, check=True)
        if not check:
            Path(str(binary) + ".new").replace(binary)
            saved = {"state": current, "compiler": current["compiler"], "flags": flags,
                     "libs": libs, "binary_sha256": digest(binary)}
            (work / "build.json").write_text(json.dumps(saved, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sources", nargs="+", type=Path)
    parser.add_argument("--fc", default=os.environ.get("FC", "gfortran"))
    parser.add_argument("--flags", default=os.environ.get("FFLAGS", "-O2 -fopenmp -ffree-line-length-none -fallow-argument-mismatch"))
    parser.add_argument("--libs", default=os.environ.get("LIBS", DEFAULT_LIBS))
    parser.add_argument("--check", action="store_true", help="GNU syntax check only")
    parser.add_argument("--verify", action="store_true", help="reject stale or missing builds without compiling")
    args = parser.parse_args()
    if args.verify:
        for source in args.sources:
            print(f"Verified: {verify(source)}")
    else:
        build(args.sources, args.fc, shlex.split(args.flags), shlex.split(args.libs), args.check)
