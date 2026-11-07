#!/usr/bin/env python3
"""Apply common sources and one variant to a clean QE 7.6 source tree."""
from pathlib import Path
import argparse
import hashlib
import json
import shutil

HERE=Path(__file__).resolve().parent


def install(target: Path, variant: str) -> None:
    manifest=json.loads((HERE/"upstream.json").read_text())
    for name, expected in manifest["files"].items():
        path=target/name
        if (expected is None and path.exists()) or (expected is not None and (not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=expected)):
            raise ValueError(f"Expected an unmodified QE 7.6 source: {path}")
    for folder in ("common",variant):
        for source in (HERE/folder).rglob("*"):
            if source.is_file():
                output=target/source.relative_to(HERE/folder)
                output.parent.mkdir(parents=True,exist_ok=True)
                shutil.copy2(source,output)
    print(f"Applied QE 7.6 {variant} sources to {target}")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("variant",choices=("magnetic","nonmagnetic"))
    parser.add_argument("target",type=Path)
    args=parser.parse_args()
    install(args.target.resolve(),args.variant)
