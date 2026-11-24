"""Read QE/Wannier preparation parameters and generate build-time Fortran values."""
from pathlib import Path
import math
import os
import shlex
import sys
try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib  # Python 3.10 hosts may use an already installed tomli.

ROOT = Path(__file__).resolve().parents[1]
PATH_KEYS = {"nscf_dir", "pseudo_dir"}


def load(path: Path | None = None) -> dict:
    if path is None:
        path = ROOT / os.environ.get("ANSATZ_PARAMETER", "parameter")
    path = Path(path).resolve()
    sections = tomllib.loads(path.read_text())
    values = {}
    for section in sections.values():
        for key, value in section.items():
            if key in values:
                raise ValueError(f"duplicate global parameter: {key}")
            values[key] = value
    positive = ("nwan", "nspin", "nr1", "nr2", "nr3", "nk1", "nk2", "nk3")
    for key in positive:
        if type(values[key]) is not int or values[key] < 1:
            raise ValueError(f"{key} must be a positive integer")
    for key, value in values.items():
        for number in value if isinstance(value, list) else [value]:
            if isinstance(number, float) and not math.isfinite(number):
                raise ValueError(f"{key} must be finite")
    if values["nspin"] not in (1, 2):
        raise ValueError("require nspin=1/2")
    if type(values["band_segments"]) is not int or not 1 <= values["band_segments"] <= 999:
        raise ValueError("band_segments must be in 1..999")
    if not 1 <= values["wan_min"] <= values["wan_max"] <= values["nwan"]:
        raise ValueError("require 1 <= wan_min <= wan_max <= nwan")
    if (not 1 <= values["first_band"] <= values["last_band"] <= values["nbnd"]
            or values["num_bands"] != values["last_band"]-values["first_band"]+1
            or values["num_bands"] < values["nwan"]):
        raise ValueError("Wannier band window, num_bands, nwan and nbnd are inconsistent")
    excluded = set()
    for item in values["exclude_bands"].replace(",", " ").split():
        bounds = [int(x) for x in item.split("-")]
        if len(bounds) > 2 or min(bounds) < 1 or max(bounds) > values["nbnd"] or bounds[-1] < bounds[0]:
            raise ValueError("invalid exclude_bands")
        excluded.update(range(bounds[0], bounds[-1]+1))
    expected = set(range(1, values["first_band"])) | set(range(values["last_band"]+1, values["nbnd"]+1))
    if excluded != expected:
        raise ValueError("exclude_bands must leave exactly first_band..last_band")
    for axis, k in zip("xyz", ("nk1", "nk2", "nk3")):
        extent = values[f"convert_size_{axis}"]
        if type(extent) is not int or extent < 0 or 2*extent+1 > values[k]:
            raise ValueError(f"invalid Wannier support on {axis}")
    if not isinstance(values["prefix"], str) or not values["prefix"].replace("_", "").isalnum():
        raise ValueError("prefix must contain letters, digits or underscores")
    if values["export_grid"] not in ("smooth", "dense") or type(values["write_v_qe"]) is not bool:
        raise ValueError("export_grid must be smooth/dense and write_v_qe must be boolean")
    for key in ("a0", "celldm3", "ecutwfc", "conv_thr", "ni_magnetization"):
        if isinstance(values[key], bool) or not isinstance(values[key], (int, float)):
            raise ValueError(f"{key} must be numeric")
    if values["a0"] <= 0 or values["celldm3"] <= 0 or values["conv_thr"] <= 0:
        raise ValueError("lattice constants and conv_thr must be positive")
    if values["ni_magnetization"] < 0:
        raise ValueError("ni_magnetization must be nonnegative")
    for key in PATH_KEYS:
        values[key] = str((path.parent / values[key]).resolve())
    return values


def fortran_module(values: dict) -> str:
    lines = ["! Generated from parameter; do not edit.", "module parameter_mod", "  implicit none"]
    for key, value in values.items():
        name = "p_" + key
        if isinstance(value, list):
            integers = all(type(x) is int for x in value)
            kind = "integer" if integers else "real(8)"
            lines.append(f"  {kind}, parameter :: {name}({len(value)}) = [ &")
            lines.extend((str(x) if integers else f"{float(x):.17e}_8") + (", &" if i < len(value)-1 else " ]")
                         for i, x in enumerate(value))
        elif isinstance(value, bool):
            lines.append(f"  logical, parameter :: {name} = .{str(value).lower()}.")
        elif isinstance(value, int):
            lines.append(f"  integer, parameter :: {name} = {value}")
        elif isinstance(value, float):
            lines.append(f"  real(8), parameter :: {name} = {value:.17e}_8")
        else:
            text = value + ("/" if key in PATH_KEYS else "")
            chunks = [text[i:i+64].replace("'", "''") for i in range(0, len(text), 64)] or [""]
            lines.append(f"  character(len=*), parameter :: {name} = &")
            lines.extend("    '" + chunk + "'" + (" // &" if i < len(chunks)-1 else "")
                         for i, chunk in enumerate(chunks))
    lines += ["contains", "  subroutine read_row_index(row)", "    integer, intent(out) :: row",
              "    integer :: ios", "    character(len=32) :: argument",
              "    call get_command_argument(1, argument)",
              "    read(argument, *, iostat=ios) row",
              "    if (ios /= 0) error stop 'First argument must be a Wannier row index'",
              "    if (row < 1 .or. row > p_nwan*2/p_nspin) error stop 'row index out of range'",
              "  end subroutine read_row_index", "  function row_text(row) result(value)",
              "    integer, intent(in) :: row", "    character(len=20) :: value",
              "    write(value, '(I0)') row", "  end function row_text"]
    return "\n".join(lines + ["end module parameter_mod", ""])


if __name__ == "__main__":
    values = load()
    if sys.argv[1:] == ["--shell"]:
        for key, value in values.items():
            if not isinstance(value, list):
                print(f"{key.upper()}={shlex.quote(str(value))}")
    else:
        for key in sys.argv[1:]:
            print(values[key])
