"""Native little-endian Fortran sequential records (four-byte markers)."""
from pathlib import Path
import struct
import numpy as np


def read_record(stream) -> bytes:
    chunks = []
    previous = False
    while True:
        marker = stream.read(4)
        if len(marker) != 4:
            raise ValueError("missing Fortran record marker")
        lead = struct.unpack("<i", marker)[0]
        data = stream.read(abs(lead))
        tail = stream.read(4)
        if len(data) != abs(lead) or len(tail) != 4:
            raise ValueError("truncated Fortran record")
        end = struct.unpack("<i", tail)[0]
        if abs(end) != abs(lead) or (end < 0) != previous:
            raise ValueError("inconsistent Fortran record markers")
        chunks.append(data)
        if lead >= 0:
            return b"".join(chunks)
        previous = True


def read_array(stream, dtype, shape: tuple) -> np.ndarray:
    data = read_record(stream)
    if len(data) != int(np.prod(shape))*np.dtype(dtype).itemsize:
        raise ValueError(f"Fortran payload does not match {shape}")
    result = np.frombuffer(data, dtype=dtype).reshape(shape, order="F")
    if not np.isfinite(result).all():
        raise ValueError("non-finite values in Fortran record")
    return result


def write_record(stream, *arrays) -> None:
    data = b"".join(a.tobytes(order="F") if isinstance(a, np.ndarray) else a for a in arrays)
    maximum = 2147483639
    for offset in range(0, max(1, len(data)), maximum):
        block = data[offset:offset+maximum]
        first, last = offset == 0, offset+len(block) == len(data)
        stream.write(struct.pack("<i", len(block) if last else -len(block)))
        stream.write(block)
        stream.write(struct.pack("<i", len(block) if first else -len(block)))


def new_outputs(path: Path, names: list[str]) -> None:
    for name in names:
        if (path/name).exists():
            raise FileExistsError(f"refusing to replace existing output: {path/name}")
    path.mkdir(parents=True, exist_ok=True)
