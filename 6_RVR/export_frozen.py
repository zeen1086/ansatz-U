"""One-time, lossless storage conversion; never used by the Fortran numerical run."""
from pathlib import Path
import argparse
import hashlib
import json
import struct
import numpy as np

FROZEN_SHA = '6fe8cb35d6271638fea0fe9227456a4f1b04910fb8cf181fe30cc212420c2ba0'


def export(source, output, expected_sha):
    source, output = Path(source), Path(output)
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    if digest != expected_sha:
        raise ValueError(f'Frozen input SHA mismatch: {digest}')
    with np.load(source, allow_pickle=False) as z:
        arrays = [np.asarray(z['rr'].T, dtype='<i4'),
                  np.asarray(z['hc_R'].transpose(1, 2, 0), dtype='<c16'),
                  np.asarray(z['p_R'].transpose(1, 2, 0), dtype='<c16'),
                  np.asarray(z['KB'], dtype='<c16'), np.asarray(z['N0'], dtype='<c16')]
        assert np.array_equal(arrays[0].T, z['rr'])
    assert [a.shape for a in arrays] == [(3,729),(30,30,729),(5,30,729),(810,810),(10,10)]
    assert all(np.isfinite(a).all() for a in arrays)
    assert np.max(abs(arrays[3]-arrays[3].conj().T)) < 1e-10
    assert np.max(abs(arrays[4]-arrays[4].conj().T)) < 1e-12
    with output.open('xb') as f:
        f.write(b'LANIO2_AU_V1    ')
        f.write(struct.pack('<8i',0x01020304,729,30,5,10,810,1,0))
        for a in arrays:
            f.write(a.tobytes(order='F'))
    with output.open('rb') as f:
        assert f.read(16) == b'LANIO2_AU_V1    '
        assert struct.unpack('<8i', f.read(32)) == (0x01020304,729,30,5,10,810,1,0)
        for a in arrays:
            b=np.frombuffer(f.read(a.nbytes),dtype=a.dtype).reshape(a.shape,order='F')
            assert np.array_equal(a,b), 'Lossless readback failed'
        assert not f.read(1)
    meta=dict(source=str(source),source_sha256=digest,output=str(output),
              output_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),bytes=output.stat().st_size,
              endian='little',integer='int32',real='IEEE binary64',complex='pair of binary64 (real,imag)',
              array_order='Fortran column-major',arrays=[dict(shape=a.shape,dtype=str(a.dtype)) for a in arrays],
              energy_unit='eV',background='hc_R=H0+fixed scalar Hv; contains NO HU',
              source_components='H0 and Hv are combined in this frozen NPZ; no individual decomposition is invented',
              projector_order='rows=[up d1..d5,down d1..d5]; Hamiltonian=(cell,orbital,up/down)',
              translations='rr=WS [-4,4]^3; KB=R3 [-1,1]^3 x,y,z then orbital',
              lossless_readback=True)
    output.with_suffix(output.suffix+'.json').write_text(json.dumps(meta,indent=2)+'\n')
    return meta

if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('source');p.add_argument('output')
    p.add_argument('--expected-sha256',default=FROZEN_SHA)
    a=p.parse_args()
    print(json.dumps(export(a.source,a.output,a.expected_sha256),indent=2))
