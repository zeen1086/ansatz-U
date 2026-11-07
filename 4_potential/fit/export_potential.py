#!/usr/bin/env python3
"""Export a fitted scalar spin kernel as the diagonal tensor consumed by RVR."""
from pathlib import Path
import argparse
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from fortran_io import read_array, read_record, write_record
from parameters import load


def export(source: Path, output: Path, n_func: int) -> None:
    with source.open('rb') as stream:
        if read_record(stream).rstrip(b' ') != b'CRI3_COM_FIT_CHANNEL_V1':
            raise ValueError('unsupported potential-fit format')
        grid=read_array(stream,'<i4',(3,))
        ni,ns,nf,_,_=read_array(stream,'<i4',(5,))
        if ns != 6 or nf != 12 or ni != int(np.prod(grid+1)):
            raise ValueError('fit dimensions disagree')
        read_array(stream,'<i4',(2,))
        read_array(stream,'<f8',(5,))
        linear=read_array(stream,'<f8',(int(ni),6))
        modulation=read_array(stream,'<f8',(int(ni),6))
        if stream.read(1): raise ValueError('unexpected trailing fit records')
    if n_func not in (18,36): raise ValueError('n_func must be 18 or 36')
    coefficients=np.concatenate((linear,modulation),axis=1) if n_func==36 else linear
    tensor=np.zeros((int(ni),3,n_func//3,3),dtype='<f8',order='F')
    for component in range(3): tensor[:,component,:,component]=coefficients
    with output.open('xb') as stream:
        for _ in range(3): write_record(stream,b'')
        write_record(stream,tensor)
        write_record(stream,b'')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source',type=Path)
    args=parser.parse_args()
    p=load()
    export(args.source,Path(p['func_dir'])/'V_func_tmp_tmp',p['n_func'])
