"""Real exported V1 -> one-way prediction -> native H/dH with explicit synthetic RVR."""
from pathlib import Path
import argparse
import json
import struct
import subprocess
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(Path(__file__).resolve().parent)]
from local_vfunc import metadata
from test_local_vfunc import background,spin_matrix,put


def main(args):
    meta=metadata(args.model);r=args.output;r.mkdir(parents=True,exist_ok=False)
    np_=meta['ncoeff'];basis='d'*64
    scalar=np.arange(1,np_+1)*1e-5;magnetic=-scalar/3;bias=.02
    kernel=np.r_[scalar,bias,magnetic]
    with (r/'B_general').open('wb') as f:
        f.write(b'LANIO2_LRVR_V2  ');f.write(struct.pack('<8i',1,np_,1,30,1,1,0,0))
        f.write(meta['model_id'].encode());f.write(basis.encode())
        f.write(np.asarray(meta['primitive_cell_A'],'<f8').tobytes(order='F'))
        alpha=json.loads((args.model/'KERNELS.json').read_text())['scalar_bias_coefficient_per_Ni']
        f.write(struct.pack('<d',alpha));f.write(np.zeros(np_+1,dtype='<f8').tobytes())
        f.write(np.zeros(4,dtype='<i4').tobytes());f.write(struct.pack('<d',0.))
        for i in range(30):
            for j in range(30):
                f.write(np.array([i==j],dtype='<c16').tobytes())
                f.write(np.asarray(kernel*(i==j),'<c16').tobytes())
    h0=np.diag(np.arange(30)*.1).astype(complex);a=np.eye(5,30,dtype=complex)*(.2+.1j)
    background(r/'pure.au',h0,a)
    put(r/'PROJECTION.json',dict(basis_id=basis,background_id='fixture-background',
                                 scope='synthetic native RVR for real-model runtime test'))
    np.savetxt(r/'kpoints.dat',[[.13,.21,0.,1.]],header='1',comments='')
    command=[sys.executable,'-B',str(ROOT/'local_vfunc.py'),'cal_Hsk','--model',str(args.model),
        '--structure',str(args.structure),'--output',str(r/'output'),'--rvr',str(r/'B_general'),
        '--background',str(r/'pure.au'),'--kpoints',str(r/'kpoints.dat'),'--binary',str(args.binary),
        '--derivative-sites','1','2','--components']
    subprocess.run(command+['--keep-intermediates'],check=True)
    d=np.load(r/'output/prediction.npz');nc=len(d['p']);nd=60*nc
    ap=np.zeros((10,60),complex);ap[:5,::2]=a;ap[5:,1::2]=a
    errors={}
    for index,label in enumerate(['value','dsite_0001','dsite_0002']):
        p=d['p'] if index==0 else d['dP'][index-1]
        n=d['N'] if index==0 else d['dN'][index-1]
        expected=np.zeros((nd,nd),complex)
        for site in range(nc):
            channels=np.r_[scalar@p[site,:,0]+(alpha*bias if index==0 else 0), magnetic@p[site,:,1:]]
            block=spin_matrix(channels[:,None,None]*np.eye(30)[None])
            u=meta['Ueff_eV']*(-n[site]+(.5*np.eye(10) if index==0 else 0))
            block+=ap.conj().T@u@ap
            if index==0:block+=np.kron(h0,np.eye(2))+meta['scalar_offset_eV']*np.eye(60)
            expected[site*60:(site+1)*60,site*60:(site+1)*60]=block
        actual=np.fromfile(r/'output'/f'{label}.k000001.h','<c16').reshape(nd,nd,order='F')
        errors[label]=float(np.max(abs(expected-actual)))
    assert max(errors.values())<1e-10,errors
    command[command.index('--output')+1]=str(r/'compact')
    subprocess.run(command,check=True)
    report=json.loads((r/'compact/HAMILTONIAN.json').read_text())
    matrices=sorted(p.name for p in (r/'output').glob('*.k*'))
    assert len(matrices)==15
    assert {p.name for p in (r/'compact').iterdir()}==set(matrices)|{'HAMILTONIAN.json'}
    assert all((r/'compact'/name).read_bytes()==(r/'output'/name).read_bytes() for name in matrices)
    for label in ('value','dsite_0001','dsite_0002'):
        assert report['native_indices'][label]==(r/'output'/f'{label}.index').read_text().splitlines()
    assert report['prediction_rule']=='N_dir_to_v_B'
    result=dict(status='PASS',errors_eV=errors,real_V1=True,prediction_rule=report['prediction_rule'],
                native_assembly=True,projection_basis='synthetic known-answer fixture',
                compact_vs_debug_identical_matrices=len(matrices),compact_indices_preserved=True,
                actual_material_Hamiltonian_validation=False)
    put(r/'RESULT.json',result);print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('model','structure','output','binary'):p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
