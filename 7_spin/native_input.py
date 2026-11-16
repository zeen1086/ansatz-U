"""Initialize a QE grid/core/pseudopotential evaluator without SCF or diagonalization."""
from pathlib import Path
import json,struct
import numpy as np
from common import CONFIG, INTERFACE

def write_density(path,values,L):
    temp=path.with_name(path.name+'.partial')
    with temp.open('wb') as f:
        for raw in [np.array([40*L,40*L,36],dtype='<i4').tobytes(),np.asarray(values,dtype='<f8').tobytes(order='F')]:
            f.write(struct.pack('<i',len(raw)));f.write(raw);f.write(struct.pack('<i',len(raw)))
    temp.replace(path)

def read_density(path,L):
    with path.open('rb') as f:
        def rec():
            n=f.read(4);size=struct.unpack('<i',n)[0];b=f.read(size);assert len(b)==size and f.read(4)==n;return b
        assert np.array_equal(np.frombuffer(rec(),'<i4'),[40*L,40*L,36]);a=np.frombuffer(rec(),'<f8').reshape(-1,4,order='F').copy();assert not f.read(1)
    return a

def input_file(folder,L):
    meta=json.loads((INTERFACE/'V_func_tmp/V_FUNC.json').read_text())
    prim=np.array(meta['primitive_cell_A']);motif=meta['primitive_motif'];names={57:'La',28:'Ni',8:'O'}
    # Two initialization bands are allocated. The artificial startup charge is overwritten
    # by the full physical density before every v_h/v_xc call; no SCF/energy result uses it.
    atoms=[]
    for x in range(L):
        for y in range(L):
            for z,fx,fy,fz in motif:
                pos=(np.array([fx+x,fy+y,fz])@prim);atoms.append(f'{names[int(z)]} '+' '.join(f'{v:.16g}' for v in pos))
    cell=prim.copy();cell[:2]*=L;pseudo=Path(CONFIG['complete_energy']['pseudo_dir'])
    text=f'''&CONTROL calculation='scf', prefix='native_grid', outdir='{folder}/tmp', pseudo_dir='{pseudo}', disk_io='none' /
&SYSTEM ibrav=0, nat={len(atoms)}, ntyp=3, ecutwfc=70, ecutrho=280, noncolin=.true., lspinorb=.false.,
 nbnd=2, tot_charge={41*L*L-1}, occupations='smearing', smearing='gaussian', degauss=0.01,
 starting_magnetization(2)=0.1, nr1={40*L}, nr2={40*L}, nr3=36 /
&ELECTRONS startingpot='atomic', startingwfc='random' /
ATOMIC_SPECIES
La 138.90547 La.upf
Ni 58.6934 Ni.upf
O 15.999 O.upf
CELL_PARAMETERS angstrom
'''+''.join(' '.join(f'{v:.16g}' for v in row)+'\n' for row in cell)+'ATOMIC_POSITIONS angstrom\n'+'\n'.join(atoms)+'\nK_POINTS automatic\n1 1 1 0 0 0\n'
    (folder/'grid.in').write_text(text)
