"""Initialize a QE grid/core/pseudopotential evaluator without SCF or diagonalization."""
from pathlib import Path
import json,sys
import numpy as np
from common import CONFIG, INTERFACE, SOURCE
sys.path.insert(0, str(SOURCE))
from tools.fortran_io import read_array, write_record

def write_density(path,values,L):
    temp=path.with_name(path.name+'.partial')
    with temp.open('wb') as f:
        write_record(f,np.array([40*L,40*L,36],dtype='<i4'))
        write_record(f,np.asarray(values,dtype='<f8'))
    temp.replace(path)

def read_density(path,L):
    with path.open('rb') as f:
        if not np.array_equal(read_array(f,'<i4',(3,)),[40*L,40*L,36]):
            raise ValueError('Density grid differs from the requested supercell')
        a=read_array(f,'<f8',(40*L*40*L*36,4)).copy()
        if f.read(1):raise ValueError('Unexpected trailing density records')
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


def check_large() -> None:
    import subprocess
    from density_map import grid_to_cells,cells_to_grid
    from common import compute_node

    compute_node()
    config=CONFIG['complete_energy'];R=Path(config['shared_work']);reference=Path(config['reference_run']);out=R/'native_check';out.mkdir(parents=True,exist_ok=False)
    rho=grid_to_cells(read_density(reference/'inputs/reference_actual/DENSITY.bin',2),2);idx=[(x%2)*2+y%2 for x in range(4) for y in range(4)]
    write_density(out/'density.bin',cells_to_grid(rho[idx],4),4);input_file(out,4)
    (out/'input_list.txt').write_text(f'1\nreference "{out}/density.bin" 0 0 1\n')
    with (out/'native.out').open('w') as f:subprocess.run([str(Path(config['native_binary'])),'grid.in'],cwd=out,stdout=f,stderr=subprocess.STDOUT,check=True)
    new=np.loadtxt(out/'reference_ENERGY.dat');ref=np.loadtxt(reference/'results/native/reference_actual_ENERGY.dat');field=grid_to_cells(read_density(out/'reference_TRUE.bin',4),4);before=grid_to_cells(read_density(reference/'results/native/reference_actual_TRUE.bin',2),2)
    report=dict(charge=float(new[3]),EH_error_Ry=float(new[0]-4*ref[0]),Exc_error_Ry=float(new[1]-4*ref[1]),field_max_Ry=float(np.max(abs(field-before[idx]))),no_SCF=True)
    assert abs(new[3]-656)<1e-6 and abs(report['EH_error_Ry'])<1e-6 and abs(report['Exc_error_Ry'])<1e-6 and report['field_max_Ry']<1e-7,report
    (out/'PASS.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=['check'])
    parser.parse_args()
    check_large()
