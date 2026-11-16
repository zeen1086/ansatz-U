from pathlib import Path
import subprocess,json
import numpy as np
from native_input import input_file,write_density,read_density
from density_map import grid_to_cells,cells_to_grid
from common import CONFIG, compute_node
compute_node()
config=CONFIG['complete_energy'];R=Path(config['work']);old=Path(config['reference_run']);out=R/'native_check';out.mkdir(parents=True,exist_ok=False)
rho=grid_to_cells(read_density(old/'inputs/reference_actual/DENSITY.bin',2),2);idx=[(x%2)*2+y%2 for x in range(4) for y in range(4)]
write_density(out/'density.bin',cells_to_grid(rho[idx],4),4);input_file(out,4)
(out/'input_list.txt').write_text(f'1\nreference "{out}/density.bin" 0 0 1\n')
with (out/'native.out').open('w') as f:subprocess.run([str(Path(config['native_binary'])),'grid.in'],cwd=out,stdout=f,stderr=subprocess.STDOUT,check=True)
new=np.loadtxt(out/'reference_ENERGY.dat');ref=np.loadtxt(old/'results/native/reference_actual_ENERGY.dat');field=grid_to_cells(read_density(out/'reference_TRUE.bin',4),4);before=grid_to_cells(read_density(old/'results/native/reference_actual_TRUE.bin',2),2)
report=dict(charge=float(new[3]),EH_error_Ry=float(new[0]-4*ref[0]),Exc_error_Ry=float(new[1]-4*ref[1]),field_max_Ry=float(np.max(abs(field-before[idx]))),no_SCF=True)
assert abs(new[3]-656)<1e-6 and abs(report['EH_error_Ry'])<1e-6 and abs(report['Exc_error_Ry'])<1e-6 and report['field_max_Ry']<1e-7,report
(out/'PASS.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
