"""Summarize completed J12 and check native energy curvature against full torque."""
from pathlib import Path
import csv
import json
import math

from common import CONFIG
ROOT = Path(CONFIG['complete_energy']['work'])

def richardson(a, b):
    return (4*a-b)/3

def check_curvature(energy, torque):
    ec, tc = [], []
    for j, angle in enumerate((.25, .5, 1.)):
        h = math.radians(angle)
        ec.append(1000*(energy[2*j+1]+energy[2*j+2]-2*energy[0])/h**2)
        tc.append(-1000*(torque[2*j+1]-torque[2*j+2])/(2*h))
    a, b = richardson(*ec[:2]), richardson(*tc[:2])
    return dict(energy_curvatures_meV=ec, torque_curvatures_meV=tc,
                energy_Richardson_meV=a, torque_Richardson_meV=b,
                discrepancy_meV=abs(a-b),
                energy_step_indicator_meV=abs(a-richardson(*ec[1:])),
                torque_step_indicator_meV=abs(b-richardson(*tc[1:])))

def selfcheck():
    angles = [0]+[math.radians(s*a) for a in (.25,.5,1.) for s in (1,-1)]
    out=check_curvature([2*x*x+3*x**4 for x in angles],[-4*x-12*x**3 for x in angles])
    assert abs(out['energy_Richardson_meV']-4000)<1e-7
    assert out['discrepancy_meV']<1e-7

def main():
    selfcheck()
    rows=[]
    for L in CONFIG['complete_energy']['sizes']:
        root=ROOT/f'L{L:02}'
        if not (root/'RESULT.json').exists():
            continue
        result=json.loads((root/'RESULT.json').read_text())
        energies=[json.loads((root/f'case_{i:02}/ENERGY.json').read_text()) for i in range(11)]
        torques=[r['torque'][r['sites'].index(1)] for r in result['torques']]
        checks=check_curvature([e['energy'] for e in energies],torques)
        checks['meaning']='Single-site curvature of complete functional; this is an energy/torque consistency check, not J12 or DFT validation.'
        (root/'ENERGY_TORQUE_CHECK.json').write_text(json.dumps(checks,indent=2)+'\n')
        rows.append(dict(L=L,J12_meV=result['J12_meV'],shell_spread_meV=result['spread_meV'],
                         max_step_indicator_meV=max(p['step_indicator_meV'] for p in result['pairs']),
                         reverse_pair_discrepancy_meV=max(abs(p['J12_meV']-p['reverse_J12_meV']) for p in result['pairs'] if p['reverse_J12_meV'] is not None),
                         energy_torque_discrepancy_meV=checks['discrepancy_meV'],mesh='2x2x8',inner_step_rad=.001,outer_steps_deg='0.25;0.5;1',V1_unchanged=True))
        with (root/'energy_components.csv').open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=['case','energy',*energies[0]['terms']]);writer.writeheader()
            for i,e in enumerate(energies):writer.writerow(dict(case=i,energy=e['energy'],**e['terms']))
    if rows:
        with (ROOT/'J12_by_size.csv').open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    print(json.dumps(rows,indent=2))

if __name__=='__main__':
    main()
