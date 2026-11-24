"""Complete trial energy, Hubbard terms and metric-consistent density response.
Energy arguments use eV; densities and fields must share a grid and gauge.
Active-space band energies require an explicit excluded-space contribution.
"""
from pathlib import Path
import argparse
import csv
import json
import math
import sys
import numpy as np
from scipy.special import erfc
from common import CONFIG, SOURCE, RY, SIGMA, quality_check
sys.path.insert(0, str(SOURCE))

def trace(a,b):
    return float(np.sum(a*b.swapaxes(-1,-2)).real)

def hubbard(n,u):
    assert n.shape[-2:]==(10,10) and np.max(abs(n-n.conj().swapaxes(-1,-2)))<1e-10
    return .5*u*(float(np.trace(n,axis1=-2,axis2=-1).real.sum())-trace(n,n)),u*(.5*np.eye(10)-n)

def integral(a,b,dv):
    assert a.shape==b.shape and np.isfinite(a).all() and np.isfinite(b).all() and dv>0
    return float(np.sum(a*b,dtype=np.float64)*dv)

def energy(band,excluded,rho,v,n,vu,eh,exc,exc_core,ion,dv,u):
    eu,_=hubbard(n,u)
    local=-integral(rho,v,dv);uc=-trace(n,vu)
    terms=dict(band=band,excluded=excluded,minus_rho_v=local,minus_N_VU=uc,EH=eh,Exc=exc,minus_Exc_core=-exc_core,EU=eu,ion=ion)
    return sum(terms.values()),terms

def torque(band_gradient,excluded_gradient,rho,v,n,vu,drho,dv_field,dn,dvu,native,dvol,u):
    # tau=-dF/dt; the nonvariational mismatch term is necessary for a fitted V1.
    _,vn=hubbard(n,u)
    terms=dict(band=-band_gradient,excluded=-excluded_gradient,
        input_work=integral(rho,dv_field,dvol),density_mismatch=integral(v-native,drho,dvol),
        Hubbard_input_work=trace(n,dvu),Hubbard_mismatch=trace(vu-vn,dn))
    return sum(terms.values()),terms

def energy_selfcheck():
    rng=np.random.default_rng(1611);shape=(11,4);rho=rng.normal(size=shape);dr=rng.normal(size=shape)
    v=rng.normal(size=shape);dv=rng.normal(size=shape);a=rng.normal(size=(1,10,10))+1j*rng.normal(size=(1,10,10));n=np.eye(10)[None]*.4+.01*(a+a.conj().swapaxes(-1,-2));dn=.02*(a+a.conj().swapaxes(-1,-2))
    vu=np.eye(10)[None]*.2;dvu=np.eye(10)[None]*.07;volume=.3;u=1.7
    # A quadratic density functional makes the complete energy derivative independently known.
    def f(t):
        r=rho+t*dr;nn=n+t*dn
        return energy(2*t+.2*t*t,3*t,r,v+t*dv,nn,vu+t*dvu,.5*integral(r,r,volume),0,0,0,volume,u)[0]
    tau,_=torque(2,3,rho,v,n,vu,dr,dv,dn,dvu,rho,volume,u)
    h=1e-5;assert abs(tau+(f(h)-f(-h))/(2*h))<1e-8
    # Gauge consistency requires band+excluded electron count to equal trial charge.
    charge=float(rho[:,0].sum()*volume);shift=np.zeros_like(v);shift[:,0]=.17
    e1=energy(0,0,rho,v,n,vu,0,0,0,0,volume,u)[0]
    e2=energy(charge*.17,0,rho,v+shift,n,vu,0,0,0,0,volume,u)[0]
    assert abs(e1-e2)<1e-12
    print('energy/torque finite-difference and scalar-gauge selfchecks PASS')


def static(path: Path):
    with path.open('rb') as f:
        assert f.read(16)==b'LANIO2_AU_V1    ';head=np.fromfile(f,'<i4',8);nr=int(head[1]);assert head[2]==30
        rr=np.fromfile(f,'<i4',3*nr).reshape(3,nr,order='F').T
        hr=np.fromfile(f,'<c16',900*nr).reshape(30,30,nr,order='F').transpose(2,0,1)
        pr=np.fromfile(f,'<c16',150*nr).reshape(5,30,nr,order='F').transpose(2,0,1)
    return rr,hr,pr

def operators(L,k,rr,hr,pr):
    cells=np.array([(x,y) for x in range(L) for y in range(L)]);n=L*L;q=np.column_stack([(cells+k[:2])/L,np.full(n,k[2])])
    phase=np.exp(2j*np.pi*q@rr.T);hq=(phase@hr.reshape(len(rr),-1)).reshape(n,30,30);pq=(phase@pr.reshape(len(rr),-1)).reshape(n,5,30)
    h0=np.zeros((60*n,60*n),complex);aa=[]
    for i in range(n):h0[60*i:60*(i+1),60*i:60*(i+1)]=np.kron(hq[i],np.eye(2))
    for cell in cells:
        p=(pq*np.exp(2j*np.pi*(q[:,:2]@cell))[:,None,None]/L).transpose(1,0,2).reshape(5,30*n)
        a=np.zeros((10,60*n),complex);a[:5,::2]=p;a[5:,1::2]=p;aa.append(a)
    return h0,np.array(aa)

def metric(psi):
    n=len(psi);S=np.zeros((60*n,60*n),complex);T=S.copy();err=0.
    for i,w in enumerate(psi):
        s=w.conj().T@w;e,u=np.linalg.eigh(s)
        if not np.isfinite(e).all() or e.min()<=0:raise ValueError('Nonpositive or nonfinite Wannier overlap')
        quality_check(e.min()>.1, f'Small Wannier overlap eigenvalue: {e.min()}')
        t=(u/np.sqrt(e))@u.conj().T;sl=slice(60*i,60*(i+1));S[sl,sl]=np.kron(s,np.eye(2));T[sl,sl]=np.kron(t,np.eye(2));err=max(err,float(np.max(abs(s-np.eye(30)))))
    assert np.max(abs(T@S@T-np.eye(60*n)))<1e-10
    return S,T,err

def fermi(e,mu):
    x=(e-mu)/SIGMA;f=.5*erfc(x);fp=-np.exp(-x*x)/(np.sqrt(np.pi)*SIGMA)
    return f,fp,float(np.sum(e*f-SIGMA*np.exp(-x*x)/(2*np.sqrt(np.pi))))

def loewner(e,f,fp):
    de=e[:,None]-e[None,:];df=f[:,None]-f[None,:];mask=abs(de)>1e-9
    k=np.empty_like(de);np.divide(df,de,out=k,where=mask);k[~mask]=np.broadcast_to((fp[:,None]+fp[None,:])/2,k.shape)[~mask]
    return k

def adjoint(e,u,T,G,A,mu,chemical_adjoint):
    f,fp,_=fermi(e,mu);a=u.conj().T@T@A@T@u
    a.flat[::len(a)+1]-=chemical_adjoint
    Z=T@u@(loewner(e,f,fp)*a)@u.conj().T@T
    return G.conj().T@Z@G

def selfcheck():
    from scipy.optimize import brentq
    rng=np.random.default_rng(1116);n=8;a=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n));h=(a+a.conj().T)/5;a=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n));dh=(a+a.conj().T)/5;a=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n));A=(a+a.conj().T)/2
    def read(x):
        e,u=np.linalg.eigh(x);mu=brentq(lambda m:fermi(e,m)[0].sum()-3,-10,10);f,fp,_=fermi(e,mu);return e,u,mu,f,fp,(u*f)@u.conj().T
    e,u,mu,f,fp,p=read(h);alpha=np.dot(fp,np.diag(u.conj().T@A@u).real)/fp.sum();Z=adjoint(e,u,np.eye(n),np.eye(n),A,mu,alpha);eps=1e-5
    fd=trace(A,(read(h+eps*dh)[-1]-read(h-eps*dh)[-1])/(2*eps));exact=trace(Z,dh)
    assert abs(fd-exact)<1e-7,(fd,exact);assert abs(np.trace(Z))<1e-10
    print('fixed-electron full feedback density adjoint finite-difference PASS')

def full_functional_selfcheck():
    from scipy.optimize import brentq
    from density_map import density,project
    rng=np.random.default_rng(501);w,_=np.linalg.qr(rng.normal(size=(25,6))+1j*rng.normal(size=(25,6)))
    w=w@np.diag(np.linspace(.98,1.02,6));S=np.kron(w.conj().T@w,np.eye(2));s,v=np.linalg.eigh(S);T=(v/np.sqrt(s))@v.conj().T
    a=(rng.normal(size=(10,12))+1j*rng.normal(size=(10,12)))/20
    h0=rng.normal(size=(12,12))+1j*rng.normal(size=(12,12));h0=(h0+h0.conj().T)/8
    f0=rng.normal(size=(25,4))/10;df=rng.normal(size=(25,4))/10;u0=np.eye(10)*.2;du=np.eye(10)*.07
    dH=project(w,df)+a.conj().T@du@a;bg=np.ones((25,4))*.001;bg[:,1:]=0
    def run(theta):
        V=project(w,f0+theta*df)+a.conj().T@(u0+theta*du)@a;H=h0+V;e,u=np.linalg.eigh(T@H@T);mu=brentq(lambda m:fermi(e,m)[0].sum()-5,-10,10);Q=T@u;f,fp,F=fermi(e,mu);P=(Q*f)@Q.conj().T
        rho=density(w,P)+bg;N=(a@P@a.conj().T)[None];EU,VU=hubbard(N,1.);energy=F-trace(P,V)+.5*np.sum(rho*rho)+EU
        A=V-project(w,rho)-a.conj().T@VU[0]@a;ae=Q.conj().T@A@Q;alpha=float(fp@np.diag(ae).real/fp.sum());Z=adjoint(e,u,T,np.eye(12),A,mu,alpha)
        return energy,trace(Z,dH)
    eps=1e-5;fd=-(run(eps)[0]-run(-eps)[0])/(2*eps);tau=run(0)[1]
    assert abs(fd-tau)<1e-7,(fd,tau)
    print('complete energy/torque with nonorthogonal density, complex N and fixed number PASS')


def richardson(a, b):
    return (4*a-b)/3

def curvature_result(definition, energy, torque):
    from tools.spin_coordinates import curvature_cases
    steps=definition['outer_steps_rad'];sites=definition['sites']
    if definition['angle_configurations'] != curvature_cases(sites,steps):
        raise ValueError('Curvature probe order differs from the fixed chart stencil')
    energy=np.asarray(energy,dtype=float);torque=np.asarray(torque,dtype=float)
    if (energy.shape!=(25,) or torque.shape!=(25,4)
            or not np.isfinite(energy).all() or not np.isfinite(torque).all()):
        raise ValueError('K4 requires 25 finite energies and four torques per case')
    h=np.asarray(steps)
    t=torque[1:].reshape(3,4,2,4)
    ks=(-(t[:,:,0,:]-t[:,:,1,:])/(2*h[:,None,None])).transpose(0,2,1)
    e=energy[1:].reshape(3,4,2)
    diagonal=(e[:,:,0]+e[:,:,1]-2*energy[0])/h[:,None]**2
    k=richardson(ks[0],ks[1]);ed=richardson(diagonal[0],diagonal[1])
    return dict(kind='K4',coordinates='fixed_exponential_chart',sites=sites,
                modes=[[s,a] for s in sites for a in (1,2)],outer_steps_rad=steps,
                inner_step_rad=definition['inner_step_rad'],F0_eV=float(energy[0]),
                tau0_eV_per_rad=torque[0].tolist(),K_eV_per_rad2=k.tolist(),
                K_steps_eV_per_rad2=ks.tolist(),
                step_indicator_eV_per_rad2=abs(k-richardson(ks[1],ks[2])).tolist(),
                antisymmetric_residual_eV_per_rad2=(k-k.T).tolist(),
                max_reciprocity_error_eV_per_rad2=float(np.max(abs(k-k.T))),
                energy_diagonal_steps_eV_per_rad2=diagonal.tolist(),
                energy_diagonal_eV_per_rad2=ed.tolist(),
                energy_torque_diagonal_error_eV_per_rad2=abs(ed-np.diag(k)).tolist(),
                energy_diagonal_step_indicator_eV_per_rad2=abs(ed-richardson(diagonal[1],diagonal[2])).tolist(),
                energies_eV=energy.tolist(),torques_eV_per_rad=torque.tolist(),
                coverage=np.ones((4,4),dtype=bool).tolist(),symmetrized=False)

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

def curvature_selfcheck():
    angles = [0]+[math.radians(s*a) for a in (.25,.5,1.) for s in (1,-1)]
    out=check_curvature([2*x*x+3*x**4 for x in angles],[-4*x-12*x**3 for x in angles])
    assert abs(out['energy_Richardson_meV']-4000)<1e-7
    assert out['discrepancy_meV']<1e-7

def summary() -> None:
    curvature_selfcheck()
    output = Path(CONFIG['complete_energy']['work'])
    rows=[];curvatures=[]
    for L in CONFIG['complete_energy']['sizes']:
        root=output/f'L{L:02}'
        if not (root/'RESULT.json').exists():
            continue
        result=json.loads((root/'RESULT.json').read_text())
        if result.get('kind')=='K4':
            labels=','.join(f'site{s}_a{a}' for s,a in result['modes'])
            np.savetxt(root/'K4_eV_per_rad2.csv',result['K_eV_per_rad2'],delimiter=',',header=labels)
            curvatures.append(dict(L=L,modes=result['modes'],K_eV_per_rad2=result['K_eV_per_rad2'],
                                   max_reciprocity_error_eV_per_rad2=result['max_reciprocity_error_eV_per_rad2'],
                                   physical_validation=False))
            continue
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
        with (output/'J12_by_size.csv').open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    print(json.dumps(rows+curvatures,indent=2))


if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Summarize complete-energy results or run offline checks.')
    parser.add_argument('command', choices=('summary', 'selfcheck'), nargs='?', default='selfcheck')
    args = parser.parse_args()
    if args.command == 'summary':
        summary()
    else:
        energy_selfcheck()
        selfcheck()
        full_functional_selfcheck()
        curvature_selfcheck()
