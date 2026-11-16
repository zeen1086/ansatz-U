"""CrI3 trial-density energy/torque, extended by full complex Dudarev N.
All energy arguments use eV; densities and fields must use the same grid/gauge.
An active Wannier band requires an explicit excluded-space contribution.
"""
import numpy as np

def trace(a,b):
    return float(np.sum(a*b.swapaxes(-1,-2)).real)

def hubbard(n,u=1.):
    assert n.shape[-2:]==(10,10) and np.max(abs(n-n.conj().swapaxes(-1,-2)))<1e-10
    return .5*u*(float(np.trace(n,axis1=-2,axis2=-1).real.sum())-trace(n,n)),u*(.5*np.eye(10)-n)

def integral(a,b,dv):
    assert a.shape==b.shape and np.isfinite(a).all() and np.isfinite(b).all() and dv>0
    return float(np.sum(a*b,dtype=np.float64)*dv)

def energy(band,excluded,rho,v,n,vu,eh,exc,exc_core,ion,dv):
    eu,_=hubbard(n)
    local=-integral(rho,v,dv);uc=-trace(n,vu)
    terms=dict(band=band,excluded=excluded,minus_rho_v=local,minus_N_VU=uc,EH=eh,Exc=exc,minus_Exc_core=-exc_core,EU=eu,ion=ion)
    return sum(terms.values()),terms

def torque(band_gradient,excluded_gradient,rho,v,n,vu,drho,dv_field,dn,dvu,native,dvol):
    # tau=-dF/dt; the nonvariational mismatch term is necessary for a fitted V1.
    _,vn=hubbard(n)
    terms=dict(band=-band_gradient,excluded=-excluded_gradient,
        input_work=integral(rho,dv_field,dvol),density_mismatch=integral(v-native,drho,dvol),
        Hubbard_input_work=trace(n,dvu),Hubbard_mismatch=trace(vu-vn,dn))
    return sum(terms.values()),terms

def selfcheck():
    rng=np.random.default_rng(1611);shape=(11,4);rho=rng.normal(size=shape);dr=rng.normal(size=shape)
    v=rng.normal(size=shape);dv=rng.normal(size=shape);a=rng.normal(size=(1,10,10))+1j*rng.normal(size=(1,10,10));n=np.eye(10)[None]*.4+.01*(a+a.conj().swapaxes(-1,-2));dn=.02*(a+a.conj().swapaxes(-1,-2))
    vu=np.eye(10)[None]*.2;dvu=np.eye(10)[None]*.07;volume=.3
    # A quadratic density functional makes the complete energy derivative independently known.
    def f(t):
        r=rho+t*dr;nn=n+t*dn
        return energy(2*t+.2*t*t,3*t,r,v+t*dv,nn,vu+t*dvu,.5*integral(r,r,volume),0,0,0,volume)[0]
    tau,_=torque(2,3,rho,v,n,vu,dr,dv,dn,dvu,rho,volume)
    h=1e-5;assert abs(tau+(f(h)-f(-h))/(2*h))<1e-8
    # Gauge consistency requires band+excluded electron count to equal trial charge.
    charge=float(rho[:,0].sum()*volume);shift=np.zeros_like(v);shift[:,0]=.17
    e1=energy(0,0,rho,v,n,vu,0,0,0,0,volume)[0]
    e2=energy(charge*.17,0,rho,v+shift,n,vu,0,0,0,0,volume)[0]
    assert abs(e1-e2)<1e-12
    print('energy/torque finite-difference and scalar-gauge selfchecks PASS')
if __name__=='__main__':selfcheck()
