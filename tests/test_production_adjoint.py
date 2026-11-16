"""Compare one-pass torque contraction against native explicit dH on known-answer data."""
from pathlib import Path
import argparse
import json
import re
import sys
import tempfile
import numpy as np
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'7_spin'))
import torque_contract
from torque_contract import contract_batch,derivative


def main(args):
    sys.path.insert(0,str(args.native_tests))
    import test_local_vfunc as native
    original=native.call;errors=[];batch_sizes=[];rng=np.random.default_rng(816)
    def compare(cfg,points,densities,d):
        nd=densities.shape[1]
        for kz in np.unique(points[:,2]):
            group=np.flatnonzero(points[:,2]==kz)
            read=torque_contract.read_array
            with patch.object(torque_contract,'read_array',wraps=read) as observed:
                values,_,stats=contract_batch(Path(cfg['rvr_path']),Path(cfg['background_path']),
                                              points[group],densities[group],[d],workers=4)
            # Two complex reads per RVR block and one complex read per small projector file.
            assert sum(c.args[1]=='<c16' for c in observed.call_args_list)==1800+len(group)
            assert stats['rvr_blocks']==900
            batch_sizes.append(len(group))
            for slot,ki in enumerate(group):
                h=np.fromfile(cfg['output_prefix']+f'.k{ki+1:06}.h','<c16').reshape(nd,nd,order='F')
                error=float(abs(values[slot,0]+np.einsum('ij,ji->',densities[ki],h)))
                if error>2e-10:raise AssertionError((d['L'],ki,error))
                errors.append(error)
    def checked(binary,*args,**kw):
        text=original(binary,*args,**kw)
        if Path(binary).name!='cal_Hsk' or not kw.get('ok',True):return text
        ctl=Path(args[0]);cfg=dict(re.findall(r"(\w+)='([^']*)'",ctl.read_text()))
        coeff=Path(cfg['coefficient_path'])
        with coeff.open('rb') as f:f.seek(16);head=np.fromfile(f,'<i4',8)
        if head[4]!=1:return text
        d=derivative(coeff);nd=60*d['L']**2;points=np.loadtxt(cfg['kpoints_path'],skiprows=1,ndmin=2)[:,:3]
        matrices=rng.normal(size=(len(points),nd,nd))+1j*rng.normal(size=(len(points),nd,nd))
        densities=(matrices+matrices.conj().swapaxes(1,2))/(2*nd)
        compare(cfg,points,densities,d)
        if d['L']==2:
            content=bytearray(Path(cfg['background_path']).read_bytes());offset=48+12*729+16*900*729
            pr=np.frombuffer(content,dtype='<c16',count=150*729,offset=offset).reshape(5,30,729,order='F')
            for index in [365,373,445]:pr[:,:,index]+=(rng.normal(size=(5,30))+1j*rng.normal(size=(5,30)))*.01
            bg=ctl.parent/'nonlocal.au';bg.write_bytes(content)
            other=dict(cfg,background_path=str(bg),output_prefix=str(ctl.parent/'nonlocal_derivative'))
            newctl=ctl.parent/'nonlocal.nml'
            newctl.write_text(ctl.read_text().replace(cfg['background_path'],other['background_path']).replace(cfg['output_prefix'],other['output_prefix']))
            original(binary,newctl)
            compare(other,points,densities,d)
        return text
    native.call=checked
    with tempfile.TemporaryDirectory() as tmp:native.run(args.binary,Path(tmp),batch_checks=True)
    r=args.real_fixture;d=[derivative(r/'output'/f'dsite_{i:04}.bin') for i in [1,2]]
    nd=60*d[0]['L']**2;rho=np.eye(nd)*.2;k=np.array([.13,.21,0.])
    actual,_,_=contract_batch(r/'B_general',r/'pure.au',k[None],rho[None],d)
    for j,i in enumerate([1,2]):
        h=np.fromfile(r/'output'/f'dsite_{i:04}.k000001.h','<c16').reshape(nd,nd,order='F')
        errors.append(float(abs(actual[0,j]+np.einsum('ij,ji->',rho,h))))
    mixed,_,_=contract_batch(r/'B_general',r/'pure.au',np.array([k,k]),
                              np.array([rho,2*rho]),d,derivatives_by_point=[d,d[::-1]])
    errors.append(float(np.max(abs(mixed[0]-actual[0]))))
    errors.append(float(np.max(abs(mixed[1]-2*actual[0,::-1]))))
    if max(errors)>2e-10 or max(batch_sizes)<2:raise AssertionError((errors,batch_sizes))
    result=dict(status='PASS',checks=len(errors),maximum_torque_error_eV=float(max(errors)),
        checked_batch_sizes=batch_sizes,observed_rvr_blocks_per_batch=900,
        scope='exact adjoint versus native explicit dH; L=1,2,3, scalar mean, complex N, nonlocal projectors and Bloch phases; real V1 feedback with synthetic RVR')
    args.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ['native-tests','binary','real-fixture','output']:p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
