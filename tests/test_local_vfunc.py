"""Native RVR/H assembly versus a direct real-grid oracle, including complex full N."""
from itertools import product
from pathlib import Path
import argparse
import json
import os
import re
import struct
import subprocess
import sys
import tempfile
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from local_vfunc import write_coeff, export_background


def put(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n')


def call(binary, *args, ok=True):
    p = subprocess.run([str(binary), *map(str, args)], capture_output=True, text=True,
                       env={'OMP_NUM_THREADS':'1', **os.environ, 'VECLIB_MAXIMUM_THREADS':'1'})
    if (p.returncode == 0) != ok:
        raise AssertionError(p.stdout+p.stderr)
    return p.stdout+p.stderr


def local(point, values, lookup, spacing, cutoff):
    r = np.linalg.norm(point)
    if r >= cutoff:
        return np.zeros(values.shape[1])
    position = point/spacing; lo = np.floor(position).astype(int); f = position-lo
    answer = np.zeros(values.shape[1])
    for shift in product((0,1), repeat=3):
        key = tuple(lo+shift)
        if key in lookup:
            answer += np.prod([f[c] if shift[c] else 1-f[c] for c in range(3)])*values[lookup[key]]
    return answer*np.cos(np.pi*r/(2*cutoff))**2


def spin_matrix(fields):
    n = fields.shape[1]; h = np.zeros((2*n,2*n), complex)
    h[::2,::2] = fields[0]+fields[3]; h[1::2,1::2] = fields[0]-fields[3]
    h[::2,1::2] = fields[1]-1j*fields[2]; h[1::2,::2] = fields[1]+1j*fields[2]
    return h


def background(path, h0, projector):
    rr = np.array(list(product(range(-4,5), repeat=3)))
    h = np.zeros((729,30,30),complex); p = np.zeros((729,5,30),complex)
    zero = np.flatnonzero(np.all(rr==0, axis=1))[0]; h[zero]=h0;p[zero]=projector
    with path.open('xb') as f:
        f.write(b'LANIO2_AU_V1    ');f.write(struct.pack('<8i',16909060,729,30,5,10,810,1,0))
        f.write(np.asarray(rr.T,'<i4').tobytes(order='F'))
        for a in (h.transpose(1,2,0),p.transpose(1,2,0),np.zeros((810,810),complex),np.zeros((10,10),complex)):
            f.write(np.asarray(a,'<c16').tobytes(order='F'))
    put(Path(str(path)+'.json'),dict(kind='pure_H0_no_local_fields_no_HU',background_id='fixture-background'))


class BackgroundExportTests(unittest.TestCase):
    def test_bare_background_excludes_source_scalar_and_occupation(self):
        with tempfile.TemporaryDirectory() as tmp:
            r=Path(tmp);source=r/'primitive.bin';q=np.array(list(product(range(8),repeat=3)))/8
            h=np.diag(np.arange(30)*.1).astype(complex);p=np.eye(5,30,dtype=complex)*(.2+.1j)
            def write(f,a,dtype):f.write(np.asarray(a,dtype=dtype).tobytes(order='F'))
            with source.open('wb') as f:
                f.write(b'LANIO2_DFT_V1    ');f.write(struct.pack('<4i',57600,512,128,30))
                write(f,q.T,'<f8');write(f,np.repeat(h[:,:,None],512,axis=2),'<c16')
                write(f,np.eye(10)*.9,'<c16');write(f,np.full(57600,1234.),'<f8')
                for _ in range(512):
                    f.write(struct.pack('<i',30));write(f,p,'<c16')
                    write(f,np.eye(30),'<c16');write(f,np.eye(30),'<c16')
            report=export_background(source,r/'pure.au')
            self.assertLess(report['H0_source_roundtrip_max_eV'],1e-10)
            self.assertEqual(len(report['background_id']),64)
            self.assertEqual(set(report),{'kind','source','background_id','transformation','H0_source_roundtrip_max_eV','excluded'})
            with (r/'pure.au').open('rb') as f:
                self.assertEqual(f.read(16),b'LANIO2_AU_V1    ');f.read(32+12*729)
                hr=np.frombuffer(f.read(16*900*729),'<c16').reshape(30,30,729,order='F')
                pr=np.frombuffer(f.read(16*150*729),'<c16').reshape(5,30,729,order='F')
            rr=np.array(list(product(range(-4,5),repeat=3)))
            center=int(np.flatnonzero(np.all(rr==0,axis=1))[0])
            self.assertLess(np.max(abs(hr[:,:,center]-h)),1e-12)
            self.assertLess(np.max(abs(pr[:,:,center]-p)),1e-12)
            with self.assertRaises(FileExistsError):export_background(source,r/'pure.au')


def run(binary: Path, destination: Path, *, batch_checks: bool=False, shared_cache_checks: bool=False):
    rng = np.random.default_rng(141126)
    with tempfile.TemporaryDirectory(prefix='V1_native_', dir=destination) as tmp:
        r = Path(tmp); model=r/'V_func_tmp';wan=r/'wan';project=r/'rvr'
        for d in (model,wan,project): d.mkdir()
        ncoef=2; nodes=np.array(list(product(range(-3,4),repeat=3)), dtype=np.int32)
        values=rng.normal(scale=.1,size=(len(nodes),5)); alpha=.7;cutoff=1.1
        np.save(model/'C_v_local.npy',values[:,:3].copy());np.save(model/'C_B_local.npy',values[:,3:].copy())
        model_id='a'*64;basis_id='b'*64
        with (model/'native_geometry.bin').open('xb') as f:
            f.write(b'LANIO2_VFUNC_V2 ');f.write(struct.pack('<8i',1,2,len(nodes),2,2,2,0,0))
            f.write(np.array([3,3,3],'<i4').tobytes());f.write(np.eye(3,dtype='<f8').tobytes(order='F'))
            f.write(np.array([.5,.5,.5,.5,.5,.5,cutoff,alpha],'<f8').tobytes())
            f.write(model_id.encode());f.write(nodes.astype('<i4').tobytes())
        source=np.array([[-1,0,0],[0,0,0],[1,0,0],[0,0,1]])
        weights=np.array([.5,1.,.7,.3])
        np.savetxt(r/'translations.dat',np.column_stack([source,weights]),header='4',comments='',fmt=['%d','%d','%d','%.17g'])
        raw=(rng.normal(size=(30,4,8))+1j*rng.normal(size=(30,4,8)))/20
        for i in range(30):
            (wan/f'wan_{i+1}').write_bytes(np.asarray(raw[i],'<c16').tobytes(order='F'))
        kz=np.array([0.,.17]);np.savetxt(r/'kz.dat',kz,header='2',comments='')
        control=r/'projection.nml'
        control.write_text("&local_projection\n"+'\n'.join(f" {k}='{v}'," for k,v in dict(
            model_path=model,wan_dir=wan,translations_path=r/'translations.dat',rvr_dir=project,
            kz_path=r/'kz.dat',basis_id=basis_id).items())+"\n nwan=30, grid_weight=1.\n/\n")
        # A stale partial is regenerated; a validated completed file is reused unchanged.
        (project/'spatial_fields.partial').write_bytes(b'interrupted preparation')
        call(binary/'rvr','prepare',control)
        cache=project/'spatial_fields';stamp=cache.stat().st_mtime_ns
        assert 'SPATIAL_FIELDS_REUSED' in call(binary/'rvr','prepare',control)
        assert cache.stat().st_mtime_ns==stamp
        call(binary/'rvr','project',control,1,1)
        first=project/'B_general_1_kz1';stamp=first.stat().st_mtime_ns
        (project/'B_general_1_kz2.partial').write_bytes(b'interrupted row')
        call(binary/'rvr','project',control)
        assert first.stat().st_mtime_ns==stamp
        rows=sorted(project.glob('B_general_*'));assert len(rows)==60
        last=project/'B_general_30_kz2';saved=last.read_bytes()
        last.rename(Path(str(last)+'.partial'));Path(str(last)+'.partial').write_bytes(saved[:100])
        call(binary/'rvr','project',control,30,2)
        assert last.read_bytes()==saved
        stamps={p:p.stat().st_mtime_ns for p in rows}
        call(binary/'rvr','project',control)
        assert all(p.stat().st_mtime_ns==stamps[p] for p in rows)
        wrong=bytearray(saved);wrong[48:112]=b'x'*64;last.write_bytes(wrong)
        assert 'provenance' in call(binary/'rvr','project',control,30,2,ok=False)
        assert last.read_bytes()==wrong
        head=struct.unpack_from('<8i',saved,16)
        bad=bytearray(saved);struct.pack_into('<d',bad,192+8*(head[2]+head[5]),float('nan'))
        last.write_bytes(bad)
        assert 'Nonfinite' in call(binary/'rvr','project',control,30,2,ok=False)
        bad=bytearray(saved);struct.pack_into('<d',bad,len(bad)-8,float('nan'));last.write_bytes(bad)
        assert 'Nonfinite' in call(binary/'rvr','project',control,30,2,ok=False)
        last.write_bytes(saved[:-1])
        assert 'byte count mismatch' in call(binary/'rvr','combine',control,ok=False)
        assert cache.is_file() and all(p.is_file() for p in rows)
        assert not (project/'B_general').exists() and (project/'B_general.partial').exists()
        last.write_bytes(saved)
        # Merging only needs cache metadata and means, never the full grid payload.
        cache_bytes=cache.read_bytes();bad=bytearray(cache_bytes)
        grid_head=struct.unpack_from('<8i',bad,16)
        struct.pack_into('<d',bad,112+8*grid_head[3],float('nan'));cache.write_bytes(bad)
        assert 'SPATIAL_FIELDS_REUSED' in call(binary/'rvr','prepare',control)
        call(binary/'rvr','combine',control,'--keep-intermediates')
        assert not (project/'B_general.partial').exists()
        reference=(project/'B_general').read_bytes()
        assert cache.is_file() and all(p.is_file() for p in rows)
        # Resume after publication and an interrupted cleanup, with the cache already gone.
        cache.unlink();rows[0].unlink();(project/'unrelated.txt').write_text('preserve this file')
        stamp=(project/'B_general').stat().st_mtime_ns
        call(binary/'rvr','combine',control)
        assert (project/'B_general').read_bytes()==reference
        assert {p.name for p in project.iterdir()}=={'B_general','unrelated.txt'}
        for command in ('prepare','project','combine'):
            assert 'B_GENERAL_REUSED' in call(binary/'rvr',command,control)
        assert (project/'B_general').stat().st_mtime_ns==stamp and not cache.exists()
        published=(project/'B_general').read_bytes()
        bad=bytearray(published);bad[48:112]=b'x'*64;(project/'B_general').write_bytes(bad)
        assert 'provenance mismatch' in call(binary/'rvr','combine',control,ok=False)
        assert (project/'B_general').read_bytes()==bad
        (project/'B_general').write_bytes(published)
        h0=rng.normal(size=(30,30))+1j*rng.normal(size=(30,30));h0=(h0+h0.conj().T)/20
        pr=(rng.normal(size=(5,30))+1j*rng.normal(size=(5,30)))/15
        background(r/'pure.au',h0,pr)
        q=np.array([[.13,.21,0.,.5],[.29,-.11,.17,.5]])
        if batch_checks:
            q=np.array([[.13,.21,0.,.1],[.29,-.11,.17,.2],[-.19,.07,0.,.3],[.41,.22,.17,.4]])
        np.savetxt(r/'k.dat',q,header=str(len(q)),comments='',fmt='%.17g')
        grid=np.array([[x/2,y/2,z/2] for z in range(2) for y in range(2) for x in range(2)])
        lookup={tuple(row):i for i,row in enumerate(nodes)}
        errors=[];batch_errors=[];cache_errors=[]
        for length in (1,2,3):
            cells=np.array(list(product(range(length),repeat=2)));nc=length**2;dim=60*nc
            coeff=rng.normal(size=(nc,ncoef,4))
            n=rng.normal(size=(nc,10,10))+1j*rng.normal(size=(nc,10,10));n=(n+n.conj().swapaxes(1,2))/30+.5*np.eye(10)
            for derivative in (False,True):
                meta=dict(ncoeff=ncoef,model_id=model_id,Ueff_eV=1.,scalar_offset_eV=-.19)
                coeffpath=r/f'coeff_{length}_{int(derivative)}';write_coeff(coeffpath,meta,length,coeff,n,derivative)
                prefix=r/f'H_{length}_{int(derivative)}'
                ctl=prefix.with_suffix('.nml')
                ctl.write_text("&local_hamiltonian\n"+'\n'.join(f" {k}='{v}'," for k,v in dict(
                    rvr_path=project/'B_general',background_path=r/'pure.au',
                    background_kind='pure_H0_no_local_fields_no_HU',coefficient_path=coeffpath,
                    kpoints_path=r/'k.dat',output_prefix=prefix).items())+"\n dump_components=.true.\n/\n")
                result=call(binary/'cal_Hsk',ctl)
                if shared_cache_checks:
                    from unittest.mock import patch
                    import sys
                    sys.path.insert(0,str(ROOT/'7_spin'))
                    import complete_energy as ce
                    from torque_contract import contract_batch,derivative as read_derivative
                    old=r/'shared_inputs';(old/'RVR').mkdir(parents=True,exist_ok=True)
                    put(old/'RVR/PROJECTION.json',dict(basis_id=basis_id))
                    with patch.object(ce,'WORK',old),patch.object(ce,'SHARED',r/'shared'),patch.dict(ce.CONFIG,model_id=model_id):
                        op=[ce.cached_operators(length,k[:3],{'fixture':model_id},ce.static(r/'pure.au')) for k in q]
                    listing=r/f'cache_L{length}.txt';listing.write_text(''.join(str(d/'native.bin')+'\n' for d in op))
                    cached=r/f'cached_{length}_{int(derivative)}';cached_ctl=cached.with_suffix('.nml')
                    cached_ctl.write_text(ctl.read_text().replace(str(prefix),str(cached)).replace(' dump_components=',
                        f" background_cache_list='{listing}',\n dump_components="))
                    call(binary/'cal_Hsk',cached_ctl)
                    for ki in range(len(q)):
                        for suffix in ('h','h0','hv','hb','hu'):
                            direct=np.fromfile(f'{prefix}.k{ki+1:06}.{suffix}',dtype='<c16')
                            reused=np.fromfile(f'{cached}.k{ki+1:06}.{suffix}',dtype='<c16')
                            cache_errors.append(float(np.max(abs(direct-reused))))
                            assert cache_errors[-1]<2e-10,cache_errors[-1]
                    if derivative:
                        d=read_derivative(coeffpath)
                        for ki,k in enumerate(q):
                            rho=np.eye(dim)[None]*.2
                            a=np.load(op[ki]/'native_projectors.npy')
                            direct,_,_=contract_batch(project/'B_general',r/'pure.au',k[None,:3],rho,[d])
                            reused,_,_=contract_batch(project/'B_general',r/'pure.au',k[None,:3],rho,[d],projectors_by_point=[a])
                            np.testing.assert_allclose(reused,direct,atol=2e-10,rtol=0)
                    if length==1 and not derivative:
                        payload=(op[0]/'native.bin').read_bytes()
                        for tag,content in [('truncated',payload[:-1]),('wrong_k',payload[:24]+struct.pack('<d',.99)+payload[32:])]:
                            (op[0]/'native.bin').write_bytes(content)
                            badctl=r/f'bad_cache_{tag}.nml'
                            badctl.write_text(cached_ctl.read_text().replace(str(cached),str(r/tag)))
                            assert 'Cached background' in call(binary/'cal_Hsk',badctl,ok=False)
                        (op[0]/'native.bin').write_bytes(payload)
                if batch_checks:
                    passes=re.findall(r'RVR_KZ_BATCH sector, points, blocks=\s*(\d+)\s+(\d+)\s+(\d+)',result)
                    assert passes==[('1','2','900'),('2','2','900')],passes
                index=np.loadtxt(f'{prefix}.index',skiprows=4)
                np.testing.assert_array_equal(index[:,0],np.arange(1,len(q)+1))
                np.testing.assert_allclose(index[:,1:5],q,rtol=0,atol=1e-15)
                points=np.concatenate([grid+[*cell,0] for cell in cells])
                field=np.zeros((len(points),4))
                for g,point in enumerate(points):
                    for site,cell in enumerate(cells):
                        # Direct periodic potential sum, independent of RVR folding/maps.
                        for tx,ty,tz in product(range(-2,3),range(-2,3),range(-2,3)):
                            delta=point-(np.array([*cell,0])+.5)-[tx*length,ty*length,tz]
                            v=local(delta,values,lookup,.5,cutoff)
                            field[g,0]+=v[:2]@coeff[site,:,0]+(0 if derivative else alpha*v[2])
                            field[g,1:]+=v[3:]@coeff[site,:,1:]
                field[:,0]-=field[:,0].mean()
                for ki,k in enumerate(q[:,:3]):
                    wave=np.zeros((len(points),30*nc),complex)
                    for ci,cell in enumerate(cells):
                        for sr,ww in zip(source,weights):
                            center=cell-sr[:2];wrapped=center%length;boundary=(center-wrapped)//length
                            site=wrapped[0]*length+wrapped[1]
                            phase=np.exp(2j*np.pi*(boundary@k[:2]-sr[2]*k[2]))
                            wave[ci*8:(ci+1)*8,site*30:(site+1)*30]+=phase*ww*raw[:,np.flatnonzero(np.all(source==sr,axis=1))[0],:].T
                    projected=np.array([wave.conj().T@(field[:,c,None]*wave) for c in range(4)])
                    hv=spin_matrix(np.array([projected[0],*np.zeros_like(projected[1:])]))
                    hb=spin_matrix(np.array([np.zeros_like(projected[0]),*projected[1:]]))
                    hu=np.zeros((dim,dim),complex)
                    for site in range(nc):
                        a=np.zeros((10,dim),complex);a[:5,site*60:(site+1)*60:2]=pr;a[5:,site*60+1:(site+1)*60:2]=pr
                        v=-n[site]+(0 if derivative else .5*np.eye(10));hu+=a.conj().T@v@a
                    bg=np.zeros((dim,dim),complex) if derivative else np.kron(np.eye(nc),np.kron(h0,np.eye(2)))
                    if not derivative: hv-=.19*np.eye(dim)
                    expected={'h0':bg,'hv':hv,'hb':hb,'hu':hu,'h':bg+hv+hb+hu}
                    for suffix,want in expected.items():
                        actual=np.fromfile(f'{prefix}.k{ki+1:06}.{suffix}',dtype='<c16').reshape(dim,dim,order='F')
                        error=float(np.max(abs(actual-want)));errors.append(error)
                        assert error<2e-10,(length,derivative,ki,suffix,error)
                    if batch_checks and length==2:
                        folder=r/f'single_{int(derivative)}_{ki}';folder.mkdir()
                        single=folder/'H';kfile=folder/'k.dat'
                        np.savetxt(kfile,np.r_[k,1.][None],header='1',comments='',fmt='%.17g')
                        single.with_suffix('.nml').write_text(ctl.read_text().replace(str(prefix),str(single))
                            .replace(str(r/'k.dat'),str(kfile)))
                        call(binary/'cal_Hsk',single.with_suffix('.nml'))
                        for suffix in expected:
                            grouped=np.fromfile(f'{prefix}.k{ki+1:06}.{suffix}',dtype='<c16')
                            separate=np.fromfile(f'{single}.k000001.{suffix}',dtype='<c16')
                            batch_errors.append(float(np.max(abs(grouped-separate))))
                            assert batch_errors[-1]<2e-12
        # Reject a prediction from another frozen model before assembly.
        bad=r/'wrong_coeff';contents=bytearray(coeffpath.read_bytes());contents[48:112]=b'c'*64;bad.write_bytes(contents)
        ctl.write_text(ctl.read_text().replace(str(coeffpath),str(bad)).replace(str(prefix),str(r/'bad_output')))
        rejection=call(binary/'cal_Hsk',ctl,ok=False)
        assert 'model mismatch' in rejection
        result=dict(status='PASS',L=[1,2,3],k_points=q[:,:3].tolist(),
            real_and_derivative_components=len(errors),maximum_matrix_error_eV=max(errors),
            full_complex_N=True,scalar_mean_same_basis=True,periodic_z_phases=True,
            wrong_model_rejected=True,compact_rvr_matches_debug=True,failed_combine_retains_intermediates=True,
            interrupted_stages_resumed=True,completed_outputs_not_rewritten=True,merge_skips_grid_payload=True,
            scope='synthetic full native chain versus direct periodic real-grid integration',
            actual_material_Hamiltonian_validation=False)
        if shared_cache_checks:
            result.update(shared_cache_components=len(cache_errors),shared_cache_max_eV=max(cache_errors))
        if batch_checks:
            result.update(batch_vs_single_components=len(batch_errors),batch_vs_single_max_eV=max(batch_errors),
                          RVR_reads_per_sector=900,interleaved_input_order_preserved=True)
        put(destination/'NATIVE_CHAIN_TEST.json',result)
        return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--shared-cache-checks',action='store_true')
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    print(json.dumps(run(args.binary,args.output,batch_checks=True,shared_cache_checks=args.shared_cache_checks),indent=2))
