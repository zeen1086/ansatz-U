"""Frozen V1 parameters and single-response predictions for the native local RVR chain."""
from __future__ import annotations

import argparse
from dataclasses import replace
import uuid
import importlib
from itertools import product
import json
from pathlib import Path
import shutil
import struct
import subprocess
import sys

import numpy as np
from tools.spin_coordinates import CHART_KEYS, chart_spins, derivative_files

SCHEMA = 'LANIO2_LOCAL_VFUNC_V2'
PREDICTION_RULE = 'N_single_response_to_v_B'
LOCAL_CODE = Path(__file__).resolve().parent/'4_potential/local'
RUNTIME = ('features.py', 'symmetry.py', 'model.py', 'workflow.py', 'predict.py')
DATA = ('C_v_local.npy', 'C_B_local.npy', 'local_grid.npz', 'KERNELS.json',
        'model/model.json', 'model/parameters.npz')



def save(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def runtime_path(directory: Path, meta: dict) -> Path:
    if meta.get('prediction_rule') != PREDICTION_RULE:
        raise ValueError('Unsupported model prediction rule')
    code = (directory/'runtime').resolve()
    if not all((code/name).is_file() for name in RUNTIME):
        raise ValueError('Incomplete model runtime')
    return code


def runtime(directory: Path, meta: dict | None = None):
    code = runtime_path(directory, metadata(directory) if meta is None else meta)
    # Frozen modules use their original absolute imports; never silently mix versions.
    for name in (Path(f).stem for f in RUNTIME):
        old = sys.modules.get(name)
        if old is not None and Path(old.__file__).resolve().parent != code:
            raise RuntimeError(f'Conflicting model module {name}; use a fresh process')
    sys.path.insert(0, str(code))
    return importlib.import_module('predict'), importlib.import_module('features')


def export(source: Path, code: Path, structure: Path, output: Path,
           scalar_offset_eV: float) -> dict:
    """Copy frozen weights and package single-response inference with native geometry."""
    source, code, output = source.resolve(), code.resolve(), output.absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(f'Refusing to replace existing model: {output}')
    if not np.isfinite(scalar_offset_eV):
        raise ValueError('Finite, fixed scalar energy offset required')
    meta = json.loads((source/'model/model.json').read_text())
    if meta['schema'] not in ('LANIO2_ACTIVE_LOCAL_GRID_FEEDBACK_V2', 'LANIO2_ACTIVE_LOCAL_GRID_SINGLE_RESPONSE_V1'):
        raise ValueError('Unsupported source model schema')
    with np.load(structure, allow_pickle=False) as z:
        cell, positions, species, ni = (z[k].copy() for k in ('cell','positions','species','ni'))
    if len(ni) != 4 or cell.shape != (3, 3):
        raise ValueError('Export geometry must be the original four-Ni cell')
    primitive = cell.copy(); primitive[:2] /= 2
    if not np.allclose(primitive, np.diag(np.diag(primitive)), rtol=0, atol=1e-10) or np.any(np.diag(primitive) <= 0):
        raise ValueError('Native local branch requires the orthorhombic LxLx1 geometry')
    fractional = positions @ np.linalg.inv(primitive)
    motif = sorted({(int(s), *np.round(f % 1, 10)) for s, f in zip(species, fractional)})
    centers = np.mod(fractional[ni], 1)
    if not np.allclose(centers, centers[0], atol=1e-10):
        raise ValueError('One equivalent Ni per primitive cell required')
    cv = np.load(source/'C_v_local.npy', mmap_mode='r')
    cb = np.load(source/'C_B_local.npy', mmap_mode='r')
    with np.load(source/'local_grid.npz') as z:
        nodes, spacing, half = (z[k].copy() for k in ('nodes','spacing_A','half'))
    kernel = json.loads((source/'KERNELS.json').read_text())
    ncoef = len(meta['p_scales']); grid = np.rint(np.diag(primitive)/spacing).astype(int)
    if (cv.shape != (len(nodes), ncoef+1) or cb.shape != (len(nodes), ncoef)
            or not cv.flags.c_contiguous or not cb.flags.c_contiguous
            or cv.dtype != np.float64 or cb.dtype != np.float64
            or not np.allclose(grid*spacing, np.diag(primitive), atol=1e-10)):
        raise ValueError('Local kernel shape, storage or commensurate grid mismatch')
    if kernel['cutoff_A'] != 6. or meta['scalar_gauge'] != 'cell_mean_zero':
        raise ValueError('Unsupported frozen kernel/gauge definition')
    identity = {f: {'bytes':(source/f).stat().st_size} for f in DATA}
    runtime_files = {f: (LOCAL_CODE if f in ('model.py', 'workflow.py', 'predict.py') else code)/f for f in RUNTIME}
    identity.update({f'runtime/{f}': {'bytes':path.stat().st_size} for f,path in runtime_files.items()})
    contract = dict(schema=SCHEMA, prediction_rule=PREDICTION_RULE, files=identity, primitive_cell_A=primitive.tolist(),
        Ni_fractional=centers[0].tolist(), primitive_motif=[list(m) for m in motif],
        primitive_grid=grid.tolist(), scalar_offset_eV=scalar_offset_eV,
        scalar_gauge='subtract full-cell grid mean, then add fixed scalar_offset_eV',
        Ueff_eV=meta['ueff'], ncoeff=ncoef, spin_order='up d1..5, down d1..5',
        Hamiltonian_order='x*L+y cell, spatial orbital, up/down',
        no_refit=True)
    model_id = uuid.uuid4().hex+uuid.uuid4().hex
    output.parent.mkdir(parents=True, exist_ok=True)
    stage = output.with_name(output.name+'.partial')
    stage.mkdir(exist_ok=False)
    for name in DATA:
        (stage/name).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source/name, stage/name)
    (stage/'runtime').mkdir()
    for name,path in runtime_files.items():
        shutil.copy2(path, stage/'runtime'/name)
    with (stage/'native_geometry.bin').open('xb') as f:
        f.write(b'LANIO2_VFUNC_V2 ')
        f.write(struct.pack('<8i', 1, ncoef, len(nodes), *grid, 0, 0))
        f.write(np.asarray(half, '<i4').tobytes())
        f.write(np.asarray(primitive, '<f8').tobytes(order='F'))
        f.write(np.asarray(centers[0]@primitive, '<f8').tobytes())
        f.write(np.asarray(spacing, '<f8').tobytes())
        f.write(struct.pack('<2d', kernel['cutoff_A'], kernel['scalar_bias_coefficient_per_Ni']))
        f.write(model_id.encode('ascii'))
        f.write(np.asarray(nodes, '<i4').tobytes(order='C'))
    contract.update(model_id=model_id, source_bundle=str(source), source_code=str(code),
                    native_geometry_bytes=(stage/'native_geometry.bin').stat().st_size)
    for name, expected in identity.items():
        if (stage/name).stat().st_size != expected['bytes']:
            raise RuntimeError(f'Export readback changed {name}')
    save(stage/'V_FUNC.json', contract)
    stage.rename(output)
    return dict(status='EXPORTED', path=str(output), model_id=model_id,
                lossless_readback=True)


def metadata(directory: Path) -> dict:
    meta = json.loads((directory/'V_FUNC.json').read_text())
    if meta['schema'] != SCHEMA:
        raise ValueError('Not a native V1 V_func_tmp package')
    if len(meta['model_id'])!=64:raise ValueError('Model ID must contain 64 characters')
    runtime_path(directory, meta)
    return meta


def verify(directory: Path) -> dict:
    meta = metadata(directory)
    for name, expected in meta['files'].items():
        path=directory/name
        if not path.is_file() or (isinstance(expected,dict) and path.stat().st_size!=expected['bytes']):
            raise ValueError(f'Missing or truncated exported model file: {name}')
    geometry=directory/'native_geometry.bin'
    if not geometry.is_file() or geometry.stat().st_size!=meta.get('native_geometry_bytes',geometry.stat().st_size):
        raise ValueError('Missing or truncated native geometry')
    return meta


def prepare_projection(args) -> dict:
    if args.output.exists():
        raise FileExistsError(args.output)
    meta = metadata(args.model)
    bg = json.loads(Path(str(args.background)+'.json').read_text())
    if bg.get('kind') != 'pure_H0_no_local_fields_no_HU':
        raise ValueError('Verified bare-H0 background required')
    lines = args.translations.read_text().splitlines()
    rows = np.array([[float(v) for v in row.split()] for row in lines[1:] if row.strip()])
    if (rows.ndim != 2 or rows.shape[1] != 4 or int(lines[0]) != len(rows) or not np.isfinite(rows).all()
            or not np.array_equal(rows[:,:3], np.rint(rows[:,:3])) or np.any(rows[:,3] <= 0)
            or len(np.unique(rows[:,:3],axis=0)) != len(rows)):
        raise ValueError('Explicit unique x y z weight rows are required')
    kz = np.array(args.kz)
    if not np.isfinite(kz).all() or len(np.unique(np.round(kz % 1,12))) != len(kz):
        raise ValueError('Distinct finite periodic kz points required')
    if args.grid_weight <= 0 or not np.isfinite(args.grid_weight) or not 1 <= args.nwan <= 30:
        raise ValueError('Invalid integration weight or nwan')
    files = {}
    for i in range(1,args.nwan+1):
        path = args.wan_dir/f'wan_{i}'
        if path.stat().st_size != 16*len(rows)*int(np.prod(meta['primitive_grid'])):
            raise ValueError(f'Wrong complex128 primitive-grid byte count: {path}')
        files[path.name] = dict(bytes=path.stat().st_size)
    contract = dict(model_id=meta['model_id'], Wannier_files=files,
                    translations_bytes=args.translations.stat().st_size, kz=kz.tolist(),
                    background_id=bg.get('background_id',str(args.background.resolve())), grid_weight=args.grid_weight,
                    nwan=args.nwan, source_background=dict(kind=bg['kind'],source=bg.get('source','')))
    basis_id = uuid.uuid4().hex+uuid.uuid4().hex
    args.output.mkdir(parents=True,exist_ok=False)
    np.savetxt(args.output/'kz.dat',kz,header=str(len(kz)),comments='',fmt='%.17g')
    settings=dict(model_path=args.model,wan_dir=args.wan_dir,translations_path=args.translations,
                  rvr_dir=args.output,kz_path=args.output/'kz.dat')
    text=['&local_projection']
    text += [f" {key}='{str(path.resolve()).replace(chr(39),chr(39)*2)}'," for key,path in settings.items()]
    text += [f" basis_id='{basis_id}', nwan={args.nwan}, grid_weight={args.grid_weight:.17g}",'/']
    (args.output/'projection.nml').write_text('\n'.join(text)+'\n')
    center=np.array(meta['Ni_fractional'])*np.diag(meta['primitive_cell_A'])
    cutoff=json.loads((args.model/'KERNELS.json').read_text())['cutoff_A']
    lattice=np.diag(meta['primitive_cell_A'])
    domain_lo=np.floor((center[:2]-cutoff)/lattice[:2]).astype(int)
    domain_hi=np.floor((center[:2]+cutoff)/lattice[:2]).astype(int)
    extent=domain_hi-domain_lo+np.ptp(rows[:,:2],axis=0).astype(int)+1
    nr=int(np.prod(extent)); ndiff=int(np.prod(2*np.ptp(rows[:,:2],axis=0).astype(int)+1))
    nc=2*meta['ncoeff']+1
    contract.update(basis_id=basis_id,status='PROJECTION_INPUT_READY_NOT_RUN',
                    native_projection=str((args.output/'projection.nml').resolve()),
                    projected_center_translations=nr,
                    projected_payload_bytes=16*len(kz)*args.nwan**2*(ndiff+nr*nr*nc))
    save(args.output/'PROJECTION.json',contract)
    return dict(status=contract['status'],projection_nml=contract['native_projection'],
                model_id=meta['model_id'],basis_id=basis_id,payload_bytes=contract['projected_payload_bytes'])


def export_background(source: Path, output: Path) -> dict:
    """Use the established DFT-only primitive.bin; exclude its local scalar and N0."""
    if output.exists():
        raise FileExistsError(output)
    def array(f, dtype, shape):
        dtype = np.dtype(dtype); count = int(np.prod(shape))
        raw = f.read(count*dtype.itemsize)
        if len(raw) != count*dtype.itemsize:
            raise ValueError('Truncated DFT-only primitive input')
        a = np.frombuffer(raw, dtype=dtype).reshape(shape, order='F')
        if not np.isfinite(a).all():
            raise ValueError('Nonfinite DFT-only primitive input')
        return a
    with source.open('rb') as f:
        if f.read(17) != b'LANIO2_DFT_V1    ' or struct.unpack('<4i', f.read(16)) != (57600,512,128,30):
            raise ValueError('Expected the registered DFT-only primitive.bin format')
        q = array(f, '<f8', (3,512)).T
        hq = array(f, '<c16', (30,30,512)).transpose(2,0,1)
        array(f, '<c16', (10,10)); array(f, '<f8', (57600,))
        pq = []
        for _ in range(512):
            nb = struct.unpack('<i', f.read(4))[0]
            if not 30 <= nb <= 1000:
                raise ValueError('Invalid outer-window count')
            p = array(f, '<c16', (5,nb)); dis = array(f, '<c16', (nb,30)); rot = array(f, '<c16', (30,30))
            pq.append(p@dis@rot)
        if f.read(1):
            raise ValueError('Unexpected trailing DFT-only data')
    iq = np.rint(q*8).astype(int) % 8
    if len(np.unique(iq, axis=0)) != 512 or np.max(abs(q*8-np.rint(q*8))) > 1e-9:
        raise ValueError('Expected the original full 8x8x8 primitive q mesh')
    if np.max(abs(hq-hq.conj().swapaxes(1,2))) > 1e-9:
        raise ValueError('Non-Hermitian bare H0 input')
    rr = np.array(list(product(range(-4,5), repeat=3)), dtype=np.int32)
    degeneracy = 2**np.sum(abs(rr)==4, axis=1)
    phase = np.exp(-2j*np.pi*(rr@q.T))/(512*degeneracy[:,None])
    hr = (phase@hq.reshape(512,-1)).reshape(729,30,30)
    pr = (phase@np.array(pq).reshape(512,-1)).reshape(729,5,30)
    back = np.einsum('qr,rij->qij', np.exp(2j*np.pi*(q@rr.T)), hr, optimize=True)
    if np.max(abs(back-hq)) > 1e-9:
        raise ValueError('Bare H0 transform roundtrip failed')
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('xb') as f:
        f.write(b'LANIO2_AU_V1    ')
        f.write(struct.pack('<8i',16909060,729,30,5,10,810,1,0))
        f.write(np.asarray(rr.T,'<i4').tobytes(order='F'))
        for a in (hr.transpose(1,2,0), pr.transpose(1,2,0), np.zeros((810,810),complex), np.zeros((10,10),complex)):
            f.write(np.asarray(a,'<c16').tobytes(order='F'))
    report = dict(kind='pure_H0_no_local_fields_no_HU', source=str(source.resolve()),
        background_id=uuid.uuid4().hex+uuid.uuid4().hex,
        transformation='original weighted WS inverse transform; p(q)=raw_projector*dis*U',
        H0_source_roundtrip_max_eV=float(np.max(abs(back-hq))),
        excluded='source scalar, source N0 and all old fitted potential kernels')
    save(Path(str(output)+'.json'), report)
    return report


def read_structure(path: Path, meta: dict, features):
    with np.load(path, allow_pickle=False) as z:
        s = features.Structure(*(z[k].copy() for k in ('cell','positions','species','ni','spins')))
    s.validate()
    cell0 = np.array(meta['primitive_cell_A']); inv = np.linalg.inv(cell0)
    scales = s.cell@inv; length = int(round(scales[0, 0]))
    if length <= 0 or not np.allclose(scales, np.diag([length,length,1]), atol=1e-9):
        raise ValueError('RVR belongs to this fixed geometry and LxLx1 cells only')
    frac = s.positions @ inv
    found = sorted((int(z), *np.round(f % [length,length,1], 8)) for z, f in zip(s.species, frac))
    expected = sorted((int(row[0]), *np.round(np.array(row[1:])+[x,y,0], 8))
                      for x in range(length) for y in range(length) for row in meta['primitive_motif'])
    if found != expected:
        raise ValueError('Atomic structure differs from the exported RVR geometry')
    delta = frac[s.ni]-meta['Ni_fractional']
    if not np.allclose(delta, np.rint(delta), atol=1e-9):
        raise ValueError('Ni positions disagree with native center convention')
    xy = np.rint(delta[:, :2]).astype(int) % length
    order = np.argsort(xy[:, 0]*length+xy[:, 1])
    return replace(s, ni=s.ni[order], spins=s.spins[order]), length, order


def write_coeff(path: Path, meta: dict, length: int, p: np.ndarray, n: np.ndarray,
                derivative: bool) -> None:
    if p.shape != (length*length, meta['ncoeff'], 4) or n.shape != (length*length, 10, 10):
        raise ValueError('Native coefficient shape mismatch')
    if not np.isfinite(p).all() or not np.isfinite(n).all() or np.max(abs(n-n.conj().swapaxes(1, 2))) > 1e-10:
        raise ValueError('Nonfinite prediction or non-Hermitian occupation/derivative')
    with path.open('xb') as f:
        f.write(b'LANIO2_COEFF_V2 ')
        f.write(struct.pack('<8i', 1, meta['ncoeff'], length, length*length, int(derivative), 0, 0, 0))
        f.write(meta['model_id'].encode('ascii'))
        f.write(struct.pack('<2d', meta['Ueff_eV'], 0. if derivative else meta['scalar_offset_eV']))
        f.write(np.asarray(p.transpose(1, 2, 0), '<f8').tobytes(order='F'))
        f.write(np.asarray(n.transpose(1, 2, 0), '<c16').tobytes(order='F'))


def evaluate(directory: Path, structure: Path, output: Path,
             derivative_sites: list[int], step: float, *, keep_intermediates: bool=False,
             tangent: bool=False) -> dict:
    meta = metadata(directory)
    if meta['schema'] != SCHEMA or step <= 0 or not np.isfinite(step):
        raise ValueError('Invalid package or derivative step')
    predict, features = runtime(directory, meta)
    bundle = predict.Bundle(directory)
    s, length, order = read_structure(structure, meta, features)
    if len(set(derivative_sites)) != len(derivative_sites) or any(i < 1 or i > len(s.ni) for i in derivative_sites):
        raise ValueError('Unique 1-based native Ni derivative sites required')
    with np.load(structure, allow_pickle=False) as data:
        if tangent:
            if not all(key in data for key in CHART_KEYS):
                raise ValueError('Tangent derivatives require the fixed background chart')
            base, sites, basis, eta = (data[key].copy() for key in CHART_KEYS)
            eta = np.asarray(eta, dtype=float)
            if not np.array_equal(order, np.arange(len(order))) or sites.tolist() != derivative_sites:
                raise ValueError('Chart and derivative sites must use the same native ordering')
            if not np.allclose(chart_spins(base, sites, basis, eta), s.spins, rtol=0, atol=1e-12):
                raise ValueError('Structure spins do not match the fixed chart')
        elif any(key in data for key in CHART_KEYS):
            raise ValueError('Chart structures require --tangent; global-y derivatives are incompatible')
    output.mkdir(parents=True, exist_ok=False)
    def value(spins):
        pred = bundle.model.predict(replace(s, spins=spins))
        return bundle.model.physical_p(pred), pred.N
    p, n = value(s.spins)
    write_coeff(output/'value.bin', meta, length, p, n, False)
    derivatives_p, derivatives_n = [], []
    modes = [(site, alpha) for site in derivative_sites for alpha in (1, 2)] if tangent else [(site, 0) for site in derivative_sites]
    for (site, alpha), filename in zip(modes, derivative_files(derivative_sites, tangent)):
        samples = []
        for angle in (step, -step, 2*step, -2*step):
            if tangent:
                displaced = eta.copy()
                displaced[derivative_sites.index(site), alpha-1] += angle
                spins = chart_spins(base, sites, basis, displaced)
            else:
                spins = s.spins.copy(); v = spins[site-1].copy()
                spins[site-1] = [np.cos(angle)*v[0]+np.sin(angle)*v[2], v[1],
                                 -np.sin(angle)*v[0]+np.cos(angle)*v[2]]
            samples.append(value(spins))
        dp, dn = ((8*(samples[0][i]-samples[1][i])-samples[2][i]+samples[3][i])/(12*step)
                  for i in (0, 1))
        if keep_intermediates:
            derivatives_p.append(dp); derivatives_n.append(dn)
        write_coeff(output/filename, meta, length, dp, dn, True)
    if keep_intermediates:
        np.savez(output/'prediction.npz', p=p, N=n, VU=meta['Ueff_eV']*(.5*np.eye(10)-n),
                 dP=np.array(derivatives_p), dN=np.array(derivatives_n),
                 derivative_sites=derivative_sites, native_from_input_order=order)
    report = dict(status='COMPLETE', model_id=meta['model_id'], L=length,
                  derivative_sites=derivative_sites, derivative_step_rad=step,
                  prediction_rule=PREDICTION_RULE,
                  native_from_input_order=order.tolist(), scalar_offset_eV=meta['scalar_offset_eV'])
    if tangent:
        report.update(coordinates='fixed_exponential_chart', derivative_modes=modes,
                      chart={key: value.tolist() for key, value in zip(CHART_KEYS, (base, sites, basis, eta))})
    else:
        report['derivative_axis'] = [0,1,0]
    save(output/'PREDICTION.json', report)
    return report


def hamiltonian(args) -> dict:
    background = json.loads(Path(str(args.background)+'.json').read_text())
    if (background.get('kind') != 'pure_H0_no_local_fields_no_HU'):
        raise ValueError('Require a verified pure-H0/projector background, never an old fitted scalar background')
    projection = json.loads((args.rvr.parent/'PROJECTION.json').read_text())
    with args.rvr.open('rb') as f:
        header=f.read(176)
    if (len(header)!=176 or header[:16]!=b'LANIO2_LRVR_V2  '
            or header[48:112].decode()!=metadata(args.model)['model_id']
            or header[112:176].decode()!=projection['basis_id']
            or projection.get('background_id',str(args.background.resolve()))!=background.get('background_id',str(args.background.resolve()))):
        raise ValueError('B_general model/basis/background binding mismatch')
    work = args.output.absolute()
    report = evaluate(args.model, args.structure, work, args.derivative_sites, args.step,
                      keep_intermediates=args.keep_intermediates, tangent=getattr(args, 'tangent', False))
    indices = {}
    intermediates = [work/'PREDICTION.json']
    dimension = 60*report['L']**2
    for coeff in [work/'value.bin', *sorted(work.glob('dsite_*.bin'))]:
        prefix = work/coeff.stem
        control = prefix.with_suffix('.nml')
        values = dict(rvr_path=args.rvr, background_path=args.background,
                      coefficient_path=coeff, kpoints_path=args.kpoints, output_prefix=prefix)
        lines = ['&local_hamiltonian']
        for key, value in values.items():
            name = str(Path(value).resolve()).replace("'", "''")
            lines.append(f" {key}='{name}',")
        lines += [" background_kind='pure_H0_no_local_fields_no_HU',",
                  f" dump_components={'.true.' if args.components else '.false.'}", '/']
        control.write_text('\n'.join(lines)+'\n')
        subprocess.run([str(args.binary.resolve()), str(control)], check=True)
        index = prefix.with_suffix('.index')
        rows = index.read_text().splitlines()
        if len(rows)<5:
            raise RuntimeError('Incomplete native Hamiltonian index; intermediates retained')
        length, nd, nk = map(int, rows[1].split()[:3])
        if (length!=report['L'] or nd!=dimension or nk!=len(rows)-4
                or [s.strip() for s in rows[2:4]]!=[report['model_id'],projection['basis_id']]):
            raise RuntimeError('Invalid native Hamiltonian index; intermediates retained')
        for ik in range(1,nk+1):
            if int(rows[ik+3].split()[0])!=ik:
                raise RuntimeError('Incomplete native k-point order; intermediates retained')
            for suffix in (('h','h0','hv','hb','hu') if args.components else ('h',)):
                if (work/f'{coeff.stem}.k{ik:06}.{suffix}').stat().st_size!=16*dimension**2:
                    raise RuntimeError('Truncated native Hamiltonian output; intermediates retained')
        indices[coeff.stem] = rows
        intermediates.extend([coeff,control,index])
    report['native_Hamiltonian_status'] = 'COMPLETE'
    report['intermediates_kept'] = args.keep_intermediates
    save(work/'HAMILTONIAN.json', {**report, 'native_indices':indices,
        'inputs':{key:str(getattr(args,key).resolve())
                  for key in ('model','structure','rvr','background','kpoints','binary')}})
    if not args.keep_intermediates:
        for path in intermediates:
            path.unlink()
    return report


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__); sub = p.add_subparsers(dest='command', required=True)
    e = sub.add_parser('export')
    for name in ('source','code','structure','output'):
        e.add_argument('--'+name, type=Path, required=True)
    e.add_argument('--scalar-offset-eV', type=float, required=True)
    v = sub.add_parser('verify'); v.add_argument('model', type=Path)
    b = sub.add_parser('export-background')
    b.add_argument('--source', type=Path, required=True); b.add_argument('--output', type=Path, required=True)
    b = sub.add_parser('prepare-projection')
    for key in ('model','wan-dir','translations','background','output'):
        b.add_argument('--'+key,type=Path,required=True)
    b.add_argument('--grid-weight',type=float,required=True);b.add_argument('--nwan',type=int,default=30)
    b.add_argument('--kz',type=float,nargs='+',required=True)
    for name in ('evaluate','cal_Hsk'):
        a = sub.add_parser(name)
        for key in ('model','structure','output'):
            a.add_argument('--'+key, type=Path, required=True)
        a.add_argument('--derivative-sites', type=int, nargs='*', default=[])
        a.add_argument('--tangent', action='store_true', help='differentiate the fixed chart stored in structure.npz')
        a.add_argument('--step', type=float, default=.001)
        a.add_argument('--keep-intermediates', action='store_true',
                       help='keep diagnostic prediction arrays and, for cal_Hsk, native inputs and indices')
        if name == 'cal_Hsk':
            for key in ('rvr','background','kpoints','binary'):
                a.add_argument('--'+key, type=Path, required=True)
            a.add_argument('--components', action='store_true')
    args = p.parse_args()
    if args.command == 'export':
        result = export(args.source, args.code, args.structure, args.output, args.scalar_offset_eV)
    elif args.command == 'verify':
        result = dict(status='VERIFIED', model_id=verify(args.model)['model_id'])
    elif args.command == 'export-background':
        result = export_background(args.source, args.output)
    elif args.command == 'prepare-projection':
        result = prepare_projection(args)
    elif args.command == 'evaluate':
        result = evaluate(args.model, args.structure, args.output, args.derivative_sites, args.step,
                          keep_intermediates=args.keep_intermediates, tangent=args.tangent)
    else:
        result = hamiltonian(args)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
