"""Read converged QE 7.6 LaNiO2 saves without losing complex Hubbard phases."""
from pathlib import Path
import re
import xml.etree.ElementTree as ET

import numpy as np

BOHR_ANGSTROM = 0.529177210903
HARTREE_EV = 27.211386245988
ORBITAL_ORDER = ['dz2', 'dxz', 'dyz', 'dx2-y2', 'dxy']
NUMBER = r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?'
COMPLEX = re.compile(r'\(\s*(' + NUMBER + r')\s*,\s*(' + NUMBER + r')\s*\)')


def _numbers(text: str | None, count: int) -> np.ndarray:
    tokens = (text or '').split()
    if len(tokens) != count or any(re.fullmatch(NUMBER, v) is None for v in tokens):
        raise ValueError(f'Expected {count} finite numeric values in QE XML')
    result = np.array([float(v.replace('D', 'e').replace('d', 'e')) for v in tokens])
    if not np.isfinite(result).all():
        raise ValueError('Nonfinite QE XML values')
    return result


def _text(root: ET.Element, path: str) -> str:
    element = root.find(path)
    if element is None or element.text is None or not element.text.strip():
        raise ValueError(f'Missing QE XML value: {path}')
    return element.text.strip()


def _boolean(root: ET.Element, path: str) -> bool:
    value = _text(root, path).lower()
    if value not in ('true', 'false'):
        raise ValueError(f'Invalid QE XML boolean: {path}')
    return value == 'true'


def _grid(root: ET.Element, name: str) -> np.ndarray:
    element = root.find('output/basis_set/' + name)
    if element is None:
        raise ValueError(f'Missing QE XML grid: {name}')
    try:
        grid = np.array([int(element.attrib[f'nr{i}']) for i in (1, 2, 3)])
    except (KeyError, ValueError) as exc:
        raise ValueError(f'Invalid QE XML grid: {name}') from exc
    if np.any(grid <= 0):
        raise ValueError(f'Nonpositive QE XML grid: {name}')
    return grid


def read_save(save: Path) -> dict:
    """Read atomic-order labels; N uses spin-major indices m + 5*spin.

    QE's native occup.txt stores ns_nc(5,5,4,nat) in Fortran order.
    XML Hubbard_ns_nc stores absolute values and cannot replace this file.
    The QE real d-harmonic order already matches V1; no orbital permutation.
    """
    save = Path(save)
    try:
        root = ET.parse(save / 'data-file-schema.xml').getroot()
    except ET.ParseError as exc:
        raise ValueError('Malformed QE data-file-schema.xml') from exc
    for element in root.iter():
        element.tag = element.tag.split('}')[-1]
    if root.attrib.get('Units') != 'Hartree atomic units':
        raise ValueError('Require QE XML in Hartree atomic units')
    if _text(root, 'input/control_variables/calculation') != 'scf':
        raise ValueError('Require an SCF save')
    if not _boolean(root, 'output/convergence_info/scf_conv/convergence_achieved'):
        raise ValueError('SCF has not converged')
    if (not _boolean(root, 'input/spin/noncolin')
            or _boolean(root, 'input/spin/spinorbit')
            or _boolean(root, 'input/spin/lsda')):
        raise ValueError('Require noncollinear nspin=4 without SOC')
    if any(_text(root, section + '/dft/functional').upper() != 'PBE' for section in ('input', 'output')):
        raise ValueError('Current V1 labels require the PBE functional')

    structure = root.find('output/atomic_structure')
    if structure is None or structure.attrib.get('nat') != '16':
        raise ValueError('Require the 16-atom, four-Ni LaNiO2 cell')
    atoms = structure.findall('atomic_positions/atom')
    if len(atoms) != 16:
        raise ValueError('Incomplete QE atomic positions')
    try:
        atoms.sort(key=lambda a: int(a.attrib['index']))
        if [int(a.attrib['index']) for a in atoms] != list(range(1, 17)):
            raise ValueError('Duplicate or missing QE atom indices')
        names = [a.attrib['name'] for a in atoms]
    except (KeyError, ValueError) as exc:
        raise ValueError('Invalid QE atom names or indices') from exc
    elements = [re.fullmatch(r'(La|Ni|O)\d*', name) for name in names]
    if any(match is None for match in elements):
        raise ValueError('Only La, Ni and O species are supported')
    species = np.array([{'La': 57, 'Ni': 28, 'O': 8}[m[1]] for m in elements])
    if [np.count_nonzero(species == z) for z in (57, 28, 8)] != [4, 4, 8]:
        raise ValueError('Require La4Ni4O8 stoichiometry')
    ni = np.flatnonzero(species == 28)
    cell = np.array([_numbers(_text(structure, 'cell/' + a), 3) for a in ('a1', 'a2', 'a3')]) * BOHR_ANGSTROM
    positions = np.array([_numbers(a.text, 3) for a in atoms]) * BOHR_ANGSTROM
    if (not np.allclose(cell, np.diag(np.diag(cell)), atol=1e-8, rtol=0)
            or np.any(np.diag(cell) <= 0)
            or not np.isclose(cell[0, 0], cell[1, 1], atol=1e-8, rtol=0)):
        raise ValueError('Require the positive axis-aligned tetragonal V1 cell')
    fractional_ni = positions[ni] / np.diag(cell)
    registered_ni = np.array([[.25, .25, .5], [.75, .25, .5], [.25, .75, .5], [.75, .75, .5]])
    matches = np.all(np.isclose(fractional_ni[:, None, :], registered_ni[None, :, :], atol=1e-8, rtol=0), axis=2)
    if not np.all(matches.sum(axis=0) == 1) or not np.all(matches.sum(axis=1) == 1):
        raise ValueError('Ni positions do not match the registered V1 2x2x1 layout')

    hubbard = root.find('output/dft/dftU')
    if hubbard is None or _text(hubbard, 'U_projection_type') != 'ortho-atomic':
        raise ValueError('Require ortho-atomic Hubbard projectors')
    if _text(hubbard, 'lda_plus_u_kind') != '0':
        raise ValueError('Require Dudarev Hubbard U without J or intersite V')
    allowed = {'lda_plus_u_kind', 'Hubbard_U', 'U_projection_type', 'Hubbard_Occ', 'Hubbard_ns_nc'}
    if any(element.tag not in allowed for element in hubbard):
        raise ValueError('Unsupported additional Hubbard channels or parameters')
    u_values = {}
    for element in hubbard.findall('Hubbard_U'):
        name = element.attrib.get('specie')
        if element.attrib.get('label') != '3d' or name in u_values:
            raise ValueError('Require one Ni-3d Hubbard U per magnetic species')
        u_values[name] = float(_numbers(element.text, 1)[0] * HARTREE_EV)
    if set(u_values) != {names[i] for i in ni} or any(abs(u - 1.) > 1e-7 for u in u_values.values()):
        raise ValueError('Require U=1 eV on every Ni-3d species and no other Hubbard shells')

    text = (save / 'occup.txt').read_text()
    pairs = COMPLEX.findall(text)
    if len(pairs) != 5 * 5 * 4 * 16 or COMPLEX.sub('', text).strip():
        raise ValueError('occup.txt must contain exactly 1600 complex values for (5,5,4,16)')
    values = np.array([complex(float(a.replace('D', 'e').replace('d', 'e')),
                               float(b.replace('D', 'e').replace('d', 'e'))) for a, b in pairs])
    if not np.isfinite(values).all():
        raise ValueError('Nonfinite Hubbard occupations')
    raw = values.reshape((5, 5, 4, 16), order='F')
    if np.max(np.abs(raw[:, :, :, species != 28])) > 1e-8:
        raise ValueError('Unexpected Hubbard occupations on non-Ni atoms')
    n = np.empty((4, 10, 10), dtype=complex)
    for j, atom in enumerate(ni):
        for s1 in range(2):
            for s2 in range(2):
                # QE stores conj(p_a)*p_b, whereas V1 uses p_a*conj(p_b).
                n[j, 5*s1:5*(s1+1), 5*s2:5*(s2+1)] = raw[:, :, 2*s2+s1, atom].T
    hermitian_error = float(np.max(np.abs(n - n.conj().swapaxes(-1, -2))))
    if hermitian_error > 1e-8:
        raise ValueError('Non-Hermitian full complex Hubbard occupations')
    eigenvalues = np.linalg.eigvalsh(n)
    traces = np.trace(n, axis1=1, axis2=2)
    if (np.max(np.abs(traces.imag)) > 1e-8 or np.any(traces.real < -1e-7)
            or np.any(traces.real > 10 + 1e-7)
            or eigenvalues.min() < -1e-7 or eigenvalues.max() > 1 + 1e-7):
        raise ValueError('Hubbard occupation trace or eigenvalues outside physical bounds')
    xml_n = {}
    for element in hubbard.findall('Hubbard_ns_nc'):
        try:
            atom = int(element.attrib['index']) - 1
        except (KeyError, ValueError) as exc:
            raise ValueError('Invalid XML Hubbard atom index') from exc
        if (atom not in ni or atom in xml_n or element.attrib.get('specie') != names[atom]
                or element.attrib.get('label') != '3d' or element.attrib.get('rank') != '2'
                or element.attrib.get('dims', '').split() != ['10', '10']
                or element.attrib.get('order') != 'F' or element.attrib.get('spin') != '1'):
            raise ValueError('Unsupported or duplicate XML Hubbard matrix metadata')
        xml_n[atom] = _numbers(element.text, 100).reshape((10, 10), order='F')
    if set(xml_n) != set(ni.tolist()):
        raise ValueError('Missing XML Hubbard matrix for a Ni atom')
    xml_n_error = max(float(np.max(np.abs(xml_n[atom] - np.abs(n[j].T)))) for j, atom in enumerate(ni))
    if xml_n_error > 1e-8:
        raise ValueError('occup.txt disagrees with the saved XML Hubbard matrix magnitudes')
    scf_error = float(_numbers(_text(root, 'output/convergence_info/scf_conv/scf_error'), 1)[0])
    steps = int(_text(root, 'output/convergence_info/scf_conv/n_scf_steps'))
    if scf_error < 0 or steps < 1:
        raise ValueError('Invalid SCF convergence metadata')
    initial_species = {}
    for element in root.findall('input/atomic_species/species'):
        name = element.attrib.get('name')
        if name is None or name in initial_species:
            raise ValueError('Invalid input atomic species')
        initial_species[name] = element
    starting_axis = None
    initial_atoms = root.findall('input/atomic_structure/atomic_positions/atom')
    if len(initial_atoms) != 16:
        raise ValueError('Missing input atomic positions for the GGA spin axis')
    for atom in initial_atoms:
        element = initial_species.get(atom.attrib.get('name'))
        if element is None:
            raise ValueError('Input atom has no matching species')
        magnetization = float(_numbers(_text(element, 'starting_magnetization'), 1)[0])
        if magnetization ** 2 <= 1e-12:
            continue
        theta = float(_numbers(_text(element, 'spin_teta'), 1)[0])
        phi = float(_numbers(_text(element, 'spin_phi'), 1)[0])
        starting_axis = np.sign(magnetization) * np.array([np.sin(theta) * np.cos(phi),
                                                         np.sin(theta) * np.sin(phi), np.cos(theta)])
        break
    if starting_axis is None:
        raise ValueError('No nonzero initial magnetic direction for the GGA spin axis')
    total_energy = float(2 * _numbers(_text(root, 'output/total_energy/etot'), 1)[0])
    return dict(cell=cell, positions=positions, species=species, ni=ni, N=n, nat=16,
                prefix=_text(root, 'input/control_variables/prefix'),
                starting_axis=starting_axis, total_energy_Ry=total_energy,
                grid_dense=_grid(root, 'fft_grid'), grid_smooth=_grid(root, 'fft_smooth'),
                metadata=dict(length_unit='angstrom', occupation_unit='dimensionless',
                              atom_order='QE atom index ascending; ni indices are zero-based',
                              orbital_order=ORBITAL_ORDER.copy(), spin_order='up five orbitals, down five orbitals',
                              occupation_source='occup.txt; Fortran ns_nc(5,5,4,nat)',
                              occupation_conversion='N block(s1,s2)=transpose(ns_nc(:,:,2*s2+s1,atom)); spin indices zero-based',
                              U_eV=1., projector='ortho-atomic', functional='PBE',
                              geometry_scope='axis-aligned tetragonal 2x2x1 cell; four registered Ni fractional sites; lattice lengths not fixed',
                              scf_steps=steps, scf_error_Ha=scf_error, xml_magnitude_error=xml_n_error,
                              hermitian_error=hermitian_error, traces=traces.real.tolist(),
                              eigenvalue_range=[float(eigenvalues.min()), float(eigenvalues.max())]))
