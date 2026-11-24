#!/usr/bin/env python3
"""Generate 26 magnetic SCF inputs using the V1 sampling rules; run no jobs."""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
import re

import numpy as np

BOHR_ANGSTROM = 0.529177210903
NI_POSITIONS = np.array([[.25, .25, .5], [.75, .25, .5], [.25, .75, .5], [.75, .75, .5]])
KEY = re.compile(r"([a-zA-Z][\w]*(?:\(\s*\d+(?:\s*,\s*\d+)*\s*\))?)\s*=")
MAGNETIC = re.compile(r"(?:starting_magnetization|angle1|angle2|constrain_species)\(\d+\)$")
CARDS = {'ATOMIC_SPECIES', 'ATOMIC_POSITIONS', 'CELL_PARAMETERS', 'K_POINTS', 'HUBBARD'}


def configurations(seed: int = 0) -> dict[str, dict]:
    """Target polar/azimuthal angles in degrees, in the registered Ni order."""
    if type(seed) is not int or seed < 0:
        raise ValueError('The random seed must be a nonnegative integer')
    rng = np.random.Generator(np.random.PCG64(seed))
    result = {'fm': dict(angles=[[0., 0.]] * 4, constrained=False)}
    for name, phi in (('ni56_pp', 0.), ('ni56_pm', 180.)):
        result[name] = dict(angles=[[5., 0.], [5., phi], [0., 0.], [0., 0.]], constrained=True)
    for i in range(1, 10):
        angles = np.column_stack((180 * rng.random(4), 360 * rng.random(4)))
        result[f'cri3_full_{i:02}'] = dict(angles=angles.tolist(), constrained=True)
    for i in range(1, 10):
        angles = np.column_stack((rng.integers(0, 90, 4), rng.integers(0, 360, 4)))
        result[f'seki_upper_{i:02}'] = dict(angles=angles.tolist(), constrained=True)
    for angle in (2, 6):
        for sign, phi in (('p', (0., 180.)), ('m', (180., 0.))):
            result[f'rel78_{sign}{angle}'] = dict(angles=[[0., 0.], [0., 0.], [angle, phi[0]], [angle, phi[1]]], constrained=True)
    result['Neel_anchor'] = dict(angles=[[0., 0.], [180., 0.], [180., 0.], [0., 0.]], constrained=True)
    return result


def mask_strings(text: str) -> str:
    masked = list(text)
    quote = None
    i = 0
    while i < len(text):
        char = text[i]
        if quote:
            masked[i] = ' '
            if char == quote:
                if i + 1 < len(text) and text[i + 1] == quote:
                    masked[i + 1] = ' '
                    i += 1
                else:
                    quote = None
        elif char in "'\"":
            quote = char
            masked[i] = ' '
        elif char == '!':
            end = text.find('\n', i)
            end = len(text) if end < 0 else end
            masked[i:end] = ' ' * (end - i)
            i = end - 1
        i += 1
    if quote:
        raise ValueError('Unterminated quoted string in SCF input')
    return ''.join(masked)


def without_comments(text: str) -> str:
    mask_strings(text)
    lines = []
    for raw in text.splitlines():
        # A comment starts at the first unquoted !, which was masked above.
        quote = None
        end = len(raw)
        i = 0
        while i < len(raw):
            c = raw[i]
            if quote:
                if c == quote:
                    if i + 1 < len(raw) and raw[i + 1] == quote:
                        i += 1
                    else:
                        quote = None
            elif c in "'\"":
                quote = c
            elif c == '!':
                end = i
                break
            i += 1
        lines.append(raw[:end])
    return '\n'.join(lines)


def number(value: str) -> float:
    result = float(value.lower().replace('d', 'e'))
    if not np.isfinite(result):
        raise ValueError('Finite numerical input required')
    return result


def integer(value: str) -> int:
    result = number(value)
    if not result.is_integer():
        raise ValueError(f'Integer required: {value}')
    return int(result)


def unquote(value: str) -> str:
    if len(value) < 2 or value[0] not in "'\"" or value[-1] != value[0]:
        raise ValueError(f'Quoted string required: {value}')
    return value[1:-1].replace(value[0] * 2, value[0])


def quoted(value: str) -> str:
    if '\n' in value or '\r' in value:
        raise ValueError('Newlines are not allowed in QE paths')
    return "'" + value.replace("'", "''") + "'"


def read_scf(text: str) -> dict:
    text = without_comments(text)
    masked = mask_strings(text)
    groups = {}
    cursor = 0
    while True:
        match = re.match(r'\s*&([A-Za-z]+)\b', masked[cursor:])
        if match is None:
            break
        start = cursor + match.end()
        end = masked.find('/', start)
        if end < 0:
            raise ValueError('Unterminated QE namelist')
        name = match[1].upper()
        if name not in ('CONTROL', 'SYSTEM', 'ELECTRONS') or name in groups:
            raise ValueError(f'Unsupported or duplicate namelist: {name}')
        body, plain = text[start:end], masked[start:end]
        assignments = list(KEY.finditer(plain))
        if not assignments or plain[:assignments[0].start()].strip(' ,\n\t'):
            raise ValueError(f'Invalid namelist assignments: {name}')
        values = {}
        for i, item in enumerate(assignments):
            key = re.sub(r'\s+', '', item[1]).lower()
            stop = assignments[i + 1].start() if i + 1 < len(assignments) else len(body)
            value = body[item.end():stop].strip().rstrip(',').strip()
            if key in values or not value or '\n' in value:
                raise ValueError(f'Duplicate or empty setting: {key}')
            # Only scalar assignments are supported; arrays must use indexed keys.
            if ',' in mask_strings(value) or '=' in mask_strings(value):
                raise ValueError(f'Unsupported value for {key}: {value}')
            if not re.fullmatch(r"'(?:[^']|'')*'|\"(?:[^\"]|\"\")*\"|\.(?:true|false)\.|[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?", value, re.I):
                raise ValueError(f'Unsupported scalar assignment: {key}={value}')
            values[key] = value
        groups[name] = values
        cursor = end + 1
    if set(groups) != {'CONTROL', 'SYSTEM', 'ELECTRONS'}:
        raise ValueError('Require CONTROL, SYSTEM and ELECTRONS namelists')
    if unquote(groups['CONTROL'].get('calculation', "'scf'")).lower() != 'scf':
        raise ValueError('The base input must be an SCF calculation')
    system = groups['SYSTEM']
    nat, ntyp = integer(system['nat']), integer(system['ntyp'])
    if nat not in (4, 16) or not 3 <= ntyp <= 6:
        raise ValueError('Require a 4- or 16-atom LaNiO2 input with 3 to 6 species')
    lines = [line.strip() for line in text[cursor:].splitlines() if line.strip()]
    cards = {}
    i = 0
    while i < len(lines):
        heading = lines[i].split()[0].upper()
        if heading not in CARDS or heading in cards:
            raise ValueError(f'Unsupported or duplicate QE card: {lines[i]}')
        option = lines[i][len(heading):].strip().strip('(){}').lower()
        i += 1
        if heading in ('ATOMIC_SPECIES', 'ATOMIC_POSITIONS', 'CELL_PARAMETERS'):
            count = {'ATOMIC_SPECIES': ntyp, 'ATOMIC_POSITIONS': nat, 'CELL_PARAMETERS': 3}[heading]
        elif heading == 'K_POINTS':
            if option == 'automatic':
                count = 1
            elif option == 'crystal':
                if i >= len(lines) or integer(lines[i]) < 1:
                    raise ValueError('A positive crystal k-point count is required')
                count = 1 + integer(lines[i])
            else:
                raise ValueError('Only automatic or crystal K_POINTS are supported')
        else:
            count = 0
            while i + count < len(lines) and lines[i + count].split()[0].upper() not in CARDS:
                count += 1
        data = lines[i:i + count]
        if len(data) != count or not data:
            raise ValueError(f'Incomplete card: {heading}')
        cards[heading] = (option, data)
        i += count
    if not {'ATOMIC_SPECIES', 'ATOMIC_POSITIONS', 'K_POINTS'} <= set(cards):
        raise ValueError('Require species, positions and k points')
    return {'groups': groups, 'cards': cards}


def prepare(base: dict, pseudo_dir: Path) -> dict:
    base = copy.deepcopy(base)
    groups, cards = base['groups'], base['cards']
    system = groups['SYSTEM']
    if system.get('lspinorb', '.false.').lower() not in ('.false.', 'false', 'f'):
        raise ValueError('The frozen V1 input requires no spin-orbit coupling')
    if 'tot_magnetization' in system:
        raise ValueError('tot_magnetization is incompatible with the generated noncollinear inputs')
    for key in system:
        if '(' in key and not MAGNETIC.fullmatch(key) and not key.startswith('celldm('):
            raise ValueError(f'Cannot safely remap indexed SYSTEM setting: {key}')
    species = {}
    for line in cards['ATOMIC_SPECIES'][1]:
        tokens = line.split()
        if len(tokens) != 3 or not re.fullmatch(r'(?:La|Ni|O)\d*', tokens[0]) or tokens[0] in species:
            raise ValueError(f'Invalid La/Ni/O species: {line}')
        if number(tokens[1]) <= 0 or Path(tokens[2]).name != tokens[2]:
            raise ValueError('Positive species mass and a pseudopotential basename required')
        species[tokens[0]] = (tokens[1], tokens[2])
    if cards['ATOMIC_POSITIONS'][0] != 'crystal':
        raise ValueError('ATOMIC_POSITIONS crystal is required')
    atoms = []
    for line in cards['ATOMIC_POSITIONS'][1]:
        tokens = line.split()
        if len(tokens) != 4 or tokens[0] not in species:
            raise ValueError(f'Unsupported atomic position: {line}')
        xyz = np.array([number(x) for x in tokens[1:]])
        atoms.append((tokens[0], xyz % 1))
    nni = sum(name.startswith('Ni') for name, _ in atoms)
    if nni not in (1, 4) or len(atoms) != 4 * nni:
        raise ValueError('Require a one-Ni primitive cell or the registered four-Ni 2x2 cell')
    if sum(name.startswith('La') for name, _ in atoms) != nni or sum(name.startswith('O') for name, _ in atoms) != 2 * nni:
        raise ValueError('Require LaNiO2 stoichiometry')
    ibrav = integer(system['ibrav'])
    if ibrav == 6:
        if 'CELL_PARAMETERS' in cards:
            raise ValueError('CELL_PARAMETERS conflicts with ibrav=6')
        a = number(system['celldm(1)']) * BOHR_ANGSTROM
        cell = np.diag([a, a, a * number(system['celldm(3)'])])
    elif ibrav == 0 and cards.get('CELL_PARAMETERS', ('', []))[0] == 'angstrom':
        cell = np.array([[number(x) for x in line.split()] for line in cards['CELL_PARAMETERS'][1]])
    else:
        raise ValueError('Require ibrav=6 with celldm, or ibrav=0 with CELL_PARAMETERS angstrom')
    if cell.shape != (3, 3) or not np.allclose(cell, np.diag(np.diag(cell)), atol=1e-10, rtol=0) or np.any(np.diag(cell) <= 0) or not np.isclose(cell[0, 0], cell[1, 1], atol=1e-10, rtol=0):
        raise ValueError('Require the positive tetragonal V1 cell')
    option, points = cards['K_POINTS']
    if option == 'automatic':
        mesh = [integer(x) for x in points[0].split()]
        if len(mesh) != 6 or min(mesh[:3]) <= 0 or any(mesh[3:]):
            raise ValueError('Require a positive, unshifted automatic k mesh')
        if nni == 1:
            if mesh[0] % 2 or mesh[1] % 2:
                raise ValueError('Primitive in-plane k mesh must be divisible by two')
            mesh[0] //= 2
            mesh[1] //= 2
        base['kpoints'] = 'K_POINTS automatic\n' + ' '.join(map(str, mesh)) + '\n'
    else:
        if nni == 1:
            raise ValueError('Use an automatic k mesh when expanding a primitive cell')
        rows = np.array([[number(x) for x in line.split()] for line in points[1:]])
        if rows.shape != (integer(points[0]), 4) or np.any(rows[:, 3] < 0) or not np.isclose(rows[:, 3].sum(), 1., atol=1e-9, rtol=0):
            raise ValueError('Invalid crystal k-point list or weights')
        base['kpoints'] = 'K_POINTS crystal\n' + '\n'.join(points) + '\n'
    hubbard_card = cards.get('HUBBARD', ('ortho-atomic', [f'U {name}-3d 1.0' for name in species if name.startswith('Ni')]))
    if hubbard_card[0] != 'ortho-atomic':
        raise ValueError('Only HUBBARD (ortho-atomic) is supported')
    hubbard = {}
    for line in hubbard_card[1]:
        row = line.split()
        if len(row) != 3 or row[0] != 'U' or not row[1].endswith('-3d') or row[1][:-3] not in species or not row[1].startswith('Ni'):
            raise ValueError(f'Only on-site Ni-3d U is supported: {line}')
        name = row[1][:-3]
        if name in hubbard or abs(number(row[2]) - 1.) > 1e-12:
            raise ValueError('Require one explicit U=1 eV entry per Ni species for frozen V1')
        hubbard[name] = row[2]
    if any(name not in hubbard for name, _ in atoms if name.startswith('Ni')):
        raise ValueError('Every Ni species needs an explicit Hubbard U')
    if nni == 1:
        atoms = [(name, (xyz + [x, y, 0]) / [2, 2, 1])
                 for name, xyz in atoms for y in range(2) for x in range(2)]
        cell[:2] *= 2
        if 'nbnd' in system:
            system['nbnd'] = str(4 * integer(system['nbnd']))
        for key in ('nr1', 'nr2', 'nr1s', 'nr2s'):
            if key in system:
                system[key] = str(2 * integer(system[key]))
        if number(system.get('tot_charge', '0')) != 0:
            raise ValueError('Primitive-cell expansion currently requires a neutral base')
    ni = []
    for target in NI_POSITIONS:
        matches = [(name, xyz) for name, xyz in atoms if name.startswith('Ni') and np.allclose(xyz, target, rtol=0, atol=1e-8)]
        if len(matches) != 1:
            raise ValueError('Ni positions do not match the registered V1 2x2 ordering')
        ni.append(matches[0])
    original_names = list(species)
    magnetizations = [system.get(f'starting_magnetization({original_names.index(name)+1})', '0.05') for name, _ in ni]
    if any(number(value) <= 0 for value in magnetizations):
        raise ValueError('Each Ni needs a positive starting magnetization')
    nonmagnetic = {}
    for element in ('La', 'O'):
        names = {name for name, _ in atoms if name.startswith(element)}
        values = {species[name] for name in names}
        if len(values) != 1:
            raise ValueError(f'Cannot merge distinct {element} species')
        nonmagnetic[element] = values.pop()
    base['species'] = [('La', *nonmagnetic['La'])] + [(f'Ni{i+1}', *species[name]) for i, (name, _) in enumerate(ni)] + [('O', *nonmagnetic['O'])]
    base['atoms'] = [('La', xyz) for name, xyz in atoms if name.startswith('La')] + [(f'Ni{i+1}', xyz) for i, (_, xyz) in enumerate(ni)] + [('O', xyz) for name, xyz in atoms if name.startswith('O')]
    base.update(cell=cell, hubbard=[hubbard[name] for name, _ in ni], magnetizations=magnetizations,
                hubbard_source='base input' if 'HUBBARD' in cards else 'added frozen V1 ortho-atomic U=1 eV')
    for _, _, filename in base['species']:
        if not (pseudo_dir / filename).is_file():
            raise FileNotFoundError(pseudo_dir / filename)
    for key in list(system):
        if MAGNETIC.fullmatch(key) or key.startswith('celldm('):
            del system[key]
    system.update(ibrav='0', nat='16', ntyp='6', noncolin='.true.', lspinorb='.false.', nosym='.true.', noinv='.true.',
                  constrained_magnetization="'atomic vector direction'", **{'lambda': '0.1'})
    system.pop('nspin', None)
    groups['CONTROL'].update(calculation="'scf'", restart_mode="'from_scratch'", outdir="'./tmp'", pseudo_dir=quoted(str(pseudo_dir)))
    groups['CONTROL'].pop('wfcdir', None)
    # A fresh generated case has no saved charge or wavefunctions to restart from.
    groups['ELECTRONS'].update(startingpot="'atomic'", startingwfc="'atomic+random'")
    return base


def render(base: dict, angles: list[list[float]], constrained: bool = True) -> tuple[str, list[list[float]]]:
    angles = np.asarray(angles, dtype=float)
    if angles.shape != (4, 2) or not np.isfinite(angles).all() or np.any(angles[:, 0] < 0) or np.any(angles[:, 0] > 180):
        raise ValueError('Require four finite polar/azimuthal angle pairs in degrees')
    groups = copy.deepcopy(base['groups'])
    system = groups['SYSTEM']
    if not constrained:
        system.pop('constrained_magnetization', None)
        system.pop('lambda', None)
    for i in range(1, 7):
        magnetic = 2 <= i <= 5
        system[f'constrain_species({i})'] = '.true.' if magnetic and constrained else '.false.'
        system[f'starting_magnetization({i})'] = base['magnetizations'][i-2] if magnetic else '0.0'
        if magnetic:
            system[f'angle1({i})'] = f'{angles[i-2, 0]:.12g}'
            system[f'angle2({i})'] = f'{angles[i-2, 1]:.12g}'
    lines = []
    for name, values in groups.items():
        lines += ['&' + name] + [f' {key}={value},' for key, value in values.items()] + ['/']
    lines += ['ATOMIC_SPECIES'] + [' '.join(row) for row in base['species']]
    lines += ['CELL_PARAMETERS angstrom'] + [' '.join(f'{x:.12g}' for x in row) for row in base['cell']]
    lines += ['ATOMIC_POSITIONS crystal'] + [name + ' ' + ' '.join(f'{x:.12g}' for x in xyz) for name, xyz in base['atoms']]
    lines += [base['kpoints'].rstrip(), 'HUBBARD (ortho-atomic)']
    lines += [f'U Ni{i+1}-3d {value}' for i, value in enumerate(base['hubbard'])]
    theta, phi = np.deg2rad(angles).T
    directions = np.column_stack((np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)))
    return '\n'.join(lines) + '\n', directions.tolist()


def generate(scf: Path, output: Path, pseudo_dir: Path | None = None, cases: list[str] | None = None, *, seed: int = 0) -> dict:
    if output.exists() or output.is_symlink():
        raise FileExistsError(f'Refusing existing output: {output}')
    source = scf.read_bytes()
    selected = configurations(seed)
    parsed = read_scf(source.decode())
    if pseudo_dir is None:
        value = parsed['groups']['CONTROL'].get('pseudo_dir')
        if value is None:
            raise ValueError('Provide pseudo_dir in the base SCF input or use --pseudo-dir')
        pseudo_dir = scf.resolve().parent / unquote(value)
    pseudo_dir = pseudo_dir.resolve()
    base = prepare(parsed, pseudo_dir)
    if cases is not None:
        if len(set(cases)) != len(cases) or not cases or any(name not in selected for name in cases):
            raise ValueError('Unknown, duplicate or empty configuration selection')
        selected = {name: selected[name] for name in cases}
    rendered, records = {}, {}
    for name, spec in selected.items():
        if not re.fullmatch(r'[A-Za-z0-9_]+', name):
            raise ValueError('Unsafe configuration name')
        angles = spec['angles']
        text, directions = render(base, angles, spec['constrained'])
        rendered[name] = text
        records[name] = {'angles_degrees': angles, 'target_directions': directions, 'constrained': spec['constrained'],
                         'scf_sha256': hashlib.sha256(text.encode()).hexdigest()}
    record = dict(source_scf=str(scf.resolve()), source_sha256=hashlib.sha256(source).hexdigest(),
                  seed=seed, rng='numpy.Generator(PCG64)', numpy_version=np.__version__,
                  sampling={'cri3_full': '9 states; four independent Ni; theta uniform [0,180), phi uniform [0,360) degrees',
                            'seki_upper': '9 states; four independent Ni; integer theta 0..89, phi 0..359 degrees'},
                  initialization='atomic; no saved density or wavefunctions copied',
                  scratch='each case uses ./tmp for charge and wavefunctions; inherited wfcdir removed',
                  hubbard_source=base['hubbard_source'],
                  Ni_order='Ni1=(1/4,1/4,1/2), Ni2=(3/4,1/4,1/2), Ni3=(1/4,3/4,1/2), Ni4=(3/4,3/4,1/2)',
                  configurations=records)
    output.mkdir(parents=True, exist_ok=False)
    for name, text in rendered.items():
        directory = output/name
        directory.mkdir()
        (directory/'scf.in').write_text(text)
    (output/'CONFIGURATIONS.json').write_text(json.dumps(record, indent=2, allow_nan=False) + '\n')
    return record


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scf', required=True, type=Path, help='preceding primitive or four-Ni SCF input')
    parser.add_argument('--output', required=True, type=Path, help='new directory for generated inputs')
    parser.add_argument('--pseudo-dir', type=Path, help='external pseudopotentials; default: pseudo_dir relative to the base SCF input')
    parser.add_argument('--seed', type=int, default=0, help='nonnegative random seed (default: 0)')
    parser.add_argument('--cases', nargs='+', choices=list(configurations()), help='default: all 26 configurations')
    args = parser.parse_args()
    result = generate(args.scf, args.output, args.pseudo_dir, args.cases, seed=args.seed)
    print(f"Generated {len(result['configurations'])} SCF inputs; no calculations or submissions performed.")
