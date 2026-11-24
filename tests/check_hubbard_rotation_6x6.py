"""Check +U spin-rotation covariance for the archived 6x6 LaNiO2 case.

Usage: python tests/check_hubbard_rotation_6x6.py /path/to/LaNiO2/.run
"""
import argparse
import hashlib
import json
import os
import struct
import sys
import tempfile
from pathlib import Path

import numpy as np


def read_n(path: Path) -> np.ndarray:
    with np.load(path, allow_pickle=False) as saved:
        return saved["N"]


def run(runs: Path, scratch: Path) -> dict:
    repo = Path(__file__).resolve().parents[1]
    config = json.loads((runs / "V1_complete_energy_L10_20261120_r01/config.json").read_text())
    original_root = Path(config["work"]).parent
    for name in ("work", "input_root", "interface"):
        config[name] = str(runs / Path(config[name]).relative_to(original_root))
    config["complete_energy"].update(work=str(scratch), shared_work=str(scratch), sizes=[6])
    config_path = scratch / "config.json"
    config_path.write_text(json.dumps(config))
    os.environ["ANSATZ_RUN_CONFIG"] = str(config_path)
    sys.path[:0] = [str(repo / "7_spin"), str(repo)]
    import complete_energy
    import complete_readout
    import local_vfunc
    import run_material

    length = 6
    archived = runs / "V1_complete_energy_remaining_J_20261116_r01"
    case_dir = archived / "coefficients/L06/case_01"
    upstream_structure = Path(config["work"]) / "J/L06/mesh2_h0.001/case_01/structure.npz"
    density = archived / "L06/case_01"
    definition = Path(config["input_root"]) / "inputs/L06/CASE.json"
    base_case = json.loads(definition.read_text())["angle_configurations"][1]
    rotated_case = dict(base_case, global_spin_rotation_deg=[0, 180, 0])
    base = run_material.spin_structure(length, base_case)
    rotated = run_material.spin_structure(length, rotated_case)
    spin_vectors = np.diag([-1.0, 1.0, -1.0])
    spinor = np.kron([[0.0, -1.0], [1.0, 0.0]], np.eye(5))
    with np.load(upstream_structure) as saved:
        for key, val in base.items():
            np.testing.assert_allclose(val, saved[key], rtol=0, atol=1e-14)
    np.testing.assert_allclose(rotated["spins"], base["spins"] @ spin_vectors, rtol=0, atol=1e-14)

    model_dir = Path(config["interface"]) / "V_func_tmp"
    ueff = complete_energy.frozen_metadata()["Ueff_eV"]
    coeff = case_dir / "coeff/value.bin"
    offset = complete_energy.value_offset(coeff, length, ueff)
    with coeff.open("rb") as f:
        f.seek(112)
        header_ueff, header_offset = struct.unpack("<2d", f.read(16))
    assert (ueff, offset) == (header_ueff, header_offset)
    predictions = []
    for label, geometry in [("original", base), ("rotated", rotated)]:
        geometry_file = scratch / f"{label}.npz"
        np.savez(geometry_file, **geometry)
        output = scratch / label
        local_vfunc.evaluate(model_dir, geometry_file, output, [], .001, keep_intermediates=True)
        complete_energy.value_offset(output / "value.bin", length, ueff)
        with np.load(output / "prediction.npz") as saved:
            predictions.append((saved["N"], saved["VU"], saved["p"]))
    n_model, v_model, p_model = predictions[0]
    with np.load(case_dir / "coeff/prediction.npz") as saved:
        archived_model_error = float(np.max(abs(n_model - saved["N"])))
        archived_potential_error = float(np.max(abs(v_model - saved["VU"])))
        np.testing.assert_allclose(saved["VU"], ueff*(.5*np.eye(10)-saved["N"]), rtol=0, atol=1e-12)
    n_model_error = float(np.max(abs(predictions[1][0] - spinor @ n_model @ spinor.T)))
    v_model_error = float(np.max(abs(predictions[1][1] - spinor @ v_model @ spinor.T)))
    scalar_error = float(np.max(abs(predictions[1][2][..., 0] - p_model[..., 0])))
    magnetic_error = float(np.max(abs(predictions[1][2][..., 1:] - p_model[..., 1:] @ spin_vectors)))
    np.testing.assert_allclose(v_model, ueff*(.5*np.eye(10)-n_model), rtol=0, atol=1e-12)
    assert max(n_model_error, v_model_error, scalar_error, magnetic_error) < 1e-10

    active = sum(read_n(density / f"density_z{z}.npz") for z in range(8))
    ref = read_n(archived / "reference/background.npz")
    rotated_active = spinor @ active @ spinor.T
    # Reproject the archived Bloch spinors themselves for the tilted Ni and one full kz slice.
    static = complete_readout.static(Path(config["work"]) / "basis/pure_H0.au")
    mu = json.loads((density / "MU.json").read_text())["mu"]
    shell = np.zeros((10, 10), complex)
    rotated_shell = np.zeros_like(shell)
    for ki, k in complete_energy.points(length, 0):
        projector = complete_readout.operators(length, k, *static)[1][0]
        with np.load(density / f"k_{ki:04}.npz") as saved:
            q = saved["Q"]
            weights = complete_readout.fermi(saved["e"], mu)[0] / 32
        rotated_q = np.empty_like(q)
        rotated_q[::2] = -q[1::2]
        rotated_q[1::2] = q[::2]
        for wave, result in [(q, shell), (rotated_q, rotated_shell)]:
            projected = projector @ wave
            result += (projected * weights) @ projected.conj().T
    sector = read_n(density / "density_z0.npz")[0]
    eigenvector_error = float(np.max(abs(shell - sector)))
    rotated_eigenvector_error = float(np.max(abs(rotated_shell - spinor @ sector @ spinor.T)))
    assert max(eigenvector_error, rotated_eigenvector_error) < 1e-12
    orig_total = complete_energy.total_occupation(active, ref, length, base_case)
    rot_total = complete_energy.total_occupation(rotated_active, ref, length, rotated_case)
    archived_total = read_n(density / "DENSITY.npz")
    total_error = float(np.max(abs(rot_total - spinor @ orig_total @ spinor.T)))
    saved_error = float(np.max(abs(orig_total - archived_total)))
    assert max(saved_error, total_error) < 1e-12
    orig_e, orig_v = complete_readout.hubbard(orig_total, ueff)
    rot_e, rot_v = complete_readout.hubbard(rot_total, ueff)
    vtotal_error = float(np.max(abs(rot_v - spinor @ orig_v @ spinor.T)))
    saved_e = json.loads((density / "ENERGY.json").read_text())["terms"]["EU"]
    assert max(abs(orig_e - saved_e), abs(rot_e - orig_e), vtotal_error) < 1e-10

    tiled = [(x % 2) * 2 + y % 2 for x in range(length) for y in range(length)]
    fixed_e, _ = complete_readout.hubbard(rotated_active + ref[tiled], ueff)
    assert abs(fixed_e - rot_e) > .03
    return {
        "case": str(case_dir), "density": str(density), "reference": str(archived / "reference/background.npz"),
        "coefficient_file": str(coeff), "upstream_structure": str(upstream_structure),
        "occupation_checkpoint": str(density / "DENSITY.npz"), "energy_checkpoint": str(density / "ENERGY.json"),
        "occupation_target": str((density / "DENSITY.npz").resolve()),
        "energy_target": str((density / "ENERGY.json").resolve()),
        "model_id": config["model_id"], "Ueff_model_eV": ueff, "Ueff_archive_header_eV": header_ueff,
        "site_1_spins": [base["spins"][0].tolist(), rotated["spins"][0].tolist()],
        "global_spin_rotation_deg": rotated_case["global_spin_rotation_deg"],
        "active_source": "sum of N from density_z0.npz through density_z7.npz",
        "rotated_active_source": "explicit spinor rotation of archived active occupation",
        "model_source": "fresh one-way predictions for both structures, compared with archived coupled coefficients",
        "model_original_vs_checkpoint_N_max": archived_model_error,
        "model_original_vs_checkpoint_V_max_eV": archived_potential_error,
        "model_rotation_N_max": n_model_error, "model_rotation_V_max_eV": v_model_error,
        "rotation_scalar_coeff_max": scalar_error, "rotation_magnetic_coeff_max": magnetic_error,
        "spinor_reprojection": {
            "site_native_1based": 1, "k_indices": [ki for ki, _ in complete_energy.points(length, 0)],
            "original_max_error_against_density_z0": eigenvector_error,
            "rotated_max_error_against_transformed_density_z0": rotated_eigenvector_error,
        },
        "original_vs_checkpoint_Ntotal_max": saved_error, "rotation_Ntotal_max": total_error,
        "rotation_Vtotal_max_eV": vtotal_error, "checkpoint_EU_eV": saved_e,
        "original_EU_eV": orig_e, "rotated_EU_eV": rot_e, "difference_eV": rot_e - orig_e,
        "fixed_background_control_EU_eV": fixed_e, "fixed_background_change_eV": fixed_e - orig_e,
        "source_sha256": {Path(m.__file__).name: hashlib.sha256(Path(m.__file__).read_bytes()).hexdigest()
                          for m in (local_vfunc, run_material, complete_energy, complete_readout)},
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runs", type=Path, help="archive root containing the V1 material runs")
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="hubbard-rotation-") as tmp:
        print(json.dumps(run(args.runs.resolve(), Path(tmp)), indent=2, allow_nan=False))
