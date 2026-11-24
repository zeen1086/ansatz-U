# Frozen V1 pseudopotentials

These are the La/Ni/O UPF files recorded for all 26 frozen V1 training states.
They are PBE, scalar-relativistic, norm-conserving ONCVPSP pseudopotentials,
generated with ONCVPSP 3.3.0. They contain nonlinear core corrections and no
spin-orbit projectors (`has_so=F`).

| File | Valence electrons | Bytes | SHA256 |
|---|---:|---:|---|
| `La.upf` | 11 | 357433 | `7a5b05be056c96292827043b6c564de17805566e789d043250c57bba1cd10c57` |
| `Ni.upf` | 18 | 257099 | `07ac84022756676cfc9f0a328f68bd6b19c25876998ca0289ffb682092306946` |
| `O.upf` | 6 | 137824 | `16c5996fd611f61de83ea416b5248967460ef668837ee05f5b8f3d3dfac79da2` |

The files were recovered unchanged from the `lanio2_u1_g2.save` bootstrap in
`V1_complete_energy_min_20261116_r01`. Their byte counts and SHA256 values match
that archive's `inputs/DENSITY_MANIFEST.json`. An independent cross-check against
`energy_mapping_20261107_r01/RAW_ENERGIES.json` matches every species entry in all
26 states selected by `local_grid_cycle30_final_20261113_r01/prepared/FROZEN.json`.
These are archived provenance checks, not a new DFT calculation.

Use this directory as QE's `pseudo_dir`; use `La.upf`, `Ni.upf`, and `O.upf` in
`ATOMIC_SPECIES`. All separately labelled Ni species use the same `Ni.upf`.
For an external `parameter` file, set `pseudo_dir` to this directory's absolute
path, because preparation paths are resolved relative to the parameter file.

The original generation inputs and GNU GPL notice remain intact in each UPF's
`PP_INFO` section. That notice requests citation of D. R. Hamann, Phys. Rev. B
88, 085117 (2013).
