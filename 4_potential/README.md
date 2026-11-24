# Run SCF preparation, label extraction and fitting

Use an existing Python 3.11+ environment with NumPy/SciPy and the configured QE 7.6 build used for SCF. Keep inputs, pseudopotentials, builds and results outside the repository; use new output directories.

```sh
PY=/absolute/existing/python
REPO=/absolute/ansatz-U
WORK=/absolute/external/work
QE=/absolute/qe-7.6
BASE_SCF=/absolute/base/scf.in
PSEUDO=/absolute/pseudopotentials
PREFIX=your_scf_prefix
mkdir -p "$WORK"
```

Set `PREFIX` to the prefix in `BASE_SCF`. The supported labels are the registered four-Ni, noncollinear, no-SOC PBE LaNiO₂ cell with ortho-atomic U=1 eV.

## 1. Generate and run magnetic SCFs

```sh
"$PY" "$REPO/4_potential/generate_configs.py" \
  --scf "$BASE_SCF" --pseudo-dir "$PSEUDO" --seed 0 \
  --output "$WORK/cases"
```

Run each generated `scf.in` from its own case directory using your usual QE job/launcher. For one single-process case:

```sh
(cd "$WORK/cases/fm" && "$QE/bin/pw.x" -in scf.in > scf.out)
```

Repeat for the other cases. Continue only with successfully converged SCFs and preserve each complete `scf.out` and `tmp/$PREFIX.save` directory.

## 2. Build the label exporter

```sh
make -C "$QE/PP/src" -f "$REPO/4_potential/Makefile.labels" \
  BUILD="$WORK/label-build"
```

The QE libraries must already be built with the project's magnetic modifications. Extraction runs on one MPI rank and does not need wavefunctions.

## 3. Extract each converged case

The save must contain `data-file-schema.xml`, `charge-density.dat` and `occup.txt`; the original complete SCF log is also required.

```sh
CASE=fm
"$PY" "$REPO/4_potential/extract_labels.py" \
  --save "$WORK/cases/$CASE/tmp/$PREFIX.save" \
  --scf-output "$WORK/cases/$CASE/scf.out" \
  --executable "$WORK/label-build/export_labels.x" \
  --output "$WORK/labels/$CASE"
```

Repeat with each converged case name, from an environment where the QE executable and its libraries can run.

## 4. Fit

```sh
cp "$REPO/4_potential/fit.toml" "$WORK/fit.toml"
```

Edit `[data].labels` to list the extracted label directories and set one positive `[data].weights` entry per case. Review the remaining fit settings, then run:

```sh
"$PY" "$REPO/4_potential/fit.py" \
  --config "$WORK/fit.toml" --output "$WORK/fit"
```

For continuation, repeat the command with `--resume`; labels, settings and code must be unchanged. Inspect `RESULT.json`: exit code 2 means the optimizer did not report success.

## 5. Export the fitted model

```sh
BUNDLE=$("$PY" -c 'import json,sys; print(json.load(open(sys.argv[1]))["bundle"])' "$WORK/fit/RESULT.json")
SCALAR_OFFSET_EV=0
"$PY" "$REPO/local_vfunc.py" export \
  --source "$BUNDLE" --code "$REPO/4_potential/local" \
  --structure "$BUNDLE/structure.npz" \
  --scalar-offset-eV "$SCALAR_OFFSET_EV" --output "$WORK/model"
```

Zero offset retains the fitted mean-zero scalar gauge. Before combining with an absolute-energy background, set `SCALAR_OFFSET_EV` to the physically matched offset; the fit does not determine it.
