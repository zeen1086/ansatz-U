#!/usr/bin/env bash
set -euo pipefail
if [[ $# != 4 ]]; then
  echo 'Usage: run_joint.sh archives.tsv fit.nml NEW_RUN_DIR EXTERNAL_BUILD_DIR' >&2
  exit 2
fi
root=$(cd "$(dirname "$0")/../.." && pwd -P)
list=$(cd "$(dirname "$1")" && pwd -P)/$(basename "$1")
control=$(cd "$(dirname "$2")" && pwd -P)/$(basename "$2")
run=$(cd "$(dirname "$3")" && pwd -P)/$(basename "$3")
build=$(cd "$(dirname "$4")" && pwd -P)/$(basename "$4")
case "$run/" in "$root/"*) echo 'RUN must be outside the source tree' >&2; exit 2;; esac
case "$build/" in "$root/"*) echo 'BUILD must be outside the source tree' >&2; exit 2;; esac
[[ ! -e "$run" ]] || { echo 'Run directory exists; choose r02 or another new name' >&2; exit 2; }
command -v unzip >/dev/null
make -C "$root" joint "BUILD=$build"
mkdir "$run"
status=FAILED
trap 'printf "%s\n" "$status" > "$run/STATUS"' EXIT
mkdir "$run/inputs"
cp "$control" "$run/fit.nml"
cp "$list" "$run/source_archives.tsv"
rows=()
count=0
while IFS=$'\t' read -r name weight anchor unit archive extra || [[ -n "${name:-}" ]]; do
  [[ -n "$name" && "$name" != \#* ]] || continue
  [[ -n "$weight" && -n "$anchor" && -n "$unit" && -n "$archive" && -z "${extra:-}" ]] || {
    echo 'Expected five TAB-separated columns: name weight anchor unit archive' >&2; exit 2;
  }
  case "$archive" in /*) ;; *) archive=$(dirname "$list")/$archive;; esac
  [[ -f "$archive" ]] || { echo "Missing DFT archive: $archive" >&2; exit 2; }
  case "$unit" in
    Ry) scalar=v_scalar_Ry.npy; magnetic=B_material_Ry.npy;;
    eV) scalar=v_scalar_eV.npy; magnetic=B_material_eV.npy;;
    *) echo 'Field unit must be Ry or eV' >&2; exit 2;;
  esac
  count=$((count+1))
  folder="$run/inputs/state$count"
  mkdir "$folder"
  for member in grid.npy atom_ids.npy N.npy; do unzip -p "$archive" "$member" > "$folder/$member"; done
  unzip -p "$archive" "$scalar" > "$folder/scalar.npy"
  unzip -p "$archive" "$magnetic" > "$folder/B.npy"
  unzip -p "$archive" moments_muB.npy > "$folder/moments.npy"
  members=$(unzip -Z1 "$archive")
  case $'\n'"$members"$'\n' in *$'\nV_U_eV.npy\n'*) unzip -p "$archive" V_U_eV.npy > "$folder/VU.npy";; esac
  rows+=("$name $weight $anchor $unit" "$folder")

done < "$list"
[[ $count -gt 0 ]] || { echo 'No DFT states listed' >&2; exit 2; }
printf '%s\n' "$count" "${rows[@]}" > "$run/dataset.manifest"
if (cd "$run" && "$build/fit_joint" fit.nml dataset.manifest) > "$run/fit.log" 2>&1; then
  status=COMPLETE
else
  code=$?
  tail -n 30 "$run/fit.log" >&2
  exit "$code"
fi
printf 'Joint fit complete: %s\n' "$run/V_func_tmp"
