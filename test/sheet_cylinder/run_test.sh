#!/bin/bash
# Sheet response of a resonant surface current in the KiLCA cylinder against the analytic
# field of the same current. Usage: run_test.sh BUILD_DIR SOURCE_DIR WORK_DIR
set -eu
build=$1 src=$2 work=$3
python3 -c 'import h5py, numpy, scipy' 2>/dev/null || { echo "h5py/numpy/scipy missing"; exit 77; }
rm -rf "$work"
"$build/scripts/mephit.bash" init -c "$src/test/sheet_cylinder/mephit.in" \
    -g "$src/data/g000001.0001_TCFP_hip" -d kilca "$work" > /dev/null
python3 "$src/test/sheet_cylinder/make_inputs.py" "$work"
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
unset DISPLAY
"$build/scripts/mephit.bash" run "$work" > "$work/run.out" 2>&1
"$build/scripts/mephit.bash" run --sheet "$work" > "$work/sheet.out" 2>&1
python3 "$src/test/sheet_cylinder/check.py" "$work"
