#!/bin/bash
# Bordered shielding solve in the KiLCA cylinder: exact ideal shielding, the same outer field as the
# pressure-jump sheet, and partial penetration through a slowly rotating resistive layer.
# Usage: run_test.sh BUILD_DIR SOURCE_DIR WORK_DIR
set -eu
build=$1 src=$2 work=$3
python3 -c 'import h5py, numpy' 2>/dev/null || { echo "h5py/numpy missing"; exit 77; }
rm -rf "$work"
"$build/scripts/mephit.bash" init -c "$src/test/bordered_cylinder/mephit.in" \
    -g "$src/data/g000001.0001_TCFP_hip" -d kilca "$work/jump" > /dev/null
python3 "$src/test/sheet_cylinder/make_inputs.py" "$work/jump"
cp -r "$work/jump" "$work/ideal"
cp -r "$work/jump" "$work/layer"
sed -i 's|^   config%runmode = 3|   config%runmode = 3\n   config%bordered_shielding = .true.|' "$work/ideal/mephit.in"
sed -i 's|^   config%runmode = 3|   config%runmode = 3\n   config%bordered_shielding = .true.\n   config%layer_model = 1\n   config%layer_omega_factor = 1d-3|' "$work/layer/mephit.in"
# uniform E_r of 3 statV/cm gives the layer a finite rotation frequency
python3 -c "import numpy as np; x = np.linspace(0, 1, 101); np.savetxt('$work/layer/Er.dat', np.c_[x, np.full_like(x, 3.0)])"
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-2}
unset DISPLAY
for case in jump ideal layer; do
    "$build/scripts/mephit.bash" run "$work/$case" > "$work/$case/run.out" 2>&1 &
done
fail=0
for job in $(jobs -p); do wait "$job" || fail=1; done
[ "$fail" = 0 ] || { tail -5 "$work"/*/run.out; exit 1; }
python3 "$src/test/bordered_cylinder/check.py" "$work"
