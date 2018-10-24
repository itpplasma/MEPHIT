NEO-EQ: 3D plasma equilibria with kinetics
==========================================

Building
--------
    cd BUILD/
    cmake ..
    make

Preparing input data
--------------------

Prepare vacuum perturbation field:
edit `PRELOAD/VACFIELD/field_divB0.inp` to point to correct `field.dat`
(e.g. copied from `/proj/plasma/RMP/DATA/ASDEX/MESH3D_30835/`)
    
    cd PRELOAD/
    make

MHD single iteration
--------------------
    cd MHD/

Check or edit `magdif.in`
    
    ./magdif_test.sh

