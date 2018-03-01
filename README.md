NEO-EQ: 3D plasma equilibria with kinetics
==========================================

Preparing input data
--------------------
Enter directory
    
    cd PRELOAD/

Prepare mesh
    
    ./prepare_mesh.sh

Prepare vacuum perturbation field:
edit `VACFIELD/field_divB0.inp` to point to correct field.dat
(e.g. copied from `/proj/plasma/RMP/DATA/ASDEX/MESH3D_30835/`)
    
    ./prepare_vacfield.sh

MHD single iteration
--------------------
Enter directory    
    
    cd MHD/

Build MHD code
    
    make

Run MHD code
    
    ./magdif_test.sh
