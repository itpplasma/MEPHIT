FC       := gfortran
#FC       := lf95
#FC       := ifort

#FFLAGS := --chk a,e,s,u,x --trace --trap -g -M OBJS
#FFLAGS := -O -M OBJS
#FFLAGS := -O -J OBJS
#FFLAGS := -O0 -check all -check noarg_temp_created -check uninit -g -traceback -fpe0 -module OBJS

FFLAGS := -Og -g -fbacktrace -J OBJS
# -fimplicit-none does not work with field_divB0.f90

AXIS_SRCS := SRC/field_divB0.f90 \
	SRC/bdivfree.f90 \
	SRC/from_nrtype.f90 \
	SRC/brent.f90 \
	SRC/rtbis.f90 \
	SRC/odeint.f \
	SRC/rkck.f \
	SRC/rkqs.f \
	SRC/spline5_RZ.f90 \
	SRC/leastsq.f90 \
	SRC/lubksb.f \
	SRC/ludcmp.f \
	SRC/polylag.f \
	SRC/twostreights.f90 \
	SRC/mods4ax.f90 \
	SRC/in_out_cw.f90 \
	SRC/inside_of_wall.f90 \
	SRC/axis.f90

TRI_SRCS := SRC/tri_mesh.f90

RC_SRCS := SRC/from_nrtype.f90 \
	SRC/orbit_mod.f90 \
	SRC/mesh_mod.f90 \
	SRC/readcarre_m.f90

WRITEOUT_SRCS := SRC/from_nrtype.f90 \
	SRC/orbit_mod.f90 \
	SRC/mesh_mod.f90 \
	SRC/writeout.f90

VACFIELD_SRCS := SRC/field_divB0.f90 \
	SRC/bdivfree.f90 \
	SRC/spline5_RZ.f90 \
	SRC/from_nrtype.f90 \
	SRC/orbit_mod.f90 \
	SRC/mesh_mod.f90 \
	SRC/vacfield.f90

AXIS_OBJS := $(patsubst SRC/%.f, OBJS/%.o, $(patsubst SRC/%.f90, OBJS/%.o, $(AXIS_SRCS)))
TRI_OBJS := $(patsubst SRC/%.f, OBJS/%.o, $(patsubst SRC/%.f90, OBJS/%.o, $(TRI_SRCS)))
RC_OBJS := $(patsubst SRC/%.f, OBJS/%.o, $(patsubst SRC/%.f90, OBJS/%.o, $(RC_SRCS)))
WRITEOUT_OBJS := $(patsubst SRC/%.f, OBJS/%.o, $(patsubst SRC/%.f90, OBJS/%.o, $(WRITEOUT_SRCS)))
VACFIELD_OBJS := $(patsubst SRC/%.f, OBJS/%.o, $(patsubst SRC/%.f90, OBJS/%.o, $(VACFIELD_SRCS)))

AXIS_EXE := axis.x
TRI_EXE := tri_mesh.x
RC_EXE := readcarre_m.x
WRITEOUT_EXE := writeout.x
VACFIELD_EXE := vacfield.x

FORTRAN_MESH_FILE := inputformaxwell.msh
FREEFEM_MESH_FILE := inputformaxwell_ext.msh
VACFIELD_FILE := VACFIELD/Bn_flux.dat

.PHONY: asdex kilca clean

asdex: $(VACFIELD_FILE) $(FORTRAN_MESH_FILE) $(FREEFEM_MESH_FILE)

kilca:
	../BUILD/bin/minimal_example.x ../MHD/magdif_minimal_example.in
	FreeFem++ extmesh.edp

clean:
	rm -f $(AXIS_OBJS) $(TRI_OBJS) $(RC_OBJS) $(WRITEOUT_OBJS) $(VACFIELD_OBJS) $(FORTRAN_MESH_FILE) $(FREEFEM_MESH_FILE)

$(VACFIELD_FILE): $(VACFIELD_EXE) $(FORTRAN_MESH_FILE) Makefile VACFIELD/field_divB0.inp
	cd VACFIELD && ./$(VACFIELD_EXE) && cd ..

$(FREEFEM_MESH_FILE): $(FORTRAN_MESH_FILE) extmesh.py Makefile
	python extmesh.py $(basename $(FORTRAN_MESH_FILE))

$(FORTRAN_MESH_FILE): $(AXIS_EXE) $(TRI_EXE) $(RC_EXE) $(WRITEOUT_EXE) Makefile field_divB0.inp
	./$(AXIS_EXE)
	./$(TRI_EXE)
	./$(RC_EXE)
	./$(WRITEOUT_EXE)

$(AXIS_EXE): $(AXIS_OBJS) Makefile
	$(FC) $(FFLAGS) -o $(AXIS_EXE) $(AXIS_OBJS)

$(TRI_EXE): $(TRI_OBJS) Makefile
	$(FC) $(FFLAGS) -o $(TRI_EXE) $(TRI_OBJS)

$(RC_EXE): $(RC_OBJS) Makefile
	$(FC) $(FFLAGS) -o $(RC_EXE) $(RC_OBJS)

$(VACFIELD_EXE): $(VACFIELD_OBJS) Makefile
	$(FC) $(FFLAGS) -o $(VACFIELD_EXE) $(VACFIELD_OBJS)

$(WRITEOUT_EXE): $(WRITEOUT_OBJS) Makefile
	$(FC) $(FFLAGS) -o $(WRITEOUT_EXE) $(WRITEOUT_OBJS)

OBJS/%.o: SRC/%.f90 Makefile
	$(FC) $(FFLAGS) -o $@ -c $<

OBJS/%.o: SRC/%.f Makefile
	$(FC) $(FFLAGS) -o $@ -c $<
