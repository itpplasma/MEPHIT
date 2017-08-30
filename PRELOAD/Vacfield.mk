#FC       = lf95 --ap
FC       = gfortran

#OPTS= --chk a,e,s,u,x --trace --trap -g -M OBJS
#OPTS= -O -M OBJS
OPTS= -O -J OBJS

OBJS =  OBJS/field_divB0.o \
	OBJS/bdivfree.o \
	OBJS/spline5_RZ.o \
	OBJS/from_nrtype.o \
	OBJS/orbit_mod.o \
	OBJS/mesh_mod.o \
	OBJS/vacfield.o 


vacfield.x: $(OBJS) Vacfield.mk
	$(FC) $(OPTS) -o vacfield.x $(OBJS)
OBJS/field_divB0.o: SRC/field_divB0.f90 Vacfield.mk
	$(FC) $(OPTS) -c SRC/field_divB0.f90
	mv field_divB0.o OBJS/
OBJS/bdivfree.o: SRC/bdivfree.f90 Vacfield.mk
	$(FC) $(OPTS) -c SRC/bdivfree.f90
	mv bdivfree.o OBJS/
OBJS/spline5_RZ.o: SRC/spline5_RZ.f90 Vacfield.mk
	$(FC) $(OPTS) -c SRC/spline5_RZ.f90
	mv spline5_RZ.o OBJS/
OBJS/from_nrtype.o: SRC/from_nrtype.f90 Writeout.mk
	$(FC) $(OPTS) -c SRC/from_nrtype.f90
	mv from_nrtype.o OBJS/
OBJS/orbit_mod.o: SRC/orbit_mod.f90 Writeout.mk
	$(FC) $(OPTS) -c SRC/orbit_mod.f90
	mv orbit_mod.o OBJS/
OBJS/mesh_mod.o: SRC/mesh_mod.f90 Writeout.mk
	$(FC) $(OPTS) -c SRC/mesh_mod.f90
	mv mesh_mod.o OBJS/
OBJS/vacfield.o: SRC/vacfield.f90 Vacfield.mk
	$(FC) $(OPTS) -c SRC/vacfield.f90
	mv vacfield.o OBJS/
