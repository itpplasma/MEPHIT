#FC       = lf95 
FC       = gfortran 
#FC       = ifort

#OPTS= -O0 --chk a,e,s,u,x --trace --trap -g -M OBJS
#OPTS= -O -M OBJS
#OPTS= -O -J OBJS
OPTS = -Og -g -fbacktrace -fimplicit-none -J OBJS
#OPTS = -O0 -check all -check noarg_temp_created -check uninit -g -traceback -fpe0 -module OBJS

OBJS =  OBJS/from_nrtype.o \
	OBJS/orbit_mod.o \
	OBJS/mesh_mod.o \
	OBJS/writeout.o

writeout.x: $(OBJS) Writeout.mk
	$(FC) $(OPTS) -o writeout.x $(OBJS)
OBJS/from_nrtype.o: SRC/from_nrtype.f90 Writeout.mk
	$(FC) $(OPTS) -c SRC/from_nrtype.f90
	mv from_nrtype.o OBJS/
OBJS/orbit_mod.o: SRC/orbit_mod.f90 Writeout.mk
	$(FC) $(OPTS) -c SRC/orbit_mod.f90
	mv orbit_mod.o OBJS/
OBJS/mesh_mod.o: SRC/mesh_mod.f90 Writeout.mk
	$(FC) $(OPTS) -c SRC/mesh_mod.f90
	mv mesh_mod.o OBJS/
OBJS/writeout.o: SRC/writeout.f90 Writeout.mk
	$(FC) $(OPTS) -c SRC/writeout.f90
	mv writeout.o OBJS/
