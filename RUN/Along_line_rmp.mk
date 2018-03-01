#MODE = PARALLEL
MODE = SEQUENCIAL
DEBUG = YES
#DEBUG = YES

ifeq ($(MODE), PARALLEL)
FC = mpif90
FD = -Cpp -D'PARALLEL'
endif

ifeq ($(MODE), SEQUENCIAL)
#FC = lf95
#FC = f95-lah
FC = gfortran
#FD = -Cpp
FD = 
endif
#OPTS = -O3 -M OBJS  #-fp-model strict
OPTS = -O -J OBJS  #-fp-model strict
ifeq ($(DEBUG), YES)
OPTS= -O0 --trace -g  -J  OBJS   
endif
#OPTZ = -O3 -module OBJS  #-fp-model strict 

OBJS =  OBJS/from_nrtype.o \
	OBJS/orbit_mod.o \
	OBJS/mesh_mod.o \
	OBJS/check_in_tri.o \
	OBJS/along_line_rmp.o

along_line_rmp.x: $(OBJS) Along_line_rmp.mk SRC/orbit_mod.f90 SRC/mesh_mod.f90
	$(FC) $(OPTS) -o along_line_rmp.x $(OBJS)
OBJS/from_nrtype.o: SRC/from_nrtype.f90 Along_line_rmp.mk
	$(FC) $(OPTS) -c SRC/from_nrtype.f90
	mv from_nrtype.o OBJS/
OBJS/orbit_mod.o: SRC/orbit_mod.f90 Along_line_rmp.mk 
	$(FC) $(OPTS) -c SRC/orbit_mod.f90
	mv orbit_mod.o OBJS/
OBJS/sub_orbit_fast.o: SRC/sub_orbit_fast.f90 Along_line_rmp.mk SRC/orbit_mod.f90 SRC/mesh_mod.f90
	$(FC) $(OPTS) -c SRC/sub_orbit_fast.f90
	mv sub_orbit_fast.o OBJS/
OBJS/mesh_mod.o: SRC/mesh_mod.f90 Along_line_rmp.mk  SRC/orbit_mod.f90
	$(FC) $(OPTS) -c SRC/mesh_mod.f90
	mv mesh_mod.o OBJS/
OBJS/check_in_tri.o: SRC/check_in_tri.f90 Along_line_rmp.mk  SRC/orbit_mod.f90
	$(FC) $(OPTS) -c SRC/check_in_tri.f90
	mv check_in_tri.o OBJS/
OBJS/along_line_rmp.o: SRC/along_line_rmp.f90 Along_line_rmp.mk SRC/orbit_mod.f90 SRC/mesh_mod.f90
	$(FC) $(FD) $(OPTS) -c SRC/along_line_rmp.f90
	mv along_line_rmp.o OBJS/
