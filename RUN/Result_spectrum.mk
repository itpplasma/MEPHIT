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
	OBJS/sub_orbit_fast.o \
	OBJS/check_in_tri.o \
	OBJS/result_spectrum.o

result_spectrum.x: $(OBJS) Result_spectrum.mk SRC/orbit_mod.f90 SRC/mesh_mod.f90
	$(FC) $(OPTS) -o result_spectrum.x $(OBJS)
OBJS/from_nrtype.o: SRC/from_nrtype.f90 K2d_lf.mk
	$(FC) $(OPTS) -c SRC/from_nrtype.f90
	mv from_nrtype.o OBJS/
OBJS/orbit_mod.o: SRC/orbit_mod.f90 K2d_lf.mk 
	$(FC) $(OPTS) -c SRC/orbit_mod.f90
	mv orbit_mod.o OBJS/
OBJS/sub_orbit_fast.o: SRC/sub_orbit_fast.f90 K2d_lf.mk SRC/orbit_mod.f90 SRC/mesh_mod.f90
	$(FC) $(OPTS) -c SRC/sub_orbit_fast.f90
	mv sub_orbit_fast.o OBJS/
OBJS/mesh_mod.o: SRC/mesh_mod.f90 K2d_lf.mk  SRC/orbit_mod.f90
	$(FC) $(OPTS) -c SRC/mesh_mod.f90
	mv mesh_mod.o OBJS/
OBJS/check_in_tri.o: SRC/check_in_tri.f90 K2d_lf.mk SRC/mesh_mod.f90 SRC/orbit_mod.f90
	$(FC) $(OPTS) -c SRC/check_in_tri.f90
	mv check_in_tri.o OBJS/
OBJS/result_spectrum.o: SRC/result_spectrum.f90 K2d_lf.mk SRC/orbit_mod.f90 SRC/mesh_mod.f90
	$(FC) $(FD) $(OPTS) -c SRC/result_spectrum.f90
	mv result_spectrum.o OBJS/
