MODE = PARALLEL
#MODE = SEQUENCIAL
DEBUG = NO
#DEBUG = YES

ifeq ($(MODE), PARALLEL)
FC = mpif90
FD = -cpp -D'PARALLEL'
endif

ifeq ($(MODE), SEQUENCIAL)
#FC = lf95
#FC = f95-lah
FC = gfortran
FD = -cpp
endif
#OPTS = -O3 -M OBJS  #-fp-model strict
OPTS = -O3 -J OBJS  #-fp-model strict
ifeq ($(DEBUG), YES)
#OPTS= -O0 -cpp -g -J OBJS -fbounds-check
OPTS= -O0 -cpp -g -J OBJS -fbounds-check -ffpe-trap=invalid,zero,overflow
endif
#OPTZ = -O3 -module OBJS  #-fp-model strict 

OBJS =  OBJS/from_nrtype.o \
	OBJS/orbit_mod.o \
	OBJS/mesh_mod.o \
	OBJS/sub_orbit_fast.o \
	OBJS/lpk_rngstream.o \
	OBJS/born_particles.o \
	OBJS/perp_step.o \
	OBJS/check_in_tri.o \
	OBJS/macrostep.o \
	OBJS/macrostep_ba.o \
	OBJS/collis.o \
	OBJS/arnoldi.o \
	OBJS/kin2d.o

kin2d.x: $(OBJS) K2d_lf.mk SRC/orbit_mod.f90 SRC/mesh_mod.f90
	$(FC) $(OPTS) -o kin2d.x $(OBJS) -llapack
OBJS/from_nrtype.o: SRC/from_nrtype.f90 K2d_lf.mk
	$(FC) $(OPTS) -c SRC/from_nrtype.f90
	mv from_nrtype.o OBJS/
OBJS/orbit_mod.o: SRC/orbit_mod.f90 K2d_lf.mk 
	$(FC) $(OPTS) -c SRC/orbit_mod.f90
	mv orbit_mod.o OBJS/
OBJS/mesh_mod.o: SRC/mesh_mod.f90 K2d_lf.mk  SRC/orbit_mod.f90
	$(FC) $(OPTS) -c SRC/mesh_mod.f90
	mv mesh_mod.o OBJS/
OBJS/sub_orbit_fast.o: SRC/sub_orbit_fast.f90 K2d_lf.mk SRC/orbit_mod.f90 SRC/mesh_mod.f90
	$(FC) $(OPTS) -c SRC/sub_orbit_fast.f90
	mv sub_orbit_fast.o OBJS/
OBJS/lpk_rngstream.o: SRC/lpk_rngstream.f90 K2d_lf.mk
	$(FC) $(OPTS) -c SRC/lpk_rngstream.f90
	mv lpk_rngstream.o OBJS/
OBJS/born_particles.o: SRC/born_particles.f90 K2d_lf.mk SRC/mesh_mod.f90
	$(FC) $(OPTS) -c SRC/born_particles.f90
	mv born_particles.o OBJS/
OBJS/perp_step.o: SRC/perp_step.f90 K2d_lf.mk SRC/mesh_mod.f90
	$(FC) $(OPTS) -c SRC/perp_step.f90
	mv perp_step.o OBJS/
OBJS/check_in_tri.o: SRC/check_in_tri.f90 K2d_lf.mk SRC/mesh_mod.f90 SRC/orbit_mod.f90
	$(FC) $(OPTS) -c SRC/check_in_tri.f90
	mv check_in_tri.o OBJS/
OBJS/macrostep.o: SRC/macrostep.f90 K2d_lf.mk SRC/mesh_mod.f90 SRC/orbit_mod.f90
	$(FC) $(OPTS) -c SRC/macrostep.f90
	mv macrostep.o OBJS/
OBJS/macrostep_ba.o: SRC/macrostep_ba.f90 K2d_lf.mk SRC/mesh_mod.f90 SRC/orbit_mod.f90
	$(FC) $(OPTS) -c SRC/macrostep_ba.f90
	mv macrostep_ba.o OBJS/
OBJS/collis.o: SRC/collis.f90 K2d_lf.mk SRC/orbit_mod.f90
	$(FC) $(OPTS) -c SRC/collis.f90
	mv collis.o OBJS/
OBJS/arnoldi.o: SRC/arnoldi.f90 K2d_lf.mk 
	$(FC) $(FD) $(OPTS) -c SRC/arnoldi.f90
	mv arnoldi.o OBJS/
OBJS/kin2d.o: SRC/kin2d.f90 K2d_lf.mk SRC/orbit_mod.f90 SRC/mesh_mod.f90 SRC/arnoldi.f90
	$(FC) $(FD) $(OPTS) -c SRC/kin2d.f90
	mv kin2d.o OBJS/
#OBJS/.o: .f90 K2d_lf.mk orbit_mod.f90
#	$(FC) $(OPTS) -c .f90
#	mv .o OBJS/
