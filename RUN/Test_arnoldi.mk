#MODE = PARALLEL
MODE = SEQUENCIAL
#DEBUG = NO
DEBUG = YES

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
OPTS= -O0 -cpp -g -J OBJS -fbounds-check
endif
#OPTZ = -O3 -module OBJS  #-fp-model strict 

OBJS =  OBJS/orbit_mod.o\
        OBJS/lpk_rngstream.o\
	OBJS/arnoldi.o \
	OBJS/test_arnoldi.o

test_arnoldi.x: $(OBJS) Test_arnoldi.mk SRC/test_arnoldi.f90 SRC/arnoldi.f90
	$(FC) $(OPTS) -o test_arnoldi.x $(OBJS) -lblas -llapack
OBJS/from_nrtype.o: SRC/from_nrtype.f90 Test_arnoldi.mk
	$(FC) $(OPTS) -c SRC/from_nrtype.f90
	mv from_nrtype.o OBJS/
OBJS/orbit_mod.o: SRC/orbit_mod.f90 Test_arnoldi.mk
	$(FC) $(OPTS) -c SRC/orbit_mod.f90
	mv orbit_mod.o OBJS/	
OBJS/lpk_rngstream.o: SRC/lpk_rngstream.f90 Test_arnoldi.mk
	$(FC) $(OPTS) -c SRC/lpk_rngstream.f90
	mv lpk_rngstream.o OBJS/
OBJS/arnoldi.o: SRC/arnoldi.f90 Test_arnoldi.mk
	$(FC) $(OPTS) -c SRC/arnoldi.f90
	mv arnoldi.o OBJS/
OBJS/test_arnoldi.o: SRC/test_arnoldi.f90 Test_arnoldi.mk SRC/arnoldi.f90 
	$(FC) $(FD) $(OPTS) -c SRC/test_arnoldi.f90
	mv test_arnoldi.o OBJS/
#OBJS/.o: .f90 K2d_lf.mk orbit_mod.f90
#	$(FC) $(OPTS) -c .f90
#	mv .o OBJS/
