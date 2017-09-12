# compiler

TOOLCHAIN = gnu
PROJLIBS = $(HOME)/CODE/Libs

ifeq ($(TOOLCHAIN),intel)
	FC := ifort
	CC := icc
	FCFLAGS = -c -g -O0 -mkl -traceback -check bounds #-implicitnone
	CCFLAGS = -c -g -O0 -mkl -traceback 
else ifeq ($(TOOLCHAIN),gnu)
	FC := gfortran
	CC := gcc
	FCFLAGS = -c -g -Og -Wall -fbacktrace -fcheck=bounds #-fimplicit-none
	CCFLAGS = -c -g -Og -Wall
	PROJLIBS := $(PROJLIBS)/GNU
else ifeq ($(TOOLCHAIN),sgi)
	FC := gfortran
	#FC := f90
	CC := gcc
	FCFLAGS = -c -g -Wall -fbacktrace -fcheck=bounds 
        #FCFLAGS = -c -g -O0 -Wall
	CCFLAGS = -std=c99 -c -g -O0 -Wall
	PROJLIBS := $(PROJLIBS)/GNU
endif


SUPERLU_DIR = $(PROJLIBS)/SuperLU_5.2.1
SUPERLU_LIB = $(SUPERLU_DIR)/lib
SUPERLU_HDR = $(SUPERLU_DIR)/SRC
SUPERLU_F90 = $(SUPERLU_DIR)/FORTRAN

SUITESPARSE_DIR = $(PROJLIBS)/SuiteSparse-4.5.5
SUITESPARSE_LIB = $(SUITESPARSE_DIR)/lib
SUITESPARSE_HDR = $(SUITESPARSE_DIR)/include
SUITESPARSE_F90 = $(SUITESPARSE_DIR)/UMFPACK/Demo

ifeq ($(OS),Windows_NT)
	SUPERLU_DIR     = /home/chral/src/SuperLU_5.2.1
	SUITESPARSE_LIB = /usr/lib
	SUITESPARSE_HDR = /usr/include/suitesparse
	SUITESPARSE_F90 = /usr/src/umfpack-5.7.4-1.src/umfpack-5.7.4/Demo
endif

# link flags
LDFLAGS = -L $(SUITESPARSE_LIB) -L $(SUPERLU_LIB) -lsuperlu -lumfpack -lamd -lcholmod \
          -lcolamd -lcamd -lccolamd  

ifeq ($(TOOLCHAIN),intel)
	LDFLAGS += -limf -lirc -mkl
else ifeq ($(TOOLCHAIN),gnu)
	ifeq ($(OS),Windows_NT)
		LDFLAGS += -lblas -llapack
	else
		LDFLAGS += -lopenblas -lpthread
	endif
else ifeq ($(TOOLCHAIN),sgi)
	LDFLAGS = -L $(SUITESPARSE_LIB) -L $(SUPERLU_LIB) -lsuperlu -lumfpack -lamd -lcholmod \
          -lcolamd -lcamd -lccolamd -lblas
endif


# source files and objects
SRC = from_nrtype.f90 spline5_RZ.f90 field_divB0.f90 bdivfree.f90 orbit_mod.f90 mesh_mod.f90 \
      sparse_mod.f90 magdif.f90 magdif_test.f90
SRCS = $(patsubst %.f90, %.o, $(SRC))

# program name
PROGRAM = magdif_test

.PHONY: all clean

all: $(PROGRAM)

$(PROGRAM): $(SRCS) c_fortran_dgssv.o c_fortran_zgssv.o umf4_f77wrapper.o umf4_f77zwrapper.o
	$(FC) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(FC) -o $@ $< $(FCFLAGS) 

c_fortran_dgssv.o: $(SUPERLU_F90)/c_fortran_dgssv.c
	$(CC) $(CCFLAGS) -I $(SUPERLU_HDR) -o $@ $^

c_fortran_zgssv.o: $(SUPERLU_F90)/c_fortran_zgssv.c
	$(CC) $(CCFLAGS) -I $(SUPERLU_HDR) -o $@ $^

umf4_f77wrapper.o: $(SUITESPARSE_F90)/umf4_f77wrapper.c
	$(CC) $(CCFLAGS) -I $(SUITESPARSE_HDR) -DDLONG -o $@ $^

umf4_f77zwrapper.o: $(SUITESPARSE_F90)/umf4_f77zwrapper.c
	$(CC) $(CCFLAGS) -I $(SUITESPARSE_HDR) -DZLONG -o $@ $^

clean:
	rm -f *.o *.mod magdif_test magdif_test.exe
