# compiler
FC := gfortran
CC := gcc

PROJLIBS = $(HOME)/CODE/Libs/
SUITESPARSE_DIR = $(PROJLIBS)SuiteSparse-4.5.3/
SUITESPARSE_LIB = $(SUITESPARSE_DIR)lib/
SUITESPARSE_HDR = $(SUITESPARSE_DIR)include/
SUITESPARSE_F90 = $(SUITESPARSE_DIR)UMFPACK/Demo/
SUITESPARSE_CFLAGS = -I $(SUITESPARSE_HDR)

# compile flags
FCFLAGS = -g -c -fdefault-real-8 -fbacktrace -fno-align-commons -fbounds-check
CFLAGS = 
# link flags
FLFLAGS = -L $(SUITESPARSE_DIR)/lib -lumfpack

# source files and objects
SRCS = $(patsubst %.f90, %.o, $(wildcard *.f90))
#       $(patsubst %.h, %.mod, $(wildcard *.h))

# program name
PROGRAM = magdif_test

all: $(PROGRAM)

$(PROGRAM): $(SRCS) umf4_f77wrapper.o umf4_f77zwrapper.o
	$(FC) $(FLFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(FCFLAGS) -o $@ $<

umf4_f77wrapper.o: $(SUITESPARSE_F90)umf4_f77wrapper.c
	$(CC) $(CFLAGS) -c $(SUITESPARSE_CFLAGS) -o $@ $^ -DDLONG

umf4_f77zwrapper.o: $(SUITESPARSE_F90)umf4_f77zwrapper.c
	$(CC) $(CFLAGS) -c $(SUITESPARSE_CFLAGS) -o $@ $^ -DZLONG

clean:
	rm -f *.o *.mod
