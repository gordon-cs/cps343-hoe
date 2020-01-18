# Makefile for 00-intro-to-hpc

# Set compilers and compiler flags
CC	= gcc
FC	= gfortran
MPICC	= mpicc
CFLAGS	= -Wall -O3 -funroll-loops
FFLAGS	= $(CFLAGS)
OMPFLAGS= -fopenmp

# List of sources and files to build from them
CSRCS	= matmat_ijk.c matmat_c1.c matmat_c2.c matmat_c3.c matmat_c4.c
FSRCS	= matmat_f77.f
BINS	= $(CSRCS:.c=) $(FSRCS:.f=)

# First (default) target
.PHONY:	all
all:	$(BINS)

# Alternate targets used to clean up directory
.PHONY:	clean
clean:	# nothing to do

.PHONY:	clobber
clobber:
	$(RM) $(BINS)

# Explicit dependencies and adjustments to variables
matmat_c1:	matmat_c1.c
matmat_c2:	matmat_c2.c
matmat_c3:	matmat_c3.c
matmat_c4:	matmat_c4.c
matmat_ijk:	matmat_ijk.c
matmat_f77:	matmat_f77.f
