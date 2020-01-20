# Makefile for 00-intro-to-hpc/matrix_matrix_mutiply

# Set compilers and compiler flags
CXX	= g++
FC	= gfortran
CXXFLAGS= -Wall -O3 -funroll-loops
FFLAGS	= $(CXXFLAGS)

# List of sources and files to build from them
CXXSRCS	= matmat_cxx.cc matmat_ijk.cc
FSRCS	= matmat_f77.f
BINS	= $(CXXSRCS:.cc=) $(FSRCS:.f=)

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
matmat_cxx:	matmat_cxx.cc
matmat_ijk:	matmat_ijk.cc
matmat_f77:	matmat_f77.f
