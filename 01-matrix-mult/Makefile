# Makefile for C and Fortran matrix-matrix multiplication programs and
# C++ matrix as 1D or 2D array program.
#
# Jonathan Senning <jonathan.senning@gordon.edu>
# Department of Mathematics and Computer Science
# Gordon College, 255 Grapevine Road, Wenham MA 01984-1899
#
# This file is released into the public domain.

# Determine if PGI C compiler is available.  If so, we assume we want to
# compile programs using the PGI compiler suite.
PGCC	:= $(shell which pgcc 2> /dev/null)

# Define variables
ifdef PGCC
	CC	= pgcc
	CXX	= pgc++
	FC	= pgfortran
	F77	= pgfortran
	F95	= pgf95
	CFLAGS	= -fast
else
	CC	= gcc
	CXX	= g++
	FC	= gfortran
	F77	= gfortran -std=legacy
	F95	= gfortran -std=f95
	CFLAGS	= -Wall -O3 -funroll-loops
endif
CXXFLAGS= $(CFLAGS)
FFLAGS	= $(CFLAGS)

BINS = matmat_f77 matmat_f95 matmat_c1 matmat_c2 matmat_c3 matmat_c4 \
       matmat_ijk timing_1D timing_2D

#
# Build targets
#
.PHONY:		all
all: $(BINS)

#
# Include compiler flags to coax gfortran to behaving like F77 or F95
#
matmat_f77: FC = $(F77)
matmat_f77: matmat_f77.f

matmat_f95: FC = $(F95)
matmat_f95: matmat_f95.f95

matmat_c1: matmat_c1.c

matmat_c2: matmat_c2.c

matmat_c3: matmat_c3.c

matmat_c4: matmat_c4.c

matmat_ijk: matmat_ijk.c

#
# LINK.C command explicitly listed since target name does not match
# source name minus extension.
#
timing_1D: CFLAGS += -DUSE_LINEAR_ARRAY
timing_1D: timing.cc
	$(LINK.C) $^ $(LOADLIBES) $(LDLIBS) -o $@

timing_2D: timing.cc
	$(LINK.C) $^ $(LOADLIBES) $(LDLIBS) -o $@

#
# Targets for cleaning up
#
.PHONY:		clean clobber
clean clobber:
	$(RM) $(BINS)

#
# Rule for using Fortran 95 free-format source files
#
%: %.f95
	$(LINK.f) $^ $(LOADLIBES) $(LDLIBS) -o $@
