# Makefile for C and Fortran matrix-matrix multiplication programs and
# C++ matrix as 1D or 2D array program.
#
# Jonathan Senning <jonathan.senning@gordon.edu>
# Department of Mathematics and Computer Science
# Gordon College, 255 Grapevine Road, Wenham MA 01984-1899
#
# This file is released into the public domain.

GNUFLAGS= -Wall -O3

FC	= gfortran
FFLAGS	= $(GNUFLAGS)

CC	= gcc
CFLAGS	= $(GNUFLAGS)

CXX	= g++
CXXFLAGS= $(GNUFLAGS)

BINS = matmat_f77 matmat_f95 matmat_c1 matmat_c2 matmat_c3 matmat_ijk \
	timing_1D timing_2D# matmat_papi

#
# Build targets
#
.PHONY:		all
all: $(BINS)

#
# Include compiler flags to coax gfortran to behaving like F77 or F95
#
matmat_f77: FFLAGS += -std=legacy
matmat_f77: matmat_f77.f

matmat_f95: FFLAGS += -std=f95
matmat_f95: matmat_f95.f

#
# C programs need librt for timing
#
matmat_c1: LDLIBS += -lrt
matmat_c1: matmat_c1.c

matmat_c2: LDLIBS += -lrt
matmat_c2: matmat_c2.c

matmat_c3: LDLIBS += -lrt
matmat_c3: matmat_c3.c

matmat_ijk: LDLIBS += -lrt
matmat_ijk: matmat_ijk.c

matmat_papi: LDLIBS += -lrt -lpapi
matmat_papi: matmat_papi.c
#
# C++ needs librt for timing.  LINK.C command explicitly listed since target
# name does not match source name minus extension.
#
timing_1D: CXXFLAGS += -DUSE_MACRO
timing_1D: LDLIBS += -lrt
timing_1D: timing.cc
	$(LINK.C) $^ $(LOADLIBES) $(LDLIBS) -o $@

timing_2D: LDLIBS += -lrt
timing_2D: timing.cc
	$(LINK.C) $^ $(LOADLIBES) $(LDLIBS) -o $@

#
# Targets for cleaning up
#
.PHONY:		clean clobber
clean clobber:
	$(RM) $(BINS)
