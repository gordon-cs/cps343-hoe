# Makefile for C and Fortran matrix-matrix multiplication programs and
# C++ matrix as 1D or 2D array program.
#
# Jonathan Senning <jonathan.senning@gordon.edu>
# Department of Mathematics and Computer Science
# Gordon College, 255 Grapevine Road, Wenham MA 01984-1899
#
# This file is released into the public domain.

#
# Compiler variables
#
CC	= gcc
CFLAGS	= -O3 -Wall -pedantic

#
# Executables produced by this Makefile
#
BINS = args args_req matrix_prod

#
# Build targets
#
.PHONY:	all
all: $(BINS)

#
# Targets for cleaning up; both do the same thing
#
.PHONY:	clean clobber
clean clobber:
	$(RM) $(BINS)

#
# LINK.c command explicitly listed since target name does not match
# source name minus extension.
#
args_req: CFLAGS += -DARGS_REQUIRED
args_req args: args.c
	$(LINK.c) $^ $(LOADLIBES) $(LDLIBS) -o $@

#
# Determine if BLAS (ATLAS BLAS or OpenBLAS) are available
#
ifneq (,$(findstring atlas,$(LMOD_FAMILY_BLAS)))
  matrix_prod:	CFLAGS += -DHAVE_BLAS -DHAVE_CBLAS
  matrix_prod:	LDLIBS += -lcblas -lf77blas -latlas -lgfortran
else ifneq (,$(findstring openblas,$(LMOD_FAMILY_BLAS)))
  matrix_prod:	CFLAGS += -DHAVE_BLAS
  matrix_prod:	LDLIBS += -lopenblas
else
  matrix_prod:	CFLAGS += -DHAVE_GSL_CBLAS
  matrix_prod:	LDLIBS += -lgslcblas

endif
