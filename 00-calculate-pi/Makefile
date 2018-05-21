# Makefile for serial and parallel pi computation examples
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
	CFLAGS	= -fast
	OMPFLAGS= -mp=allcores
else
	CC	= gcc
	CXX	= g++
	CFLAGS	= -Wall -O3 -funroll-loops
	OMPFLAGS= -fopenmp
endif
CXXFLAGS= $(CFLAGS)
MPICC	= mpicc
MPICXX	= mpic++

# List of sources and files to build from them
SRCS	= pi_serial.cc pi_omp.cc pi_omp_dyn.cc pi_mpi.cc
OBJS	= $(SRCS:.cc=.o) # replaces .cc extension with .o
BINS	= $(SRCS:.cc=) # replaces .cc extension with empty string

# First (default) target
all: $(BINS)

# Alternate targets used to clean up directory
clean:
	$(RM) wtime.o $(OBJS)

clobber:
	$(MAKE) clean
	$(RM) $(BINS)

# Explicit dependencies and adjustments to variables
pi_serial:	pi_serial.cc wtime.o

pi_omp:		CXXFLAGS += $(OMPFLAGS)
pi_omp:		pi_omp.cc

pi_omp_dyn:	CXXFLAGS += $(OMPFLAGS)
pi_omp_dyn:	pi_omp_dyn.cc

pi_mpi:		CXX = $(MPICXX)
pi_mpi:		pi_mpi.cc

wtime.o:	wtime.c wtime.h
