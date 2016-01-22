# Makefile for serial and parallel pi computation examples
#
# Jonathan Senning <jonathan.senning@gordon.edu>
# Department of Mathematics and Computer Science
# Gordon College, 255 Grapevine Road, Wenham MA 01984-1899
#
# This file is released into the public domain.

# Define variables
CC	= gcc
CFLAGS	= -Wall -02
CXX	= g++
CXXFLAGS= $(CFLAGS)

# List of sources and files to build from them
SRCS	= pi_serial.cc pi_omp.cc pi_omp_dyn.cc pi_mpi.cc
BINS	= $(SRCS:.cc=) # replaces .cc extension with empty string

# First (default) target
all: $(BINS)

# alternate targets used to clean up directory
clean:
	$(RM) wtime.o

clobber:
	$(MAKE) clean
	$(RM) $(BINS)

# Explicit dependencies and adjustments to variables
pi_serial:	CXXFLAGS += -lrt
pi_serial:	pi_serial.cc wtime.o

pi_omp:		CXXFLAGS += -fopenmp
pi_omp:		pi_omp.cc

pi_omp_dyn:	CXXFLAGS += -fopenmp
pi_omp_dyn:	pi_omp_dyn.cc

pi_mpi:		CXX = mpic++
pi_mpi:		pi_mpi.cc

wtime.o:	wtime.c wtime.h
