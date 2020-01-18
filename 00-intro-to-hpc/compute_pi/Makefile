# Makefile for serial and parallel pi computation examples

# Set compilers and compiler flags
CC	= gcc
CXX	= g++
MPICC	= mpicc
MPICXX	= mpic++
CFLAGS	= -Wall -O3 -funroll-loops
CXXFLAGS= $(CFLAGS)
OMPFLAGS= -fopenmp

# List of sources and files to build from them
BINS	= pi_serial pi_omp pi_mpi
OBJS	= wtime.o

# First (default) target
.PHONY:	all
all:	$(BINS)

# Alternate targets used to clean up directory
.PHONY:	clean
clean:
	$(RM) $(OBJS)

.PHONY:	clobber
clobber:
	$(RM) $(OBJS) $(BINS)

# Explicit dependencies and adjustments to variables
pi_serial:	pi_serial.cc wtime.o

pi_omp:		CXXFLAGS += $(OMPFLAGS)
pi_omp:		pi_omp.cc

pi_mpi:		CXX = $(MPICXX)
pi_mpi:		pi_mpi.cc

wtime.o:	wtime.c wtime.h
