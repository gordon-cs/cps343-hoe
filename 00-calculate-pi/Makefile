# Makefile for serial and parallel pi computation examples
#
# Jonathan Senning <jonathan.senning@gordon.edu>
# Department of Mathematics and Computer Science
# Gordon College, 255 Grapevine Road, Wenham MA 01984-1899
#
# This file is released into the public domain.

CC	= gcc
CXX	= g++
CFLAGS	= -Wall -O2
CXXFLAGS= $(CFLAGS)

SRCS	= pi_serial.cc pi_omp.cc pi_omp_dyn.cc pi_mpi.cc
BINS	= $(SRCS:.cc=)

.PHONY:		all
all: $(BINS)

pi_serial:	LDLIBS += -lrt
pi_serial:	pi_serial.cc wtime.o

pi_omp:		CXXFLAGS += -fopenmp
pi_omp:		pi_omp.cc

pi_omp_dyn:	CXXFLAGS += -fopenmp
pi_omp_dyn:	pi_omp_dyn.cc

pi_mpi:		CXX = mpic++
pi_mpi:		pi_mpi.cc

wtime.o:	wtime.c wtime.h

.PHONY:		clean clobber
clean clobber:
	$(RM) $(BINS) wtime.o
