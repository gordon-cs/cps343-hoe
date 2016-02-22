# Makefile for POSIX threads and semaphores examples
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
CXX	= g++
CFLAGS	= -O2 -Wall
CXXFLAGS= $(CFLAGS)
LDLIBS	= -lpthread -lm

#
# Build targets
#
BINS	= thread bounded-buffer bank matmat_pt_static matmat_pt_dynamic
OBJS	= randomInt.o

.PHONY: all
all:	$(BINS)

#
# Targets for cleaning up
#
.PHONY: clean clobber
clean clobber:
	$(RM) $(BINS) $(OBJS)

#
# Dependencies for individual programs and objects
#
thread:		thread.cc
bounded-buffer:	bounded-buffer.cc
bank:		bank.cc randomInt.o
matmat_pt_static:	matmat_pt_static.cc
matmat_pt_dynamic:	matmat_pt_dynamic.cc

randomInt.o:	randomInt.cc randomInt.h
