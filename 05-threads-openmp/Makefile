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
CFLAGS	= -O2 -Wall -pedantic
CXXFLAGS= $(CFLAGS)
LDLIBS	=

#
# Sources
#
SRCS	= thread.cc bounded-buffer.cc bank.cc matmat_pt_static.cc \
		matmat_pt_dynamic.cc matmat_serial.cc \
		helmholtz_jacobi.cc
#
# Build targets
#
BINS	= $(SRCS:.cc=)
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
# Non-default dependencies for individual programs and objects
#
randomInt.o:	randomInt.cc randomInt.h
bank:		bank.cc randomInt.o

bank:			LDLIBS += -lpthread
thread:			LDLIBS += -lpthread
bounded-buffer:		LDLIBS += -lpthread
matmat_pt_dynamic:	CXXFLAGS += -DN=1000
matmat_pt_dynamic:	LDLIBS += -lpthread
matmat_pt_static:	CXXFLAGS += -DN=1000
matmat_pt_static:	LDLIBS += -lpthread

