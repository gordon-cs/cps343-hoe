# Makefile for serial and parallel pi computation examples
#
# Jonathan Senning <jonathan.senning@gordon.edu>
# Department of Mathematics and Computer Science
# Gordon College, 255 Grapevine Road, Wenham MA 01984-1899
#
# This file is released into the public domain.

# Define variables
CC	:= gcc
NVCC	:= nvcc
OPTFLAGS:= -O3
NVFLAGS	:= -Xptxas -v
CFLAGS	:= $(OPTFLAGS)

# List of sources and files to build from them
CUSRCS	= matmul.cu
BINS	= matmul08 matmul16 matmul32

# First (default) target
all: $(BINS)

# alternate targets used to clean up directory
# (pi_serial.o is created by PGI compiler)
clean:
	$(RM) wtime.o

clobber:
	$(MAKE) clean
	$(RM) $(BINS)

# Explicit dependencies and adjustments to variables
matmul08:	CFLAGS += -DBS=8
matmul08:	matmul.cu wtime.o
	$(NVCC) $(NVFLAGS) $(CFLAGS) $^ -o $@

matmul16:	CFLAGS += -DBS=16
matmul16:	matmul.cu wtime.o
	$(NVCC) $(NVFLAGS) $(CFLAGS) $^ -o $@

matmul32:	CFLAGS += -DBS=32
matmul32:	matmul.cu wtime.o
	$(NVCC) $(NVFLAGS) $(CFLAGS) $^ -o $@

wtime.o:	CC = gcc
wtime.o:	CFLAGS = $(OPTFLAGS)
wtime.o:	wtime.c wtime.h
