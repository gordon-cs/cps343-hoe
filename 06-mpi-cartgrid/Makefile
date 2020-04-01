CXX	:= mpic++
CXXFLAGS:= -Wall -g -O3 -funroll-loops
LDLIBS	:= -lrt

# Project files

SRCS	:= cart.cc cart-hdf5.cc cart-nb.cc \
		laplace-seq.cc laplace-mpi.cc laplace-mpi-nb.cc
BINS 	:= $(SRCS:.cc=)
HDRS	:= cartesian-mpi.h
OBJS	:= cartesian-mpi.o

# ============================================================================
# ============================== TARGET SECTION ==============================
# ============================================================================

# Program targets

.PHONY:	all
all:	$(BINS)

cart:			cart.cc
cart-nb:		cart-nb.cc

cart-hdf5:		LDLIBS += -lhdf5
cart-hdf5:		cart-hdf5.cc

laplace-seq:		CXX = g++
laplace-seq:		laplace-seq.cc

laplace-mpi:		laplace-mpi.cc $(OBJS)
laplace-mpi-nb: 	laplace-mpi-nb.cc $(OBJS)

cartesian-mpi.o:	cartesian-mpi.cc $(HDRS)

# Supporting targets

.PHONY: clean
clean:
	$(RM) $(OBJS)

.PHONY: clobber
clobber:
	$(MAKE) clean
	$(RM) $(BINS)
