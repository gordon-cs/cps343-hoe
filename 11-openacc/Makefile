CC	:= pgcc
CXX	:= pgc++
CFLAGS	:= -fast
CXXFLAGS:= $(CFLAGS)
ACCOPTS	:= -Minfo=accel -ta=tesla:cc30,managed -ta=tesla:cc50,managed

CCSRCS	:= pm_seq.cc pm_acc_kernels.cc pm_acc_loop.cc

OBJS	:= wtime.o readMatrixSeq.o $(CCSRCS:.cc=.o)
BINS	:= $(CCSRCS:.cc=)
LDLIBS	:= -L$(shell echo $(LIBRARY_PATH) | sed -e 's/:/\n/g' | grep pgi) -lhdf5

.PHONY:		all
all: $(BINS)

pm_seq:		pm_seq.cc wtime.o readMatrixSeq.o

pm_acc_kernels:	CXXFLAGS += $(ACCOPTS)
pm_acc_kernels:	pm_acc_kernels.cc wtime.o readMatrixSeq.o

pm_acc_loop:	CXXFLAGS += $(ACCOPTS)
pm_acc_loop:	pm_acc_loop.cc wtime.o readMatrixSeq.o

wtime.o:		wtime.c wtime.h

readMatrixSeq.o:	readMatrixSeq.cc

.PHONY:		clean clobber
clean clobber:
	$(RM) $(BINS) $(OBJS)
