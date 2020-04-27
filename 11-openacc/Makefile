CC	:= pgcc
CXX	:= pgc++
CFLAGS	:= -fast
CXXFLAGS:= $(CFLAGS)

CCSRCS	:= pm_seq.cc

OBJS	:= wtime.o readMatrixSeq.o $(CCSRCS:.cc=.o)
BINS	:= $(CCSRCS:.cc=)
LDLIBS	:= -L$(shell echo $(LIBRARY_PATH) | sed -e 's/:/\n/g' | grep pgi) -lhdf5

.PHONY:		all
all: $(BINS)

pm_seq:		pm_seq.cc wtime.o readMatrixSeq.o

wtime.o:		wtime.c wtime.h

readMatrixSeq.o:	readMatrixSeq.cc

.PHONY:		clean clobber
clean clobber:
	$(RM) $(BINS) $(OBJS)
