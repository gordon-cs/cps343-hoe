CC	:= nvc
CXX	:= nvc++
CFLAGS	:= -fast
CXXFLAGS:= $(CFLAGS)

CCSRCS	:= pm_seq.cc

OBJS	:= $(CCSRCS:.cc=.o) readMatrixSeq.o wtime.o
BINS	:= $(CCSRCS:.cc=)
LDLIBS	:= -lhdf5

.PHONY:		all
all: $(BINS)

pm_seq:			pm_seq.o wtime.o readMatrixSeq.o
	$(LINK.cc) $^ $(LOADLIBS) $(LDLIBS) -o $@

pm_seq.o:		pm_seq.cc wtime.h

readMatrixSeq.o:	readMatrixSeq.cc

wtime.o:		wtime.c wtime.h

.PHONY:		clean clobber
clean clobber:
	$(RM) $(BINS) $(OBJS)
