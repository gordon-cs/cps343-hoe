CXX	:= mpic++
CXXFLAGS:= -Ofast -DSHOW_TIMING_DATA
LDLIBS	:=

# Project files

SRCS	:= psrs_qsort.cc psrs_radixsort.cc

BINS 	:= $(SRCS:.cc=)
HDRS	:=
OBJS	:=

# ============================================================================
# ============================== TARGET SECTION ==============================
# ============================================================================

# Program targets

.PHONY:	all
all:	$(BINS)

# Supporting targets

.PHONY: clean
clean:
	$(RM) $(OBJS)

.PHONY: clobber
clobber:
	$(MAKE) clean
	$(RM) $(BINS)
