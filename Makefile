# Makefile for the popbam program
ifndef CC
CC=               gcc
endif
CFLAGS=           -D_FILE_OFFSET_BITS=64
ifndef CXX
CXX=              g++
endif
CXXFLAGS=         -D_FILE_OFFSET_BITS=64 -std=c++0x
C_RELEASE_FLAGS=  -O2
C_DEBUG_FLAGS=    -Wall -ggdb -DDEBUG
C_PROFILE_FLAGS=  -Wall -O2 -pg
CSOURCES=          bam_aux.c bam.c sam_header.c bam_import.c sam.c \
		   bam_plbuf.c
CXXSOURCES=        popbam.cpp pop_utils.cpp pop_sample.cpp pop_tree.cpp \
                   pop_snp.cpp pop_nucdiv.cpp pop_ld.cpp pop_sfs.cpp \
                   pop_diverge.cpp pop_haplo.cpp getopt.cpp pop_options.cpp \
                   pop_fasta.cpp
OBJS=              popbam.o pop_utils.o pop_sample.o pop_tree.o \
                   pop_snp.o pop_nucdiv.o pop_ld.o pop_sfs.o pop_diverge.o \
                   pop_haplo.o getopt_pp.o pop_options.o pop_fasta.o \
                   bam_aux.o bam.o sam_header.o bam_import.o sam.o \
                   bam_plbuf.o
PROG=              popbam
LIBFLAGS=          -lhts -lz -lm
INSTALL_DIR=       /usr/local/bin

all: CXXFLAGS +=  $(C_RELEASE_FLAGS)
all: CFLAGS += $(C_RELEASE_FLAGS)
all: $(PROG)

.PHONY: debug
debug: CXXFLAGS += $(C_DEBUG_FLAGS)
debug: CFLAGS += $(C_DEBUG_FLAGS)
debug: $(PROG)

.PHONY: release
release: CXXFLAGS += $(C_RELEASE_FLAGS)
release: CFLAGS += $(C_RELEASE_FLAGS)
release: $(PROG)

.PHONY: profile
profile: CXXFLAGS += $(C_PROFILE_FLAGS)
profile: CFLAGS += $(C_PROFILE_FLAGS)
profile: $(PROG)

$(PROG): $(OBJS)
	$(CXX) -o $@ $(OBJS) $(LIBFLAGS)

%.o : %.c
	$(CC) $(CFLAGS) -c -o $@ $<

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

.PHONY: install
install: $(PROG)
	install -m 755 $(PROG) $(INSTALL_DIR)

.PHONY: clean
clean:
	rm -rf $(PROG) $(OBJS)

