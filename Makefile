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
CXXSOURCES=        pop_main.cpp pop_utils.cpp pop_sample.cpp pop_diverge.cpp \
                   pop_error.cpp
OBJS=              pop_main.o pop_utils.o pop_sample.o pop_diverge.o \
                   pop_error.o bam_aux.o bam.o sam_header.o bam_import.o \
                   sam.o bam_plbuf.o
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

