# Compilers
CC = clang

BLOCK_SIZE ?= 128
NUM_CBLOCKS ?= 1

# Compiler flags
PARAM_DEFINES = -DBLOCK_SIZE=$(BLOCK_SIZE) -DNUM_CBLOCKS=$(NUM_CBLOCKS)
CFLAGS=-fompss-2 -Ofast -ffast-math $(PARAM_DEFINES)

# Kernel flags
KFLAGS=$(CFLAGS)

FPGA_CLOCK=100
FPGA_HWRUNTIME=pom
FPGA_MEMPORT_WIDTH=512
PROGRAM_=heat$(AIT_SUFFIX)
FROM_STEP ?= HLS
TO_STEP ?= design

OBJ=misc.o solver_dum.o

all: heat

misc.o: src/misc.c
	$(CC) $(CFLAGS) -c -o $@ $^

solver_dum.o: src/solver_dummy.c
	$(CC) $(CFLAGS) -c -o $@ $^

heat: src/main.c $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^
