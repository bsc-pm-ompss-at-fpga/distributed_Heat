# Compilers
CC = clang

BLOCK_SIZE ?= 128
NUM_CBLOCKS ?= 1

# Compiler flags
PARAM_DEFINES = -DBLOCK_SIZE=$(BLOCK_SIZE) -DNUM_CBLOCKS=$(NUM_CBLOCKS)
CFLAGS=-fompss-2 -Ofast -ffast-math $(PARAM_DEFINES)

FPGA_CLOCK=100
FPGA_MEMPORT_WIDTH=512
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

ait:
	ait -b alveo_u55c -c $(FPGA_CLOCK) -n heat -v --disable_board_support_check --wrapper_version 13 --disable_spawn_queues --placement_file u55c_placement_$(NUM_CBLOCKS).json --floorplanning_constr all --slr_slices all --regslice_pipeline_stages 1:1:1 --interconnect_regslice all --enable_pom_axilite --max_deps_per_task=7 --max_args_per_task=4 --max_copies_per_task=2 --picos_tm_size=32 --picos_dm_size=162 --picos_vm_size=162 --from_step $(FROM_STEP) --to_step $(TO_STEP)
