# distributed_Heat

**Name**: Heat Simulation  
**Contact Person**: OmpSs@FPGA Team, ompss-fpga-support@bsc.es  
**License Agreement**: GPL 3.0  
**Platform**: OmpSs@FPGA+IMP+OMPIF

## Description

The Heat benchmark simulates the propagation of heat on a 2D surfate following this formula on each position of the matrix:

![image](https://github.com/bsc-pm-ompss-at-fpga/distributed_Heat/assets/17345627/8fa75bd6-f204-464e-b3fa-79657ad8b48f)

$i,j$ are the row and column index respectively, and $k$ is the step number.
This benchmark uses a Gauss-Seidel method, which means the same matrix is used as input and output.
Therefore, in order to parallelize this applications we have to look at the dependencies the sequential version creates.
As the formula tells, each position is calculated with the value of the previous step for the right and down neighbors, and the value of the current step of the up and left neighbors.

This creates a pattern of dependencies that acts as a wave.
First, only the top-leftmost element can be updated, then, its right and down neighbors, and so on.
The following video shows this effect with several steps overlapping:

![test](https://github.com/bsc-pm-ompss-at-fpga/distributed_Heat/assets/17345627/721da55f-a0d1-4c5d-9918-cd0a3052ca22)

The numbers show the anti-diagonal identifier.
All elements of the same anti-diangonal can always be updated in parallel.
The colors represent the elements that are updated for each step.
Therefore, at the same time, only the even or odd anti-diagonals can be updated in parallel, and each anti-diagonal updates the values of a different step.

## Parallelization with tasks

There is only one kernel, which updates the position sequentially in a fixed-size block of the matrix.
Therefore, the same wave pattern applies to the tasks: a task depends on the four tasks that update the neighbor blocks.


## Parallelization with Implicit Message Passing

This one is very similar to classic MPI implementations of the Gauss-Seidel Heat.
The idea is to distribute the matrix evenly accross all ranks by consecutive rows.
Then, we add a single extra row at the top and bottom to hold the values calculated by the other ranks.
Therefore, before executing a task, the rank may need a row with the values that of the current or previous iteration.
If these are calculated by other rank, the there is a data movement using OMPIF_Send/Recv.
We can see this better in the following image:

![heat_animation (1)](https://github.com/bsc-pm-ompss-at-fpga/distributed_Heat/assets/17345627/485e3584-91eb-4058-a6e4-5be848b67ec0)

Here we see how we split the matrix in three segments and distribute each tile to a different rank.
Then, before starting the current step, every rank needs rows from their bottom neighbor.
Each OMPIF_Send/Recv only sends the row of a block, therefore there are as many tasks as blocks in a single row.
After the data exchange, the top rank can start, and as soon as the blocks are updated, it can send the updated rows to the bottom neighbors.

