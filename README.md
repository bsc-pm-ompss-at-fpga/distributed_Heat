# Distributed Heat

**Name**: Heat Simulation  
**Contact Person**: OmpSs@FPGA Team, ompss-fpga-support@bsc.es  
**License Agreement**: GPL 3.0  
**Platform**: OmpSs@FPGA+IMP+OMPIF

## Description

The Heat benchmark simulates the propagation of heat on a 2D surfate following this formula on each position of the matrix:

![image](https://github.com/bsc-pm-ompss-at-fpga/distributed_Heat/assets/17345627/8fa75bd6-f204-464e-b3fa-79657ad8b48f)

$i,j$ are the row and column index respectively, and $k$ is the step number.
`A` is the matrix, with dimensions $(n+2) \times (m+2)$.
The padding rows and columns are used to set the initial heat points, and to avoid reading out of range.
Thus, the range of $i,j$ is $i \in [1, n]$ and $j \in [1, m]$.

<img src="https://github.com/bsc-pm-ompss-at-fpga/distributed_Heat/assets/17345627/ba3e6400-cbd5-46ff-86d8-cb7b1f541781" width="256" height="256" /> ![heat_animation](https://github.com/bsc-pm-ompss-at-fpga/distributed_Heat/assets/17345627/e0770ef6-e25c-4623-9933-7cd53b3c7c0d)

The left figure shows the initial matrix with the padding rows and columns.
Blue color means temperature 0, and the more read, the higher the temperature is.
Thus, initially the matrix is set to 0, and the padding rows and columns have the inital values with temperatures different from 0.
In the right image we see an animation of the heat propagation over time.
Every frame shows the matrix in a different timestep.
Both images are at different resolutions to improve visualization.
The left one is 32x32, while the animation is set at 256x256.

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

## FPGA implementation

Although the formula to update the matrix is simple, you will see that file `computeBlock.cpp` is far from being understandable, at least on the first read.
The code just implements the update of a single 256x256 block.
However, in order to get some performance we have to get the code a bit messy.
Everything important is in function `computeBlock_moved`.
You will see this variable declaration
``` C++
	constexpr int ITER_LAT = 76;
	constexpr int ROWPADDING = 1;
	constexpr int LPADDING = 64 / sizeof(double);
	constexpr int RPADDING = 1;
	double local[BS + ROWPADDING*2][BS + LPADDING + RPADDING];
````

`local` stores the block in BRAM, which is must faster to access and lets you configure the amount of ports for parallel reads/writes.
`BS` is the block size, in this case 256.
However, `local` has extra rows and columns, defined by the `ROWPADDING`, `LPADDING` and `RPADDING` constants.
The padding is not symmetric due to how the code accesses memory.
A common technique to increase memory bandwidth is to use a high bit-width data bus to read and write memory, so we can read/write multiple data elements in parallel.
For example, in the Heat case we use a port of 256-bit (the width of the Alveo U55C HBM), which is equivalent to loading 4 doubles.
However, to do that we access memory through an pointer to an `ap_uint<256>` type, therefore the addresses must be aligned to 32 bytes.
This is why we add a padding of 64 bytes on the left and right columns of the matrix, so we can support memory bus widths up to 512 bits (the usual case of FPGA DDR4 memories).
For the `local` array, we only need to add 64 bytes of padding in the left column, because we read by rows and there is no alignment constraint for BRAMs.
Both matrix and `local` are stored by consecutive rows, but the row size of the matrix (set by the user at runtime) is bigger than `local` (256).
This prevents us from doing a single consecutive copy when loading a block in `local`.
The next image shows the memory layout of both `local` and the input matrix.

<p align=center>
    <img src="https://github.com/bsc-pm-ompss-at-fpga/distributed_Heat/assets/17345627/67f7e82b-5cf8-4520-8b37-b2b15f5162be" width=514 height=250>
</p>

The left figure is the block stored in `local`, while the right figure is the whole matrix.
The orange rows are the contents of the block that will be updated by the task.
The green rows are part of the padding of the matrix, while the blue rows are part of the padding of `local`, but are not in the padding of the matrix, thus they are updated by another task.
In this example, the code will do 9 loads, one for each row including padding, and after updating the block, it will store 7 rows without padding.

Now lets see what the `ITER_LAT` is about.
Up until now we have seen how a block is loaded from memory to BRAMs in `local`, and stored back to memory.
The update part is divided in three loops: `upper_triangle_no_collapse`, `middle_collapse`, `lower_triangle_no_collapse`.
The main idea to get performance is to pipeline the loop that updates positions, and unroll it if possible so we can update several positions every cycle.
However, we have seen in the description section that there are many dependencies between neighbor elements, so a simple double loop won't work.
To simplify the code, the following examples do not use padding and use square matrices (which is true in our case because the block size is always square), so the $i,j$ range is from 0 to $n-1$.
```C++
for (int i = 0; i < n; ++i)
  for (int j = 0; j < n; ++j)
    A[i][j] = ...
```
If we try to pipeline the innermost loop, we will get a carried dependency warning, because the results of one iteration are needed by the next one, preventing the loop from reaching II=1.
For that, we need to traverse the block in anti-diagonal order, since all elements in the same anti-diagonal can be updated in parallel:
```C++
upper_triangle:
for (int d = 1; d <= n; ++d)
  for (int i = d-1, j = 0; i >= 0 || j < d; --i, ++j)
    A[i][j] = ...
lower_triangle:
for (int d = 1; d < n; ++d)
  for (int i = n-1, j = d; i >= d || j < n; --i, ++j)
    A[i][j] = ...
```
![heat_animation (5)](https://github.com/bsc-pm-ompss-at-fpga/distributed_Heat/assets/17345627/d7141507-273b-45a4-9f46-11ff72f25488)

This image shows how the upper (red) and lower (green) triangles formed by the anti-diagonals.
The numbers show the anti-diagonal ($d$, starting at 1) each matrix element belongs to in their respective triangle.

The first loop updates all the anti-diagonals for the upper triangle (including the main anti-diagnoal when $d == n$), and the second loop updates the remaining anti-diagonals of the lower triangle.
There are two conditions in both innsermost loops that are completetely equivalent, I show them for completeness but only one is enough.
In fact, both loops can be merged, and the result is:
```C++
for (int d = 1; d < 2*n; ++d)
  for (int i = d <= n ? d-1 : n-1,
           j = d <= n ? 0   : d-n;
       i >= 0 && j < n; --i, ++j)
    A[i][j] = ...
```

Now we can pipeline the innermost loop `ij` with II=1, however this loop can't get flattened because it's not a trivial loop (the iteration count is not constant).
Thus, between two iterations of the outermost `d` loop we have to pay the latency of the `ij` loop.
We can flatten manually both loops into one, however there is a dependency between iterations when `d` is different.
The loop has a latency $L$, which means that on a certain cycle $t$, the value of $i,j$ in that iteration is updated at cycle $t+L$, assuming II=1.

On a certain diagonal $d$ of the upper triangle, any $i,j$ iteration starting at cycle $t$ has to finish by cycle $t+d$, because in this cycle the $i+1,j$ iteration starts in diagonal $d+1$, which needs the updated $i,j$ value.
Therefore, if we pipeline the collpased loop, we have to ensure that $L <= d$ for all values of $d$ produced by the loop.
Similarly, in the lower triangle, any $i,j$ iteration, except when $i == n-1$, of diagonal $d$ starting at cycle $t$ needs to finish before the starting cycle of $i+1,j$ in diagonal $d+1$, which happens in cycle $t+n-d-1$.
When $i == n-1$, the dependency is with iteration $i,j+1$ which happens at cycle $t+n-d$.
Since the first condition is more restrictive, to avoid dependencies we need $L < d$.
If we want to use a single loop for both triangles, then we have to use the condition $L < d$.

![heat_animation (4)](https://github.com/bsc-pm-ompss-at-fpga/distributed_Heat/assets/17345627/81e4ebc2-68be-45e2-a345-e922d5291197)

This image shows the starting cycle of each iteration for diagonal 3 of the upper triangle and 2 of the lower triangle.
We can see that for the first case, at cycle 3 we need the value of $(0,1)$ calculated, while in the second case we need the value of $(3,3)$ by cycle 2.

The condition $L < d$ determines which anti-diagonals are calculated in a single pipelined, collapsed loop.
The trailing anti-diagonals are calculated with the double loop without collapse.
The value of $L$ is stored in `ITER_LAT` in the code, this is why there are three loops.
In summary, the following image shows which parts of the block are updated by each loop.

![heat_animation (2)](https://github.com/bsc-pm-ompss-at-fpga/distributed_Heat/assets/17345627/6ae51381-4129-41e4-9eb8-0f6bb611ff0c)

However, there is one last plot twist.
We can in fact unroll manually the loops to get double performance.
This is a tradeoff, because if we calculate two elements per cycle, the amount of anti-diagonals we can calculate in the collapsed loop gets reduced, and `ITER_LAT=2L`.
It depends on the actual latency $L$, which is reported by Vitis HLS, and the block size.
You can actually calculate if it is better to unroll, or even find the sweet spot for the block size, but the formulas are too big to put here.

## How to compile

Sadly, there is no official support in the clang compiler for OMPIF and IMP, so we have to split manually the FPGA and the host part. The host code is under `src`. You will see a `solver.c` and `solver_dummy.c`.
The first is just the sequential version of the Heat benchmark, and the second is used to generate the executable, since the compiler needs the definition of a task even if it is not a CPU task.
In execution time, the code will search for the FPGA accelerators, which are implemented in hls under the `hls` directory.

So, first you need an OmpSs@FPGA installation with broadcaster support for both in Nanos6 and clang. Also, you need the AIT version that implements Ethernet subsystem and OMPIF support. For the moment these versions are not official nor public, you can ask for access to ompss-fpga-support@bsc.es.

In summary, to compile the host executable run `make`.
To generate the bitstream run `make ait`
By default, it generates the bitstream from the HLS phase, but you can change the start and finish phases with the Makefile variables `TO_STEP` and `FROM_STEP`.
Besides that, there are other variables:
* FPGA_CLOCK: Frequency in MHz at which the accelerators will run.
* FPGA_MEMPORT_WIDTH: Data bit-width of the memory port for all the accelerators. More bit-width may provide more bandwidth (depending on the FPGA memory path), at the cost of more resource usage. **IMPORTANT** This variable is intended to be used in the task pragma. Since there is no compiler support, if you want to change the port width you have to modify the hls code `computeBlock.cpp`. By default it is 256.
* BLOCK_SIZE: Size of one of the dimensions of the block matrices, which are square. **IMPORTANT** The block size affects both the host and hls code, so if you modify the variable, you have to apply the changes manually in both `solve.cpp` and `computeBlock.cpp`.
* NUM_CBLOCKS: Number of accelerators that compute a single block. **IMPORTANT** If you want to change this variable, change the `num_instances` field of the `ait_extracted.json` file. Usually this file is generated by clang, but since we do not have that support with OMPIF or IMP, we have to modify the file manually.
