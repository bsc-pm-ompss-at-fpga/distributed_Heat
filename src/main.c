#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "heat.h"

#include <nanos6/debug.h>
#include <nanos6/distributed.h>

int main(int argc, char **argv)
{
	HeatConfiguration conf;
	int ret = readConfiguration(argc, argv, &conf);
	if (ret == 2) {
		return 0;
	}
	else if (ret) {
		return 1;
	}

	int nranks = nanos6_dist_num_devices();
	if (nranks <= 0) {
		fprintf(stderr, "No devices found!\n");
		return 1;
	}

	if (conf.rows%BLOCK_SIZE != 0 || conf.cols%BLOCK_SIZE != 0 || (conf.rows/BLOCK_SIZE)%nranks != 0) {
		fprintf(stderr, "Incompatible configuration\n");
		return 1;
	}

	refineConfiguration(&conf, nranks*conf.rbs, conf.cbs);
	if (conf.verbose) printConfiguration(&conf);

	int64_t rows = conf.rows+VPADDING;
	int64_t cols = conf.cols+HPADDING;

	int err = initialize(&conf, rows, cols);
	assert(!err);

	const int nrb = conf.rows/conf.rbs+2;
	const int ncb = conf.cols/conf.cbs+2;
	char* representatives = malloc(nrb*ncb);

	size_t rowsPerDevice = conf.rows/nranks;

	nanos6_dist_map_address(conf.matrix, (rowsPerDevice+VPADDING)*cols*sizeof(double));
	nanos6_dist_map_address(representatives, nrb*ncb);

	double start, end;
	start = getTime();
	nanos6_dist_memcpy_to_device(0, conf.matrix, cols*sizeof(double), (UPADDING-1)*cols*sizeof(double), (UPADDING-1)*cols*sizeof(double));
	size_t* sizes = (size_t*)malloc(nranks*sizeof(size_t));
	size_t* srcOffsets = (size_t*)malloc(nranks*sizeof(size_t));
	size_t* dstOffsets = (size_t*)malloc(nranks*sizeof(size_t));
	for (int i = 0; i < nranks; ++i) {
		sizes[i] = rowsPerDevice*cols*sizeof(double);
		srcOffsets[i] = (UPADDING+rowsPerDevice*i)*cols*sizeof(double);
		dstOffsets[i] = UPADDING*cols*sizeof(double);
	}
	nanos6_dist_scatterv(conf.matrix, sizes, srcOffsets, dstOffsets);
	free(sizes);
	free(srcOffsets);
	free(dstOffsets);
	nanos6_dist_memcpy_to_device(nranks-1, conf.matrix, cols*sizeof(double), (rows-DPADDING)*cols*sizeof(double), (UPADDING+rowsPerDevice)*cols*sizeof(double));
	end = getTime();
	double bandwidth = ((conf.rows+2)*cols*sizeof(double))/(end-start);
	fprintf(stderr, "Copy bandwidth %gMB/s\n", bandwidth/1e6);

	// Solve the problem
	start = getTime();
	solve(conf.matrix, representatives, conf.rbs, conf.cbs, rows, cols, conf.timesteps);
	#pragma oss taskwait
	end = getTime();

	int64_t totalElements = conf.rows*conf.cols;
	double throughput = (totalElements*conf.timesteps)/(end-start);
	throughput = throughput/1000000.0;

	int threads = nanos6_get_num_cpus();

	if (conf.parse) {
		printf("%e\n", end - start);
	}
	else {
		fprintf(stdout, "rows, %ld, cols, %ld, rows/rank, %ld, total, %ld, total/rank, %ld, rbs, %d, "
				"cbs, %d, ranks, %d, threads, %d, timesteps, %d, time, %g, Mupdates/s, %f\n",
				conf.rows, conf.cols, conf.rows/nranks, totalElements, totalElements/nranks,
				conf.rbs, conf.cbs, nranks, threads, conf.timesteps, end-start, throughput);
	}

	if (conf.generateImage || conf.compareRef || conf.createRef) {
		for (int i = 0; i < nranks; ++i) {
			nanos6_dist_memcpy_from_device(i, conf.matrix, rowsPerDevice*cols*sizeof(double), UPADDING*cols*sizeof(double), (UPADDING+rowsPerDevice*i)*cols*sizeof(double));
		}
	}

	nanos6_dist_unmap_address(conf.matrix);
	nanos6_dist_unmap_address(representatives);
	free(representatives);

	err = 0;
	if (conf.compareRef) {
		err |= compareOutputs(conf.refFileName, conf.rows, conf.cols, LPADDING, RPADDING, UPADDING, conf.matrix);
	}
	if (conf.createRef) {
		err |= writeOutput(conf.outFileName, conf.rows, conf.cols, LPADDING, RPADDING, UPADDING, conf.matrix);
	}
	if (conf.generateImage) {
		err |= writeImage(conf.imageFileName, conf.matrix, rows, cols);
	}

	err |= finalize(&conf);

	return err;
}

