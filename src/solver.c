#include <nanos6/distributed.h>

#include "heat.h"

static inline void gaussSeidelSolver(int64_t rows, int64_t cols, int rbs, int cbs, int nrb, int ncb, double* M, char* reps)
{
#pragma HLS inline
	for (int R = 1; R < nrb-1; ++R) {
		for (int C = 1; C < ncb-1; ++C) {
			#pragma oss taskcall \
					in(reps[(R-1)*ncb + C]) in(reps[(R+1)*ncb + C]) \
					in(reps[R*ncb + C-1]) in(reps[R*ncb + C+1]) \
					inout(reps[R*ncb + C])
			computeBlock(cols, (unsigned long long int)(M + ((R-1)*rbs + UPADDING - 1)*cols + (C-1)*cbs));
		}
	}
}

void solve(double *matrix, char* reps, const int rbs, const int cbs, int rows, int cols, int timesteps)
{
#pragma HLS inline
	const int nrb = (rows-VPADDING)/rbs+2;
	const int ncb = (cols-HPADDING)/cbs+2;

	for (int t = 0; t < timesteps; ++t) {
		gaussSeidelSolver(rows, cols, rbs, cbs, nrb, ncb, matrix, reps);
	}
	#pragma oss taskwait
}

