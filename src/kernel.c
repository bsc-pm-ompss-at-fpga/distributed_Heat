#include <stdint.h>
#include "heat.h"

#include <nanos6/fpga_device.h>

void computeBlock(const int cols,
		unsigned long long int M)
{
#pragma HLS inline
	//The real latency is 37 cycles, but the loops are unrolled by a factor of two.
	//A diagonal of 73 elements uses 73/2=36.5 iterations which take 37 cycles to complete.
	const int ITER_LAT = 73;
	const int B_ROWPADDING = 1;
	const int B_LPADDING = LPADDING;
	const int B_RPADDING = 3;

	double local[BS+B_ROWPADDING*2][BS+B_LPADDING+B_RPADDING];
#pragma HLS array_partition variable=local cyclic factor=4 dim=0

	load_local:
	for (int i = 0; i < BS+B_ROWPADDING*2; ++i) {
		nanos6_fpga_memcpy_wideport_in(&local[i][0], M + (i*cols)*sizeof(double), BS+B_LPADDING+B_RPADDING);
	}

	//fl: first left, fd: first down, c1: common 1, c2: common 2, su: second up, sr: second right
	upper_triangle_no_collapse:
	for (int d = 1; d < ITER_LAT; ++d) {
		int r = (d-1)-1 + B_ROWPADDING;
		int c = B_LPADDING;
		double fl = local[r+1][c-1];
		double fd = local[r+2][c];
		for (; c < d + B_LPADDING; r -= 2, c+=2) {
#pragma HLS pipeline II=1
#pragma HLS dependence variable=local inter false
			double c1 = local[r][c];
			double c2 = local[r+1][c+1];
			double su = local[r-1][c+1];
			double sr = local[r][c+2];
			local[r+1][c] = 0.25 * (fl + fd + c1 + c2);
			if (c+1 < d + B_LPADDING)
				local[r][c+1] = 0.25 * (c1 + c2 + su + sr);
			fl = su; fd = sr;
		}
	}

	{
		int r, c;
		int d = ITER_LAT;
		r = d - 2 + B_ROWPADDING;
		c = B_LPADDING;
		upper_triangle_pipe:
		while (d <= BS) {
#pragma HLS pipeline II=1
#pragma HLS dependence variable=local inter false
			double fl = local[r+1][c-1];
			double fd = local[r+2][c];
			double c1 = local[r][c];
			double c2 = local[r+1][c+1];
			double su = local[r-1][c+1];
			double sr = local[r][c+2];
			local[r+1][c] = 0.25 * (fl + fd + c1 + c2);
			if (c+1 < d+B_LPADDING)
				local[r][c+1] = 0.25 * (c1 + c2 + su + sr);
			r -= 2, c += 2;
			if (c >= d + B_LPADDING)
			{
				++d;
				c = B_LPADDING;
				r = d - 2 + B_ROWPADDING;
			}
		}
		d = 1;
		r = (BS-1)-1 + B_ROWPADDING;
		c = d + B_LPADDING;
		lower_triangle_collapse:
		while (d <= BS - ITER_LAT)
		{
#pragma HLS pipeline II=1
#pragma HLS dependence variable=local inter false
			double fl = local[r+1][c-1];
			double fd = local[r+2][c];
			double c1 = local[r][c];
			double c2 = local[r+1][c+1];
			double su = local[r-1][c+1];
			double sr = local[r][c+2];
			local[r+1][c] = 0.25 * (fl + fd + c1 + c2);
			if (r-1 >= d-1 + B_ROWPADDING)
				local[r][c+1] = 0.25 * (c1 + c2 + su + sr);
			r -= 2, c += 2;
			if (r < d-1 + B_ROWPADDING)
			{
				++d;
				r = (BS-1)-1 + B_ROWPADDING;
				c = d + B_LPADDING;
			}
		}
	}

	lower_triangle_no_collapse:
	for (int d = BS - ITER_LAT + 1; d < BS;  ++d)
	{
		int r = (BS-1)-1 + B_ROWPADDING;
		int c = d + B_LPADDING;
		double fl = local[r+1][c-1];
		double fd = local[r+2][c];
		for (; r >= d-1 + B_ROWPADDING; (r -= 2,  c += 2))
		{
#pragma HLS pipeline II=1
#pragma HLS dependence variable=local inter false
			double c1 = local[r][c];
			double c2 = local[r+1][c+1];
			double su = local[r-1][c+1];
			double sr = local[r][c+2];
			local[r+1][c] = 0.25 * (fl + fd + c1 + c2);
			if (r-1 >= d-1 + B_ROWPADDING)
				local[r][c+1] = 0.25 * (c1 + c2 + su + sr);
			fl = su; fd = sr;
		}
	}

	//remeber to fix the number of iterations of the nanos6_fpga_memcpy_wideport
	store_local:
	for (int i = B_ROWPADDING; i < BS+B_ROWPADDING; ++i) {
		nanos6_fpga_memcpy_wideport_out(&local[i][B_LPADDING], M + (i*cols + B_LPADDING)*sizeof(double), BS);
	}
}

double computeBlockResidual(const int64_t rows, const int64_t cols,
		const int rstart, const int rend,
		const int cstart, const int cend,
		double M[rows][cols])
{
	double sum = 0.0;
	for (int r = rstart; r <= rend; ++r) {
		for (int c = cstart; c <= cend; ++c) {
			const double value = 0.25*(M[r-1][c] + M[r+1][c] + M[r][c-1] + M[r][c+1]);
			const double diff = value - M[r][c];
			sum += diff*diff;
			M[r][c] = value;
		}
	}
	return sum;
}
