///////////////////
// Automatic IP Generated by OmpSs@FPGA compiler
///////////////////
// The below code is composed by:
//  1) User source code, which may be under any license (see in original source code)
//  2) OmpSs@FPGA toolchain code which is licensed under LGPLv3 terms and conditions
///////////////////
// Top IP Function: computeBlock
// Accel. type hash: 5050920982
// Num. instances: 1
// Wrapper version: 13
///////////////////

#include <hls_stream.h>
#include <ap_int.h>
#include <ap_axi_sdata.h>

static ap_uint<64> __mcxx_taskId;
template<class T>
union __mcxx_cast {
	unsigned long long int raw;
	T typed;
};
struct mcxx_inaxis {
	ap_uint<64> data;
};
typedef ap_axiu<64, 1, 1, 2> mcxx_outaxis;

template<class T>
void nanos6_fpga_memcpy_wideport_in(T * dst, const unsigned long long int addr, const unsigned int num_elems, ap_uint<256>* mcxx_memport) {
#pragma HLS inline
	for (int i = 0; i < (num_elems-1)/(sizeof(ap_uint<256>)/sizeof(T))+1; ++i) {
#pragma HLS pipeline II=1
		ap_uint<256> tmpBuffer;
		tmpBuffer = *(mcxx_memport + addr/sizeof(ap_uint<256>) + i);
		for (int j = 0; j < (sizeof(ap_uint<256>)/sizeof(T)); ++j) {
			if (i*(sizeof(ap_uint<256>)/sizeof(T))+j >= num_elems)
				continue;
			__mcxx_cast<T> cast_tmp;
			cast_tmp.raw = tmpBuffer((j+1)*sizeof(T)*8-1, j*sizeof(T)*8);
			dst[i*(sizeof(ap_uint<256>)/sizeof(T))+j] = cast_tmp.typed;
		}
	}
}
template<class T>
void nanos6_fpga_memcpy_wideport_out(T * src, const unsigned long long int addr, const unsigned int num_elems, ap_uint<256>* mcxx_memport) {
#pragma HLS inline
	for (int i = 0; i < (num_elems-1)/(sizeof(ap_uint<256>)/sizeof(T))+1; ++i) {
#pragma HLS pipeline II=1
		ap_uint<256> tmpBuffer;
		for (int j = 0; j < (sizeof(ap_uint<256>)/sizeof(T)); ++j) {
			__mcxx_cast<T> cast_tmp;
			cast_tmp.typed = src[i*(sizeof(ap_uint<256>)/sizeof(T))+j];
			tmpBuffer((j+1)*sizeof(T)*8-1, j*sizeof(T)*8) = cast_tmp.raw;
		}
		*(mcxx_memport + addr/sizeof(ap_uint<256>) + i) = tmpBuffer;
	}
}

static constexpr int BS = 256;
void computeBlock_moved(const int cols, unsigned long long int M, ap_uint<256> *mcxx_memport)
{
#pragma HLS inline
	constexpr int ITER_LAT = 76;
	constexpr int ROWPADDING = 1;
	constexpr int LPADDING = 64 / sizeof(double);
	constexpr int RPADDING = 1;
	double local[BS + ROWPADDING*2][BS + LPADDING + RPADDING];

#pragma HLS array_partition variable=local cyclic factor=4 dim=0
	load_local : for (int i = 0; i < BS + ROWPADDING * 2;  ++i)
	{
		nanos6_fpga_memcpy_wideport_in(&local[i][0], M + i * cols * sizeof(double), BS + LPADDING + RPADDING, mcxx_memport);
	}
	upper_triangle_no_collapse : for (int d = 1; d < ITER_LAT;  ++d)
	{
		int r = d - 1 - 1 + ROWPADDING;
		int c = LPADDING;
		double fl = local[r + 1][c - 1];
		double fd = local[r + 2][c];
		for (; c < d + LPADDING; (r -= 2, c += 2))
		{
#pragma HLS pipeline II=1
#pragma HLS dependence variable=local inter false
			double c1 = local[r][c];
			double c2 = local[r + 1][c + 1];
			double su = local[r - 1][c + 1];
			double sr = local[r][c + 2];
			local[r + 1][c] = 2.50000000000000000000000000000000000000000000000000000e-01 * (fl + fd + c1 + c2);
			if (c + 1 < d + LPADDING)
			{
				local[r][c + 1] = 2.50000000000000000000000000000000000000000000000000000e-01 * (c1 + c2 + su + sr);
			}
			fl = su;
			fd = sr;
		}
	}
	{
		int r;
		int c;
		int d = ITER_LAT;
		r = d - 2 + ROWPADDING;
		c = LPADDING;
		middle_collapse : while (d <= 2*BS - ITER_LAT)
		{
#pragma HLS pipeline II=1
#pragma HLS dependence variable=local inter false
			const double fl = local[r + 1][c - 1];
			const double fd = local[r + 2][c];
			const double c1 = local[r][c];
			const double c2 = local[r + 1][c + 1];
			const double su = local[r - 1][c + 1];
			const double sr = local[r][c + 2];
			local[r + 1][c] = 2.50000000000000000000000000000000000000000000000000000e-01 * (fl + fd + c1 + c2);
			if (d <= BS ? c + 1 < d + LPADDING : r - 1 >= (d-BS) - 1 + ROWPADDING)
			{
				local[r][c + 1] = 2.50000000000000000000000000000000000000000000000000000e-01 * (c1 + c2 + su + sr);
			}
			r -= 2; c += 2;
			if (d <= BS ? c >= d + LPADDING : r < (d-BS) - 1 + ROWPADDING)
			{
				++d;
				r = d <= BS ? d - 2 + ROWPADDING : BS - 1 - 1 + ROWPADDING;
				c = d <= BS ? LPADDING : (d-BS) + LPADDING;
			}
		}
	}
	lower_triangle_no_collapse : for (int d = BS - ITER_LAT + 1; d < BS;  ++d)
	{
		int r = BS - 1 - 1 + ROWPADDING;
		int c = d + LPADDING;
		double fl = local[r + 1][c - 1];
		double fd = local[r + 2][c];
		for (; r >= d - 1 + ROWPADDING; (r -= 2, c += 2))
		{
#pragma HLS pipeline II=1
#pragma HLS dependence variable=local inter false
			double c1 = local[r][c];
			double c2 = local[r + 1][c + 1];
			double su = local[r - 1][c + 1];
			double sr = local[r][c + 2];
			local[r + 1][c] = 2.50000000000000000000000000000000000000000000000000000e-01 * (fl + fd + c1 + c2);
			if (r - 1 >= d - 1 + ROWPADDING)
			{
				local[r][c + 1] = 2.50000000000000000000000000000000000000000000000000000e-01 * (c1 + c2 + su + sr);
			}
			fl = su;
			fd = sr;
		}
	}
	store_local : for (int i = ROWPADDING; i < BS + ROWPADDING;  ++i)
	{
		nanos6_fpga_memcpy_wideport_out(&local[i][LPADDING], M + (i * cols + LPADDING) * sizeof(double), BS, mcxx_memport);
	}
}

void mcxx_write_out_port(const ap_uint<64> data, const ap_uint<3> dest, const ap_uint<1> last, hls::stream<mcxx_outaxis>& mcxx_outPort) {
#pragma HLS inline
	mcxx_outaxis axis_word;
	axis_word.data = data;
	axis_word.dest = dest;
	axis_word.last = last;
	mcxx_outPort.write(axis_word);
}

void computeBlock_wrapper(hls::stream<ap_uint<64> >& mcxx_inPort, hls::stream<mcxx_outaxis>& mcxx_outPort, ap_uint<256>* mcxx_memport) {
#pragma HLS interface ap_ctrl_none port=return
#pragma HLS interface axis port=mcxx_inPort
#pragma HLS interface axis port=mcxx_outPort
#pragma HLS interface m_axi port=mcxx_memport
	mcxx_inPort.read(); //command word
	__mcxx_taskId = mcxx_inPort.read();
	ap_uint<64> __mcxx_parent_taskId = mcxx_inPort.read();
	int cols;
	unsigned long long int M;
	{
#pragma HLS protocol fixed
		{
			ap_uint<8> mcxx_flags_0;
			ap_uint<64> mcxx_offset_0;
			mcxx_flags_0 = mcxx_inPort.read()(7,0);
			ap_wait();
			__mcxx_cast<int> mcxx_arg_0;
			mcxx_arg_0.raw = mcxx_inPort.read();
			cols = mcxx_arg_0.typed;
		}
		ap_wait();
		{
			ap_uint<8> mcxx_flags_1;
			ap_uint<64> mcxx_offset_1;
			mcxx_flags_1 = mcxx_inPort.read()(7,0);
			ap_wait();
			__mcxx_cast<unsigned long long int> mcxx_arg_1;
			mcxx_arg_1.raw = mcxx_inPort.read();
			M = mcxx_arg_1.typed;
		}
		ap_wait();
	}
	computeBlock_moved(cols, M, mcxx_memport);
	{
#pragma HLS protocol fixed
		ap_uint<64> header = 0x03;
		ap_wait();
		mcxx_write_out_port(header, 0, 0, mcxx_outPort);
		ap_wait();
		mcxx_write_out_port(__mcxx_taskId, 0, 0, mcxx_outPort);
		ap_wait();
		mcxx_write_out_port(__mcxx_parent_taskId, 0, 1, mcxx_outPort);
		ap_wait();
	}
}
