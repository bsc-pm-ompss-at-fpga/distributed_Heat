///////////////////
// Automatic IP Generated by OmpSs@FPGA compiler
///////////////////
// The below code is composed by:
//  1) User source code, which may be under any license (see in original source code)
//  2) OmpSs@FPGA toolchain code which is licensed under LGPLv3 terms and conditions
///////////////////
// Top IP Function: solve
// Accel. type hash: 8150371898
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
struct __fpga_copyinfo_t {
	unsigned long long int copy_address;
	unsigned char arg_idx;
	unsigned char flags;
	unsigned int size;
};
struct __data_owner_info_t {
	unsigned long long int size;
	unsigned char owner;
};
void mcxx_task_create(const ap_uint<64> type, const ap_uint<8> instanceNum, const ap_uint<8> numArgs, const unsigned long long int args[], const ap_uint<8> numDeps, const unsigned long long int deps[], const ap_uint<8> numCopies, const __fpga_copyinfo_t copies[], int numDataOwners, __data_owner_info_t data_owners[], hls::stream<mcxx_outaxis>& mcxx_outPort, unsigned char ompif_rank, unsigned char ompif_size, unsigned char owner);
void mcxx_task_create(const ap_uint<64> type, const ap_uint<8> instanceNum, const ap_uint<8> numArgs, const unsigned long long int args[], const ap_uint<8> numDeps, const unsigned long long int deps[], const ap_uint<8> numCopies, const __fpga_copyinfo_t copies[], hls::stream<mcxx_outaxis>& mcxx_outPort);
void mcxx_taskwait(hls::stream<ap_uint<8> >& mcxx_spawnInPort, hls::stream<mcxx_outaxis>& mcxx_outPort);
template <typename T>
struct __mcxx_ptr_t {
	unsigned long long int val;
	__mcxx_ptr_t(unsigned long long int val) : val(val) {}
	__mcxx_ptr_t() {}
	inline operator __mcxx_ptr_t<const T>() const {
		return __mcxx_ptr_t<const T>(val);
	}
	template <typename V> inline __mcxx_ptr_t<T> operator+(V const val) const {
		return __mcxx_ptr_t<T>(this->val + val*sizeof(T));
	}
	template <typename V> inline __mcxx_ptr_t<T> operator-(V const val) const {
		return __mcxx_ptr_t<T>(this->val - val*sizeof(T));
	}
	template <typename V> inline operator V() const {
		return (V)val;
	}
};
typedef enum {
	OMPIF_INT = 0,
	OMPIF_DOUBLE = 1,
	OMPIF_FLOAT = 2
} OMPIF_Datatype;
typedef enum {
	OMPIF_COMM_WORLD
} OMPIF_Comm;
void OMPIF_Send(const void *data, unsigned int size, int destination, const ap_uint<8> numDeps, const unsigned long long int deps[], hls::stream<mcxx_outaxis>& mcxx_outPort);
void OMPIF_Recv(void *data, unsigned int size, int source, const ap_uint<8> numDeps, const unsigned long long int deps[], hls::stream<mcxx_outaxis>& mcxx_outPort);

unsigned char calc_owner(unsigned int R, unsigned int nrb, unsigned int nranks, unsigned int blocks_per_rank) {
#pragma HLS inline
	if (R == 0) {
		return 0;
	}
	else if (R == nrb-1) {
		return nranks-1;
	}
	else {
		return (R-1)/blocks_per_rank;
	}
}

#define LPADDING (64/sizeof(double))
#define RPADDING LPADDING
#define HPADDING (LPADDING+RPADDING)
#define UPADDING 1
#define DPADDING 1
#define VPADDING (UPADDING+DPADDING)

constexpr unsigned int BS = 256;

typedef long int __int64_t;
typedef __int64_t int64_t;
typedef unsigned char __uint8_t;
typedef __uint8_t uint8_t;
static void gaussSeidelSolver_moved(int64_t cols, int nrb, int ncb, __mcxx_ptr_t<double> M, __mcxx_ptr_t<char> reps, unsigned char __ompif_rank, unsigned char __ompif_size, hls::stream<mcxx_outaxis>& mcxx_outPort)
{
#pragma HLS inline
	const int anrb = nrb-2;
	const int ancb = ncb-2;
	const unsigned int blocks_per_rank = anrb/(unsigned int)__ompif_size;
	diag_loop:
	//for (int R = 1; R < nrb-1; ++R)
	for (int d = 0; d < anrb+ancb-1; ++d)
	{
		RC_loop:
		//for (int C = 1; C < ncb-1; ++C)
		for (int R = 1 + (d < anrb ? d : anrb-1),
				 C = 1 + (d < anrb ? 0 : d-anrb+1);
				     (R >= 1 + 0) && (C < 1 + ancb);
				 --R, ++C)
		{
//#pragma HLS pipeline II=26
			{
				unsigned long long int __mcxx_args[2L];
				unsigned long long int __mcxx_deps[7L];
				const int localR = R - __ompif_rank*blocks_per_rank - 1;
				__mcxx_cast<int> cast_param_0;
				cast_param_0.typed = cols;
				__mcxx_args[0] = cast_param_0.raw;
				__mcxx_cast<unsigned long long int> cast_param_1;
				cast_param_1.typed = (unsigned long long int)(M + (UPADDING-1 + localR*BS)*cols + (C-1)*BS);
				__mcxx_args[1] = cast_param_1.raw;
				__mcxx_ptr_t<char> __mcxx_dep_0;
				__mcxx_dep_0 = M + (UPADDING + localR*BS - 1)*cols + LPADDING + (C-1)*BS;
				__mcxx_deps[0] = 1LLU << 58 | __mcxx_dep_0.val;
				__mcxx_ptr_t<char> __mcxx_dep_1;
				__mcxx_dep_1 = M + (UPADDING + (localR+1)*BS)*cols + LPADDING + (C-1)*BS;
				__mcxx_deps[1] = 3LLU << 58 | __mcxx_dep_1.val;
				__mcxx_ptr_t<char> __mcxx_dep_2 = M + (UPADDING + localR*BS)*cols + LPADDING + (C-1)*BS;
				__mcxx_deps[2] = 1LLU << 58 | __mcxx_dep_2.val;
				__mcxx_ptr_t<char> __mcxx_dep_3 = M + (UPADDING + (localR+1)*BS - 1)*cols + LPADDING + (C-1)*BS;
				__mcxx_deps[3] = 3LLU << 58 | __mcxx_dep_3.val;
				__mcxx_ptr_t<char> __mcxx_dep_4;
				__mcxx_dep_4 = reps + R*ncb + C-1;
				__mcxx_deps[4] = 1LLU << 58 | __mcxx_dep_4.val;
				__mcxx_ptr_t<char> __mcxx_dep_5;
				__mcxx_dep_5 = reps + R*ncb + C+1;
				__mcxx_deps[5] = 1LLU << 58 | __mcxx_dep_5.val;
				__mcxx_ptr_t<char> __mcxx_dep_6;
				__mcxx_dep_6 = reps + R*ncb + C;
				__mcxx_deps[6] = 3LLU << 58 | __mcxx_dep_6.val;
				__data_owner_info_t data_owners[2];
				const __data_owner_info_t data_owner_0 = {.size = BS*sizeof(double), .owner = calc_owner(R-1, nrb, __ompif_size, blocks_per_rank)};
				data_owners[0] = data_owner_0;
				const __data_owner_info_t data_owner_1 = {.size = BS*sizeof(double), .owner = calc_owner(R+1, nrb, __ompif_size, blocks_per_rank)};
				data_owners[1] = data_owner_1;
				mcxx_task_create(5050920982LLU, 255, 2, __mcxx_args, 7, __mcxx_deps, 0, 0, 2, data_owners, mcxx_outPort, __ompif_rank, __ompif_size, (R-1)/blocks_per_rank);
			}
			;
		}
	}
}
void solve_moved(__mcxx_ptr_t<double> matrix, __mcxx_ptr_t<char> reps, const int rbs, const int cbs, int rows, int cols, int timesteps, unsigned char __ompif_rank, unsigned char __ompif_size, hls::stream<ap_uint<8> >& mcxx_spawnInPort, hls::stream<mcxx_outaxis>& mcxx_outPort)
{
#pragma HLS inline
	const int nrb = ((rows - VPADDING) / BS) + 2;
	const int ncb = ((cols - HPADDING) / BS) + 2;
	for (int t = 0; t < timesteps;  ++t)
	{
		gaussSeidelSolver_moved(cols, nrb, ncb, matrix, reps, __ompif_rank, __ompif_size, mcxx_outPort);
	}
	mcxx_taskwait(mcxx_spawnInPort, mcxx_outPort);
}

void mcxx_write_out_port(const ap_uint<64> data, const ap_uint<3> dest, const ap_uint<1> last, hls::stream<mcxx_outaxis>& mcxx_outPort) {
#pragma HLS inline
	mcxx_outaxis axis_word;
	axis_word.data = data;
	axis_word.dest = dest;
	axis_word.last = last;
	mcxx_outPort.write(axis_word);
}

void solve_wrapper(hls::stream<ap_uint<64> >& mcxx_inPort, hls::stream<mcxx_outaxis>& mcxx_outPort, hls::stream<ap_uint<8> >& mcxx_spawnInPort, unsigned char ompif_rank, unsigned char ompif_size) {
#pragma HLS interface ap_ctrl_none port=return
#pragma HLS interface axis port=mcxx_inPort
#pragma HLS interface axis port=mcxx_outPort
#pragma HLS interface axis port=mcxx_spawnInPort
#pragma HLS stable variable=ompif_rank
#pragma HLS stable variable=ompif_size
	mcxx_inPort.read(); //command word
	__mcxx_taskId = mcxx_inPort.read();
	ap_uint<64> __mcxx_parent_taskId = mcxx_inPort.read();
	__mcxx_ptr_t<double> matrix;
	__mcxx_ptr_t<char> reps;
	int rbs;
	int cbs;
	int rows;
	int cols;
	int timesteps;
	{
#pragma HLS protocol fixed
		{
			ap_uint<8> mcxx_flags_0;
			ap_uint<64> mcxx_offset_0;
			mcxx_flags_0 = mcxx_inPort.read()(7,0);
			ap_wait();
			mcxx_offset_0 = mcxx_inPort.read();
			matrix.val = mcxx_offset_0;
		}
		ap_wait();
		{
			ap_uint<8> mcxx_flags_1;
			ap_uint<64> mcxx_offset_1;
			mcxx_flags_1 = mcxx_inPort.read()(7,0);
			ap_wait();
			mcxx_offset_1 = mcxx_inPort.read();
			reps.val = mcxx_offset_1;
		}
		ap_wait();
		{
			ap_uint<8> mcxx_flags_3;
			ap_uint<64> mcxx_offset_3;
			mcxx_flags_3 = mcxx_inPort.read()(7,0);
			ap_wait();
			__mcxx_cast<int> mcxx_arg_3;
			mcxx_arg_3.raw = mcxx_inPort.read();
			rbs = mcxx_arg_3.typed;
		}
		ap_wait();
		{
			ap_uint<8> mcxx_flags_4;
			ap_uint<64> mcxx_offset_4;
			mcxx_flags_4 = mcxx_inPort.read()(7,0);
			ap_wait();
			__mcxx_cast<int> mcxx_arg_4;
			mcxx_arg_4.raw = mcxx_inPort.read();
			cbs = mcxx_arg_4.typed;
		}
		ap_wait();
		{
			ap_uint<8> mcxx_flags_5;
			ap_uint<64> mcxx_offset_5;
			mcxx_flags_5 = mcxx_inPort.read()(7,0);
			ap_wait();
			__mcxx_cast<int> mcxx_arg_5;
			mcxx_arg_5.raw = mcxx_inPort.read();
			rows = mcxx_arg_5.typed;
		}
		ap_wait();
		{
			ap_uint<8> mcxx_flags_6;
			ap_uint<64> mcxx_offset_6;
			mcxx_flags_6 = mcxx_inPort.read()(7,0);
			ap_wait();
			__mcxx_cast<int> mcxx_arg_6;
			mcxx_arg_6.raw = mcxx_inPort.read();
			cols = mcxx_arg_6.typed;
		}
		ap_wait();
		{
			ap_uint<8> mcxx_flags_7;
			ap_uint<64> mcxx_offset_7;
			mcxx_flags_7 = mcxx_inPort.read()(7,0);
			ap_wait();
			__mcxx_cast<int> mcxx_arg_7;
			mcxx_arg_7.raw = mcxx_inPort.read();
			timesteps = mcxx_arg_7.typed;
		}
		ap_wait();
	}
	solve_moved(matrix, reps, rbs, cbs, rows, cols, timesteps, ompif_rank, ompif_size, mcxx_spawnInPort, mcxx_outPort);
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

void mcxx_task_create(const ap_uint<64> type, const ap_uint<8> instanceNum, const ap_uint<8> numArgs, const unsigned long long int args[], const ap_uint<8> numDeps, const unsigned long long int deps[], const ap_uint<8> numCopies, const __fpga_copyinfo_t copies[], hls::stream<mcxx_outaxis>& mcxx_outPort) {
#pragma HLS inline
	const ap_uint<2> destId = 2;
	ap_uint<64> tmp;
	tmp(15,8)  = numArgs;
	tmp(23,16) = numDeps;
	tmp(31,24) = numCopies;
	mcxx_write_out_port(tmp, destId, 0, mcxx_outPort);
	mcxx_write_out_port(__mcxx_taskId, destId, 0, mcxx_outPort);
	tmp(47,40) = instanceNum;
	tmp(33,0)  = type(33,0);
	mcxx_write_out_port(tmp, destId, 0, mcxx_outPort);
	for (ap_uint<4> i = 0; i < numDeps(3,0); ++i) {
#pragma HLS unroll
		mcxx_write_out_port(deps[i], destId, numArgs == 0 && numCopies == 0 && i == numDeps-1, mcxx_outPort);
	}
	for (ap_uint<4> i = 0; i < numCopies(3,0); ++i) {
#pragma HLS unroll
		mcxx_write_out_port(copies[i].copy_address, destId, 0, mcxx_outPort);
		tmp(7,0) = copies[i].flags;
		tmp(15,8) = copies[i].arg_idx;
		tmp(63,32) = copies[i].size;
		mcxx_write_out_port(tmp, destId, numArgs == 0 && i == numCopies-1, mcxx_outPort);
	}
	for (ap_uint<4> i = 0; i < numArgs(3,0); ++i) {
#pragma HLS unroll
		mcxx_write_out_port(args[i], destId, i == numArgs-1, mcxx_outPort);
	}
}

void mcxx_task_create(const ap_uint<64> type, const ap_uint<8> instanceNum, const ap_uint<8> numArgs, const unsigned long long int args[], const ap_uint<8> numDeps, const unsigned long long int deps[], const ap_uint<8> numCopies, const __fpga_copyinfo_t copies[], int numDataOwners, __data_owner_info_t data_owners[], hls::stream<mcxx_outaxis>& mcxx_outPort, unsigned char ompif_rank, unsigned char ompif_size, unsigned char owner) {
#pragma HLS inline
	const ap_uint<2> destId = 2;
	const unsigned char rank = ompif_rank;

	auto_sendrecv:
	for (ap_uint<4> i = 0; i < numDataOwners; ++i) {
#pragma HLS unroll
		const unsigned char data_owner = data_owners[i].owner;
		const unsigned int size = data_owners[i].size;
		const bool is_in = (deps[i] >> 58) & 0x1;
		if (owner != rank && data_owner == rank && is_in) {
			const unsigned long long addr = deps[i] & 0x00FFFFFFFFFFFFFF;
			const unsigned long long int dep[2] = {addr | (1LLU << 58), 0x0000100000000000LLU | (3LLU << 58)};
			OMPIF_Send((void*)addr, size, owner, 2, dep, mcxx_outPort);
		}
		else if (owner == rank && data_owner != rank && is_in) {
			const unsigned long long addr = deps[i] & 0x00FFFFFFFFFFFFFF;
			const unsigned long long int dep[2] = {addr | (2LLU << 58), 0x0000200000000000LLU | (3LLU << 58)};
			OMPIF_Recv((void*)addr, size, data_owner, 2, dep, mcxx_outPort);
		}
	}

	if (owner == rank) {
		mcxx_task_create(type, instanceNum, numArgs, args, numDeps, deps, numCopies, copies, mcxx_outPort);
	}

	/*for (ap_uint<4> i = 0; i < numDataOwners; ++i) {
#pragma HLS unroll
		const unsigned char data_owner = data_owners[i].owner;
		const unsigned int size = data_owners[i].size;
		const bool is_out = (deps[i] >> 59) & 0x1;
		if (owner == rank && data_owner == 255 && is_out) { //broadcast
			const unsigned long long addr = deps[i] & 0x00FFFFFFFFFFFFFF;
			const unsigned long long int dep[2] = {addr | (1LLU << 58), 0x0000100000000000LLU | (3LLU << 58)};
			unsigned char j = rank + 1 == ompif_size ? 0 : rank+1;
			OMPIF_Bcast((void*)addr, size, 2, dep, mcxx_outPort);
		}
		else if (owner != rank && data_owner == 255 && is_out) {
			const unsigned long long addr = deps[i] & 0x00FFFFFFFFFFFFFF;
			const unsigned long long int dep[2] = {addr | (2LLU << 58), 0x0000200000000000LLU | (3LLU << 58)};
			OMPIF_Recv((void*)addr, size, owner, 2, dep, mcxx_outPort);
		}
	}*/
}

void mcxx_taskwait(hls::stream<ap_uint<8> >& mcxx_spawnInPort, hls::stream<mcxx_outaxis>& mcxx_outPort) {
#pragma HLS inline
	ap_wait();
	mcxx_write_out_port(__mcxx_taskId, 3, 1, mcxx_outPort);
	ap_wait();
	mcxx_spawnInPort.read();
	ap_wait();
}

void OMPIF_Send(const void *data, unsigned int size, int destination, const ap_uint<8> numDeps, const unsigned long long int deps[], hls::stream<mcxx_outaxis>& mcxx_outPort) {
#pragma HLS inline
	ap_uint<64> command;
	command(7,0) = 0; //SEND
	command(15,8) = 0; //tag
	command(23,16) = destination;
	command(63, 24) = (unsigned long long int)data;
	unsigned long long int args[2] = {command, (unsigned long long int)size};
	mcxx_task_create(4294967299LU, 0xFF, 2, args, numDeps, deps, 0, 0, mcxx_outPort);
}
void OMPIF_Recv(void *data, unsigned int size, int source, const ap_uint<8> numDeps, const unsigned long long int deps[], hls::stream<mcxx_outaxis>& mcxx_outPort) {
#pragma HLS inline
	ap_uint<64> command;
	command(7,0) = 0; //RECV
	command(15,8) = 0;
	command(23,16) = source;
	command(63, 24) = (unsigned long long int)data;
	unsigned long long int args[2] = {command, (unsigned long long int)size};
	mcxx_task_create(4294967300LU, 0xFF, 2, args, numDeps, deps, 0, 0, mcxx_outPort);
}
