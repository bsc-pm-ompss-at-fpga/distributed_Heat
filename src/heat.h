#ifndef HEAT_H
#define HEAT_H

#include <stdbool.h>
#include <stdint.h>

#define IGNORE_RESIDUAL ((double) -1.0)
#define DEFAULT_DELTA ((double) 0.00005)
#define DEFAULT_BS BLOCK_SIZE
#define MAX_STRING_SIZE 100

#define ROUND(a, b) ((((a) + (b) - 1) / (b)) * (b))
#define LPADDING (64/sizeof(double))
#define RPADDING LPADDING
#define HPADDING (LPADDING+RPADDING)
#define UPADDING 1
#define DPADDING 1
#define VPADDING (UPADDING+DPADDING)

#ifndef BLOCK_SIZE
#error BLOCK_SIZE not defined
#endif

#ifndef NUM_CBLOCKS
#error NUM_CBLOCKS not defined
#endif

static const int BS = BLOCK_SIZE;
static const int NUM_COMPUTE_BLOCKS = NUM_CBLOCKS;

typedef struct {
	float row;
	float col;
	float range;
	float temperature;
} HeatSource;

typedef struct {
	int timesteps;
	int convergenceTimesteps;
	double delta;
	int64_t rows;
	int64_t cols;
	int rbs;
	int cbs;
	double *matrix;
	int numHeatSources;
	HeatSource *heatSources;
	char confFileName[MAX_STRING_SIZE];
	char imageFileName[MAX_STRING_SIZE];
	bool generateImage;
	bool warmup;
	bool verbose;
	bool createRef;
	bool compareRef;
	char refFileName[MAX_STRING_SIZE];
	char outFileName[MAX_STRING_SIZE];
	bool parse;
} HeatConfiguration;


int compareOutputs(const char* refFile, int rows, int cols, int padLeft, int padRight, int padUp, double* matrix);
int writeOutput(const char* refFile, int rows, int cols, int padLeft, int padRight, int padUp, double* matrix);
int initialize(HeatConfiguration *conf, int64_t rows, int64_t cols);
int finalize(HeatConfiguration *conf);
int writeImage(const char *fileName, double *matrix, int64_t rows, int64_t cols);
int readConfiguration(int argc, char **argv, HeatConfiguration *conf);
void refineConfiguration(HeatConfiguration *conf, int64_t rowValue, int64_t colValue);
void printConfiguration(const HeatConfiguration *conf);
void initializeMatrix(const HeatConfiguration *conf, double *matrix, int64_t rows, int64_t cols);
double getTime();

#pragma oss task device(broadcaster) inout(matrix[0]) in(reps[0])
void solve(double *matrix, char* reps, const int rbs, const int cbs, int rows, int cols, int timesteps);
//#pragma oss task device(fpga) num_instances(NUM_COMPUTE_BLOCKS)
void computeBlock(const int cols, unsigned long long int M);
double computeBlockResidual(const int64_t rows, const int64_t cols, const int rstart, const int rend, const int cstart, const int cend, double M[rows][cols]);

#endif // HEAT_H
