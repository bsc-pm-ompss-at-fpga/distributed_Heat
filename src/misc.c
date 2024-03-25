#define _POSIX_C_SOURCE 200809L

#include <assert.h>
#include <float.h>
#include <getopt.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "heat.h"

int compareOutputs(const char* refFile, int rows, int cols, int padLeft, int padRight, int padUp, double* matrix) {
	FILE* f = fopen(refFile, "r");
	if (f == NULL) {
		fprintf(stderr, "Could not open file %s\n", refFile);
		return 1;
	}

	double* refMatrix = (double*)malloc(rows*cols*sizeof(double));

	size_t n = fread(refMatrix, sizeof(double), rows*cols, f);
	fclose(f);

	if (n != rows*cols) {
		fprintf(stderr, "Could not read reference matrix\n");
		return 1;
	}

	int totalCols = padLeft+cols+padRight;

	for (int i = padUp, iRef = 0; i < padUp+rows; ++i, ++iRef) {
		for (int j = padLeft, jRef = 0; j < padLeft+cols; ++j, ++jRef) {
			int idx = i*totalCols + j;
			int idxRef = iRef*cols + jRef;
			double val = matrix[idx];
			double ref = refMatrix[idxRef];
			double err = ref == 0 ? val : fabs(val-ref)/ref;
			if (err > 1e-10) {
				double diff = fabs(val-ref); 
				fprintf(stderr, "Validation error, pos (%d, %d), expected %g found %g, diff %g err %g\n", iRef, jRef, ref, val, diff, err);
				return 1;
			}
		}
	}

	free(refMatrix);

	printf("Validation success\n");
	return 0;
}

int writeOutput(const char* refFile, int rows, int cols, int padLeft, int padRight, int padUp, double* matrix) {
	FILE* f = fopen(refFile, "w");
	if (f == NULL) {
		fprintf(stderr, "Could not open file %s\n", refFile);
		return 1;
	}

	int totalCols = padLeft+cols+padRight;

	for (int i = padUp; i < padUp+rows; ++i) {
		size_t n = fwrite(matrix + i*totalCols + padLeft, sizeof(double), cols, f);
		if (n != cols) {
			fprintf(stderr, "Could not write reference matrix\n");
			return 1;
		}
	}

	fclose(f);

	return 0;
}

int initialize(HeatConfiguration *conf, int64_t rows, int64_t cols)
{
	int pagesize = sysconf(_SC_PAGESIZE);
	assert(pagesize > 0);

	int err = posix_memalign((void **)&conf->matrix, pagesize, rows*cols*sizeof(double));
	if (err || conf->matrix == NULL) {
		fprintf(stderr, "Error: Memory cannot be allocated!\n");
		exit(1);
	}

	initializeMatrix(conf, conf->matrix, rows, cols);

	return 0;
}

int finalize(HeatConfiguration *conf)
{
	assert(conf->matrix != NULL);
	free(conf->matrix);
	conf->matrix = NULL;

	return 0;
}

int writeImage(const char *imageFileName, double *matrix, int64_t rows, int64_t cols)
{
	// RGB table
	unsigned int red[1024], green[1024], blue[1024];

	// Prepare the RGB table
	int n = 1023;
	for (int i = 0; i < 256; i++) {
		red[n] = 255; green[n] = i; blue[n] = 0;
		n--;
	}

	for (int i = 0; i < 256; i++) {
		red[n] = 255-i; green[n] = 255; blue[n] = 0;
		n--;
	}

	for (int i = 0; i < 256; i++) {
		red[n] = 0; green[n] = 255; blue[n] = i;
		n--;
	}

	for (int i = 0; i < 256; i++) {
		red[n] = 0; green[n] = 255-i; blue[n] = 255;
		n--;
	}

	// Find minimum and maximum
	double min = DBL_MAX;
	double max = -DBL_MAX;
	for (int r = UPADDING-1; r < rows-DPADDING+1; ++r) {
		for (int c = LPADDING-1; c < cols-RPADDING+1; ++c) {
			if (matrix[r*cols+c] > max)
				max = matrix[r*cols+c];
			if (matrix[r*cols+c] < min)
				min = matrix[r*cols+c];
		}
	}

	FILE *file = fopen(imageFileName, "w");
	if (file == NULL) {
		fprintf(stderr, "Error: Unable to create image file %s\n", imageFileName);
		exit(1);
	}

	fprintf(file, "P3\n");
	fprintf(file, "%ld %ld\n", cols-HPADDING+2, rows-VPADDING+2);
	fprintf(file, "255\n");

	for (int r = UPADDING-1; r < rows-DPADDING+1; ++r) {
		for (int c = LPADDING-1; c < cols-RPADDING+1; ++c) {
			int k = 0;
			if (max-min != 0) {
				k = (int)(1023.0*(matrix[r*cols+c]-min)/(max-min));
			}
			fprintf(file, "%d %d %d  ", red[k], green[k], blue[k]);
			if (c == cols-1) fprintf(file, "\n");
		}
	}

	fclose(file);

	return 0;
}

static void printUsage(int argc, char **argv)
{
	fprintf(stdout, "Usage: %s <-s size> | <-r rows -c cols> <-t timesteps> [OPTION]...\n", argv[0]);
	fprintf(stdout, "Parameters:\n");
	fprintf(stdout, "  -s, --size=SIZE          use SIZExSIZE matrix as the surface\n");
	fprintf(stdout, "  -r, --rows=ROWS          use ROWS as the number of rows of the surface\n");
	fprintf(stdout, "  -c, --cols=COLS          use COLS as the number of columns of the surface\n");
	fprintf(stdout, "  -t, --timesteps=TS       use TS as the number of timesteps\n");
	fprintf(stdout, "Optional parameters:\n");
	fprintf(stdout, "  -b, --bs=BS              use BS as the number of rows and columns of each block (default: %d)\n", DEFAULT_BS);
	fprintf(stdout, "  -R, --rbs=BS             use BS as the number of rows of each block (overrides -b option)\n");
	fprintf(stdout, "  -C, --cbs=BS             use BS as the number of columns of each block (overrides -b option)\n");
	fprintf(stdout, "  -d, --delta=DELTA        use DELTA as the residual threshold (default: %f)\n", DEFAULT_DELTA);
	fprintf(stdout, "  -f, --sources-file=NAME  get the heat sources from the NAME configuration file (default: heat.conf)\n");
	fprintf(stdout, "  -W, --no-warmup          do not perform warmup timestep (warmup enabled by default)\n");
	fprintf(stdout, "  -o, --output=NAME        save the computed matrix to the PPM file named NAME.ppm and disable warmup (disabled by default)\n");
	fprintf(stdout, "  -v, --verbose            display additional information (disabled by default)\n");
	fprintf(stdout, "  -O, --ref-output=NAME    write the resulting matrix to a reference file for future validation runs\n");
	fprintf(stdout, "  -V, --check-output=NAME  validate resulting matrix comparing with the matrix in the specified reference file\n");
	fprintf(stdout, "  -P, --parse              display only the time in seconds\n");
	fprintf(stdout, "  -h, --help               display this help and exit\n\n");
}

static void setDefaultConfiguration(HeatConfiguration *conf)
{
	conf->timesteps = 0;
	conf->convergenceTimesteps = -1;
	conf->delta = DEFAULT_DELTA;
	conf->rows = 0;
	conf->cols = 0;
	conf->rbs = DEFAULT_BS;
	conf->cbs = DEFAULT_BS;
	conf->matrix = NULL;
	conf->numHeatSources = 0;
	conf->heatSources = NULL;
	strcpy(conf->confFileName, "heat.conf");
	strcpy(conf->imageFileName, "heat.ppm");
	conf->generateImage = false;
	conf->warmup = true;
	conf->verbose = false;
	conf->createRef = false;
	conf->compareRef = false;
	conf->parse = false;
}

static int readParameters(int argc, char **argv, HeatConfiguration *conf)
{
	static struct option long_options[] = {
		{"size",         required_argument,  0, 's'},
		{"rows",         required_argument,  0, 'r'},
		{"cols",         required_argument,  0, 'c'},
		{"timesteps",    required_argument,  0, 't'},
		{"bs",           required_argument,  0, 'b'},
		{"rbs",          required_argument,  0, 'R'},
		{"cbs",          required_argument,  0, 'C'},
		{"delta",        required_argument,  0, 'd'},
		{"sources-file", required_argument,  0, 'f'},
		{"output",       required_argument,  0, 'o'},
		{"no-warmup",    no_argument,        0, 'W'},
		{"verbose",      no_argument,        0, 'v'},
		{"ref-output",   required_argument,  0, 'O'},
		{"check-output", required_argument,  0, 'V'},
		{"parse",        no_argument,        0, 'P'},
		{"help",         no_argument,        0, 'h'},
		{0, 0, 0, 0}
	};

	int c, index;
	int bs = DEFAULT_BS;
	int rbs = 0, cbs = 0;

	while ((c = getopt_long(argc, argv, "ho:f:s:r:c:t:b:R:C:d:WvO:V:P", long_options, &index)) != -1) {
		switch (c) {
			case 'h':
				printUsage(argc, argv);
				return 2;
			case 'v':
				conf->verbose = true;
				break;
			case 'f':
				if (strlen(optarg) >= MAX_STRING_SIZE) {
					fprintf(stderr, "Error: Configuration name is too long!\n");
					return 1;
				}
				strcpy(conf->confFileName, optarg);
				break;
			case 'o':
				conf->generateImage = true;
				conf->warmup = false;
				if (strlen(optarg) >= MAX_STRING_SIZE) {
					fprintf(stderr, "Error: Image name is too long!\n");
					return 1;
				}
				strcpy(conf->imageFileName, optarg);
				break;
			case 'W':
				conf->warmup = false;
				break;
			case 's':
				conf->rows = atoi(optarg);
				conf->cols = atoi(optarg);
				assert(conf->rows > 0 && conf->cols > 0);
				break;
			case 'r':
				conf->rows = atoi(optarg);
				assert(conf->rows > 0);
				break;
			case 'c':
				conf->cols = atoi(optarg);
				assert(conf->cols > 0);
				break;
			case 't':
				conf->timesteps = atoi(optarg);
				assert(conf->timesteps > 0);
				break;
			case 'b':
				bs = atoi(optarg);
				assert(bs > 0);
				break;
			case 'R':
				rbs = atoi(optarg);
				assert(rbs > 0);
				break;
			case 'C':
				cbs = atoi(optarg);
				assert(cbs > 0);
				break;
			case 'd':
				conf->delta = atof(optarg);
				assert(conf->delta > 0.0);
				break;
			case 'O':
				if (strlen(optarg) >= MAX_STRING_SIZE) {
					fprintf(stderr, "Error: Output file name is too long!\n");
					return 1;
				}
				conf->createRef = true;
				strcpy(conf->outFileName, optarg);
				break;
			case 'V':
				if (strlen(optarg) >= MAX_STRING_SIZE) {
					fprintf(stderr, "Error: Reference file name is too long!\n");
					return 1;
				}
				conf->compareRef = true;
				strcpy(conf->refFileName, optarg);
				break;
			case 'P':
				conf->parse = true;
				break;
			default:
				return 1;
		}
	}

	conf->rbs = (rbs == 0) ? bs : rbs;
	conf->cbs = (cbs == 0) ? bs : cbs;

	if (!conf->rows || !conf->cols || !conf->timesteps || !conf->rbs || !conf->cbs) {
		printUsage(argc, argv);
		return 1;
	}

	return 0;
}

int readConfiguration(int argc, char **argv, HeatConfiguration *conf)
{
	char line[MAX_STRING_SIZE];

	setDefaultConfiguration(conf);

	// Read the execution parameters
	if (readParameters(argc, argv, conf) != 0) {
		return 1;
	}

	FILE *file = fopen(conf->confFileName, "r");
	if (file == NULL) {
		fprintf(stderr, "Error: Unable to open configuration file %s!\n", conf->confFileName);
		return 1;
	}

	if (!fgets(line, MAX_STRING_SIZE, file)) {
		fprintf(stderr, "Error: Configuration file is not correct!\n");
		return 1;
	}

	int n = sscanf(line, "%d", &(conf->numHeatSources));
	if (n != 1) {
		fprintf(stderr, "Error: Configuration file not correct!\n");
		return 1;
	}
	if (conf->numHeatSources < 1) {
		fprintf(stderr, "Error: Configuration file must have at least one heat source!\n");
		return 1;
	}

	conf->heatSources = (HeatSource *) malloc(sizeof(HeatSource)*conf->numHeatSources);
	assert(conf->heatSources != NULL);

	for (int i = 0; i < conf->numHeatSources; i++) {
		if (!fgets(line, MAX_STRING_SIZE, file)) {
			fprintf(stderr, "Error: Configuration file is not correct!\n");
			return 1;
		}

		n = sscanf(line, "%f %f %f %f",
			&(conf->heatSources[i].row),
			&(conf->heatSources[i].col),
			&(conf->heatSources[i].range),
			&(conf->heatSources[i].temperature));

		if (n != 4) {
			fprintf(stderr, "Error: Configuration file not correct!\n");
			return 1;
		}
	}

	fclose(file);

	return 0;
}

void refineConfiguration(HeatConfiguration *conf, int64_t rowValue, int64_t colValue)
{
	assert(conf->rows > 0);
	assert(conf->cols > 0);
	assert(conf->timesteps > 0);

	if (conf->rows%rowValue) {
		fprintf(stderr, "Warning: The number of rows (%ld) is not divisible by %ld. Rounding it...\n", conf->rows, rowValue);
		// Make the number of rows divisible by the value
		conf->rows = ROUND(conf->rows, rowValue);
	}
	if (conf->cols%colValue) {
		fprintf(stderr, "Warning: The number of cols (%ld) is not divisible by %ld. Rounding it...\n", conf->cols, colValue);
		// Make the number of cols divisible by the value
		conf->cols = ROUND(conf->cols, colValue);
	}
}

void printConfiguration(const HeatConfiguration *conf)
{
	fprintf(stdout, "Rows x Cols       : %ld x %ld\n", conf->rows, conf->cols);
	fprintf(stdout, "Block size        : %d x %d\n", conf->rbs, conf->cbs);
	fprintf(stdout, "Timesteps         : %d\n", conf->timesteps);
	fprintf(stdout, "Delta             : %f\n", conf->delta);
	fprintf(stdout, "Num. heat sources : %d\n", conf->numHeatSources);

	for (int i = 0; i < conf->numHeatSources; i++) {
		fprintf(stdout, "  %2d: (%2.2f, %2.2f) %2.2f %2.2f\n", i+1,
			conf->heatSources[i].row,
			conf->heatSources[i].col,
			conf->heatSources[i].range,
			conf->heatSources[i].temperature
		);
	}
}

void initializeMatrix(const HeatConfiguration *conf, double *matrix, int64_t rows, int64_t cols)
{
	// Set all elements to zero
	memset(matrix, 0, rows*cols*sizeof(double));

	int onepad_rows = rows-VPADDING+2;
	int onepad_cols = cols-HPADDING+2;

	for (int i = 0; i < conf->numHeatSources; i++) {
		const HeatSource *src = &(conf->heatSources[i]);

		// Initialize top row
		for (int c = LPADDING-1; c < cols-RPADDING+1; ++c) {
			int coord_c = c-LPADDING+1;
			double dist = sqrt(pow((double)coord_c/(double)onepad_cols-src->col, 2) + pow(src->row, 2));
			if (dist <= src->range) {
				matrix[(UPADDING-1)*cols + c] += (src->range-dist)/src->range*src->temperature;
			}
		}

		// Initialize bottom row
		for (int c = LPADDING-1; c < cols-RPADDING+1; ++c) {
			int coord_c = c-LPADDING+1;
			double dist = sqrt(pow((double)coord_c/(double)onepad_cols-src->col, 2) + pow(1-src->row, 2));
			if (dist <= src->range) {
				matrix[(rows-DPADDING)*cols+c] += (src->range-dist)/src->range*src->temperature;
			}
		}

		// Initialize left column
		for (int r = UPADDING; r < rows-DPADDING; ++r) {
			int coord_r = r-UPADDING+1;
			double dist = sqrt(pow(src->col, 2) + pow((double)coord_r/(double)onepad_rows-src->row, 2));
			if (dist <= src->range) {
				matrix[r*cols + LPADDING-1] += (src->range-dist)/src->range*src->temperature;
			}
		}

		// Initialize right column
		for (int r = UPADDING; r < rows-DPADDING; ++r) {
			int coord_r = r-UPADDING+1;
			double dist = sqrt(pow(1-src->col, 2) + pow((double)coord_r/(double)onepad_rows-src->row, 2));
			if (dist <= src->range) {
				matrix[r*cols+cols-RPADDING] += (src->range-dist)/src->range*src->temperature;
			}
		}
	}
}

double getTime()
{
	struct timespec tv;
	clock_gettime(CLOCK_MONOTONIC, &tv);
	return tv.tv_sec+1e-9*tv.tv_nsec;
}
