#include <stdio.h>
#include <stdlib.h>
#include <windows.h>
#include <time.h>

// ASM kernel
extern void matvec_ymm_avx2_asm(int n, float* A, float* x, float* y);


// Helper functions
void print_array(int, float[]);
void init_data( float*, float*, int);
void matvec_c(int, float*, float*, float*);
double time_run(void (*func)(int, float*, float*, float*),
	int n, float* A, float* X, float* Y);

// Macros for n
#define TWO_TO_TEN 1 << 10
#define TWO_TO_THIRTEEN 1 << 13
#define TWO_TO_FOURTEEN 1 << 14

// Macro for total number of runs
#define RUNS 30

int main() {
	int n = TWO_TO_TEN;


	size_t bytesA = (size_t)n * (size_t)n * sizeof(float);
	size_t bytesV = (size_t)n * sizeof(float);


	float* A = (float*)malloc(bytesA);
	float* X = (float*)malloc(bytesV);
	float* Y_c = (float*)malloc(bytesV);
	float* Y_asm = (float*)malloc(bytesV);

	double exec_time_c[RUNS], exec_time_asm[RUNS];
	double exec_time_ave_c = 0, exec_time_ave_asm = 0;

	init_data( A, X,n);

	printf("Vector size n: %d\n", n);

	int curr_run;
	for (curr_run = 0; curr_run < RUNS; curr_run++) {
		printf("\n=== Run %02d of %d ===\n", curr_run + 1, RUNS);

		printf("X: ");
		print_array(16, X);
		puts("");

		// Run and get execution time of each kernel function
		exec_time_c[curr_run] = time_run(matvec_c, n, A, X, Y_c);
		exec_time_asm[curr_run] = time_run(matvec_ymm_avx2_asm, n, A, X, Y_asm);
		exec_time_ave_c += exec_time_c[curr_run];
		exec_time_ave_asm += exec_time_asm[curr_run];

		printf("kernel (C  ):");
		print_array(16, Y_c);
		printf("kernel (ymm):");
		print_array(16, Y_asm);


		printf("Execution time (C  ): %f ms\n", exec_time_c[curr_run]);
		printf("Execution time (ymm): %f ms\n", exec_time_asm[curr_run]);
	}

	exec_time_ave_c /= RUNS;
	exec_time_ave_asm /= RUNS;

	puts("All runs finished successfully with equal output.");
	printf("\nAverage execution time (C  ): %f ms\n", exec_time_ave_c);
	printf("Average execution time (ymm): %f ms\n", exec_time_ave_asm);

	return 0;
}

void matvec_c(int n, float* A, float* x, float* y) {
	for (int i = 0; i < n; ++i) {
		double s = 0.0;                   
		float* row = A + (size_t)i * n;
		for (int j = 0; j < n; ++j)
			s += (double)row[j] * (double)x[j];
		y[i] = (float)s;
	}
}

// Print the first n elements of arr
void print_array(int n, float arr[]) {
	int i;
	for (i = 0; i < n; i++)
		printf("%.1f ", arr[i]);
	printf("\n");
}

void init_data(float* A, float* x, int n) {
	for (int i = 0; i < n; ++i) {
		float* row = A + (size_t)i * n;
		for (int j = 0; j < n; ++j) {
			row[j] = 1.0f / (float)(i + j + 1);
		}
	}
	for (int j = 0; j < n; ++j) x[j] = 1.0f + 0.001f * (float)j;
}



double time_run(void (*func)(int,float*, float*, float*),
	int n, float* A, float* X, float* Y) {

	LARGE_INTEGER frequency, start, end;
	double exec_time;

	QueryPerformanceFrequency(&frequency);
	QueryPerformanceCounter(&start);

	func(n, A, X, Y);

	QueryPerformanceCounter(&end);
	exec_time = (double)(end.QuadPart - start.QuadPart)*100000 / frequency.QuadPart;

	return exec_time;
}


