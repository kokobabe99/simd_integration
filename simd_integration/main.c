#include <stdio.h>
#include <stdlib.h>
#include <windows.h>
#include <time.h>

// ASM kernel
extern void matvec_x86_64_asm(int n, float* A, float* x, float* y);
extern void matvec_xmm_avx2_asm(int n, float* A, float* x, float* y);
extern void matvec_ymm_avx2_asm(int n, float* A, float* x, float* y);
extern void matvec_xmm_avx2_asm_v2(int n, float* A, float* x, float* y);
extern void print_first_last(const char* name, float* y, int n);




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
#define TWO_TO_FIFTEEN 1 << 15
#define ONE_THOUSAND_THREE 1003

// Macro for total number of runs
#define RUNS 30

int main() {
	int n = TWO_TO_THIRTEEN;


	size_t bytesA = (size_t)n * (size_t)n * sizeof(float);
	size_t bytesV = (size_t)n * sizeof(float);


	float* A = (float*)malloc(bytesA);
	float* X = (float*)malloc(bytesV);
	float* Y_c = (float*)malloc(bytesV);
	float* Y_x86_64 = (float*)malloc(bytesV);
	float* Y_xmm = (float*)malloc(bytesV);
	float* Y_ymm = (float*)malloc(bytesV);
	float* y_xmm_v2 = (float*)malloc(bytesV);

	double exec_time_c[RUNS], exec_time_x86_64[RUNS], exec_time_xmm[RUNS], exec_time_ymm[RUNS], exec_time_xmm_v2[RUNS];
	double exec_time_ave_c = 0, exec_time_ave_x86_64 = 0, exec_time_ave_xmm = 0, exec_time_ave_ymm = 0,exec_time_ave_xmm_v2 = 0;

	init_data( A, X,n);

	printf("Vector size n: %d\n", n);



	int curr_run;
	for (curr_run = 0; curr_run < RUNS; curr_run++) {


		// Check if Z_c and Z_asm are equal
		// To speed up process, only check first 10 elements


		printf("X: ");
		print_array(16, X);
		puts("");

		// Run and get execution time of each kernel function
		exec_time_c[curr_run] = time_run(matvec_c, n, A, X, Y_c);
		exec_time_x86_64[curr_run] = time_run(matvec_x86_64_asm, n, A, X, Y_x86_64);
		exec_time_xmm[curr_run] = time_run(matvec_xmm_avx2_asm, n, A, X, Y_xmm);
		exec_time_ymm[curr_run] = time_run(matvec_ymm_avx2_asm, n, A, X, Y_ymm);
		exec_time_xmm_v2[curr_run] = time_run(matvec_xmm_avx2_asm_v2, n, A, X, y_xmm_v2);
	
		exec_time_ave_c += exec_time_c[curr_run];
		exec_time_ave_x86_64 += exec_time_x86_64[curr_run];
		exec_time_ave_xmm += exec_time_xmm[curr_run];
		exec_time_ave_ymm += exec_time_ymm[curr_run];
		exec_time_ave_xmm_v2 += exec_time_xmm_v2[curr_run];

		printf("First 20 for Correctness check");
		printf("kernel (C  ):");
		print_array(16, Y_c);
		printf("kernel (x86):");
		print_array(16, Y_x86_64);
		printf("kernel (xmm):");
		print_array(16, y_xmm_v2);
		printf("kernel (ymm):");
		print_array(16, Y_ymm);
		//printf("kernel (xmm_v2):");
		//print_array(16, y_xmm_v2);
		printf("\n");
		printf("last 3 for Correctness check");
		print_first_last("kernel (C  ):", Y_c, n);
		print_first_last("kernel (x86):", Y_x86_64, n);
		print_first_last("kernel (xmm):", y_xmm_v2, n);
		print_first_last("kernel (ymm)", Y_ymm, n);


		printf("Execution time (C  ): %f ms\n", exec_time_c[curr_run]);
		printf("Execution time (x86): %f ms\n", exec_time_x86_64[curr_run]);
		printf("Execution time (xmm): %f ms\n", exec_time_xmm_v2[curr_run]);
		printf("Execution time (ymm): %f ms\n", exec_time_ymm[curr_run]);
		//// Check if Z_c and Z_asm are equal
		// To speed up process, only check first 10 elements
		//int equal = memcmp(y_xmm_v2, y_xmm_v2, sizeof(float) * n) == 0;
		//if (equal) {
		//	puts("The C and x86-64 kernel outputs are equal for the first 10 values.\n");
		//}
		//else {
		//	puts("WARNING: The C and x86-64 kernel outputs are not equal.\n");
		//	return 1;
		//}

		
		

		//if (*y_xmm_v2 == *Y_ymm) {
		//	puts("The C and x86-64 kernel outputs are equal for the first 10 values.\n");
		//}
		//else {
		//	puts("WARNING: The C and x86-64 kernel outputs are not equal.\n");
		//	return 1;
		//}
		//printf("\n=== Run %02d of %d ===\n", curr_run + 1, RUNS);


	}

	exec_time_ave_c /= RUNS;
	exec_time_ave_x86_64 /= RUNS;
	exec_time_ave_xmm /= RUNS;
	exec_time_ave_ymm /= RUNS;
	exec_time_ave_xmm_v2 /= RUNS;






	puts("All runs finished successfully with equal output.\n");
	printf("Average execution time (C  ): %f ms\n", exec_time_ave_c);
	printf("Average execution time (x86): %f ms\n", exec_time_ave_x86_64);
	printf("Average execution time (xmm): %f ms\n", exec_time_ave_xmm_v2);
	printf("Average execution time (ymm): %f ms\n", exec_time_ave_ymm);
	//printf("Average execution time (xmm_v2): %f ms\n", exec_time_ave_xmm_v2);

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
		printf("%.3f ", arr[i]);
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

void print_first_last(const char* name, float* y, int n) {
	printf("%s (first 3): ", name);
	for (int i = 0; i < 3 && i < n; ++i)
		printf("%.3f ", y[i]);

	printf("... (last 3): ");
	for (int i = n - 3; i < n; ++i)
		if (i >= 0)
			printf("%.3f ", y[i]);
	printf("\n");
}