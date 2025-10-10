#include <stdio.h>
#include <stdlib.h>
#include <windows.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>


#define EPSILON 1e-5f 

// ASM kernel
extern void matvec_x86_64_asm(int n, float* A, float* x, float* y);
extern void matvec_xmm_avx2_asm(int n, float* A, float* x, float* y);
extern void matvec_ymm_avx2_asm(int n, float* A, float* x, float* y);
extern void matvec_xmm_avx2_asm_v2(int n, float* A, float* x, float* y);
extern void print_first_last(const char* name, float* y, int n);
extern boolean compare_results(const char* name_a, const float* A,
	const char* name_b, const float* B, int n);

// Helper functions
void print_array(int, float[]);
void init_data( float*, float*, int);
void matvec_c(int, float*, float*, float*);
double time_run(void (*func)(int, float*, float*, float*),
	int n, float* A, float* X, float* Y);
double compute_geometric_mean(double reference[], double exec_times[]);

// Macros for n
#define TWO_TO_TEN 1 << 10
#define TWO_TO_THIRTEEN 1 << 13
#define TWO_TO_FOURTEEN 1 << 14
#define ONE_THOUSAND_THREE 1003

// Macro for total number of runs
#define RUNS 30

int main() {
	int n = ONE_THOUSAND_THREE;


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

		//printf("First 20 for Correctness check\n");
		//printf("kernel (C  ):");
		//print_array(16, Y_c);
		//printf("kernel (x86):");
		//print_array(16, Y_x86_64);
		//printf("kernel (xmm):");
		//print_array(16, y_xmm_v2);
		//printf("kernel (ymm):");
		//print_array(16, Y_ymm);
		//printf("\n");
		//printf("last 3 for Correctness check\n");
		//print_first_last("kernel (C  ):", Y_c, n);
		//print_first_last("kernel (x86):", Y_x86_64, n);
		//print_first_last("kernel (xmm):", y_xmm_v2, n);
		//print_first_last("kernel (ymm):", Y_ymm, n);


		printf("Execution time (C  ): %f ms\n", exec_time_c[curr_run]);
		printf("Execution time (x86): %f ms\n", exec_time_x86_64[curr_run]);
		printf("Execution time (xmm): %f ms\n", exec_time_xmm_v2[curr_run]);
		printf("Execution time (ymm): %f ms\n", exec_time_ymm[curr_run]);

		boolean res = compare_results("C", Y_c, "X86",Y_x86_64, n);
		if (res) {
			puts("The C and X86 kernel outputs are equal\n");
		}
		else {
			puts("WARNING: The C and X86 kernel outputs are not equal.\n");
			return 1;
		}

		res = compare_results("C", Y_c, "xmm", y_xmm_v2, n);
		if (res) {
			puts("The C and XMM kernel outputs are equal\n");
		}
		else {
			puts("WARNING: The C and XMM kernel outputs are not equal.\n");
			return 1;
		}

		res = compare_results("C", Y_c, "ymm", Y_ymm, n);
		if (res) {
			puts("The C and ymm kernel outputs are equal\n");
		}
		else {
			puts("WARNING: The C and ymm kernel outputs are not equal.\n");
			return 1;
		}
	
	}

	exec_time_ave_c /= RUNS;
	exec_time_ave_x86_64 /= RUNS;
	exec_time_ave_xmm /= RUNS;
	exec_time_ave_ymm /= RUNS;
	exec_time_ave_xmm_v2 /= RUNS;

	double exec_time_geom_c = compute_geometric_mean(exec_time_c, exec_time_c);
	double exec_time_geom_x86_64 = compute_geometric_mean(exec_time_c,exec_time_x86_64);
	double exec_time_geom_xmm = compute_geometric_mean(exec_time_c, exec_time_xmm_v2);
	double exec_time_geom_ymm = compute_geometric_mean(exec_time_c, exec_time_ymm);



	puts("All runs finished successfully with equal output.\n");
	printf("Average len %d execution time (C  ): %f ms\n",n,exec_time_ave_c);
	printf("Average len %d execution time (x86): %f ms\n",n,exec_time_ave_x86_64);
	printf("Average len %d execution time (xmm): %f ms\n",n,exec_time_ave_xmm_v2);
	printf("Average len %d execution time (ymm): %f ms\n",n,exec_time_ave_ymm);

	puts("\nUsing C kernel as reference for geometric mean:");
	printf("Geometric mean len %d execution time (C  ): %f, %fx faster\n", n, exec_time_geom_c, 1/exec_time_geom_c);
	printf("Geometric mean len %d execution time (x86): %f, %fx faster\n", n, exec_time_geom_x86_64, 1/exec_time_geom_x86_64);
	printf("Geometric mean len %d execution time (xmm): %f, %fx faster\n", n, exec_time_geom_xmm, 1/exec_time_geom_xmm);
	printf("Geometric mean len %d execution time (ymm): %f, %fx faster\n", n, exec_time_geom_ymm, 1/exec_time_geom_ymm);
	//printf("Average execution time (xmm_v2): %f ms\n", exec_time_ave_xmm_v2);

	return 0;
}

void matvec_c(int n, const float* A, const float* x, float* y) {
	for (int i = 0; i < n; ++i) {
		const float* row = A + (size_t)i * n;
		double s = 0.0;  // accumulate in double for numerical stability
		for (int j = 0; j < n; ++j) {
			s += (double)row[j] * (double)x[j];
		}
		y[i] = (float)s;
	}
}

void init_data(float* A, float* x, int n) {
	// Matrix A: A[i*n + j] = 1 / ((i+1) + (j+1) - 1) = 1 / (i + j + 1)
	for (int i = 0; i < n; ++i) {
		float* row = A + (size_t)i * n;
		for (int j = 0; j < n; ++j) {
			row[j] = (float)(1.0 / (double)(i + j + 1));
		}
	}

	// Vector x: x[j] = sin(j*0.01) * cos(j*0.007) + 1.0
	for (int j = 0; j < n; ++j) {
		x[j] = (float)(sin((double)j * 0.01) * cos((double)j * 0.007) + 1.0);
	}
}


// Print the first n elements of arr
void print_array(int n, float arr[]) {
	int i;
	for (i = 0; i < n; i++)
		printf("%.4f ", arr[i]);
	printf("\n");
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
		printf("%.4f ", y[i]);

	printf("... (last 3): ");
	for (int i = n - 3; i < n; ++i)
		if (i >= 0)
			printf("%.4f ", y[i]);
	printf("\n");
}


static inline boolean nearly_equal(float a, float b, float eps) {
	float diff = fabsf(a - b);
	float scale = fmaxf(1.0f, fmaxf(fabsf(a), fabsf(b)));
	return diff <= eps * scale;
}

boolean compare_results(const char* name_a,  const float* A,
	const char* name_b, const float* B, int n) {
	int tail_start = n > 3 ? n - 3 : 0;

	printf("kernel (%s) first3: ", name_a);
	for (int i = 0; i < 3 && i < n; ++i) printf("%.4f ", A[i]);
	printf(" | last3: ");
	for (int i = tail_start; i < n; ++i) printf("%.4f ", A[i]);
	printf("\n");

	printf("kernel (%s) first3: ", name_b);
	for (int i = 0; i < 3 && i < n; ++i) printf("%.4f ", B[i]);
	printf(" | last3: ");
	for (int i = tail_start; i < n; ++i) printf("%.4f ", B[i]);
	printf("\n");

	for (int i = 0; i < n; ++i) {
		if (!nearly_equal(A[i], B[i], EPSILON)) {
			printf("Mismatch at %d: %s=%.8f, %s=%.8f (eps=%g)\n",
				i, name_a, A[i], name_b, B[i], EPSILON);
			return false;	
		}
	}
	printf("%s vs %s: numerically equal within epsilon=%g\n",
		name_a, name_b, EPSILON);
	return true;
}

double compute_geometric_mean(double reference[], double exec_times[]) {
	double product = 1;
	for (int i = 0; i < RUNS; i++) {
		product *= exec_times[i] / reference[i];
	}
	return pow(product, 1.0 / RUNS);
}