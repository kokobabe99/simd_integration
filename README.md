# simd_integration

## Project Overview

---

contains the `kernels` in

- `C program`
- `x86-64 assembly language`
- x86 `SIMD` AVX2 assembly language using `XMM register`
- x86 `SIMD` AVX2 assembly language using `YMM register`

Each kernel is to perform a `matrix vector product`.matrix-vector product `y = A·x (single-precision)` in RUN times.

### Project Layout

```
     ├── docs.                     (Screenshots of output)
     ├── LICENSE.txt
     ├── README.md
     ├── simd_integration
     │   ├── asmfunc1.asm          (x86)
     │   ├── asmfunc3.asm          (SIMD XMM)
     │   ├── asmfunc4.asm          (SIMD YMM)
     │   ├── main.c (main program)
     │   ├── simd_integration.vcxproj
     │   └── simd_integration.vcxproj.filters
     └── simd_integration.sln
```

### Enviroment

```
     Tools:Visual Studio,NASM
     Platform:aws 32G 64bit window 2022
```

## Main Module

---

### What is matrix vector product ?

A **matrix–vector product** (also known as **matrix–vector multiplication**) is a fundamental operation in linear algebra where a matrix \( A \) multiplies a vector \( x \) to produce another vector \( y \).

### Example

![martix_product](/docs/matrix_product.png)

### (C)

```C
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
```

### Results

```
Matrix sizes: n ∈ {1003,2^10,2^13,2^14} RUN in 30 times

#define RUNS 10
```

| n            |      C (ms) |   x86 (ms) |   XMM (ms) |  YMM (ms) | x86 Speedup VS C | XMM Speedup VS C | YMM Speedup VS C |
| ------------ | ----------: | ---------: | ---------: | --------: | ---------------: | ---------------: | ---------------: |
| 2^10 (1024)  |     430.762 |    156.672 |     30.243 |    30.947 |            2.75× |           14.24× |           13.92× |
| 2^13 (8192)  |  25,124.121 |  7,881.562 |  2,400.416 | 2,114.746 |            3.19× |           10.47× |           11.88× |
| 2^14 (16384) | 100,416.107 | 31,844.125 | 10,230.686 | 8,981.545 |            3.15× |            9.82× |           11.18× |
| 1003         |     450.037 |    133.948 |     36.630 |    26.653 |            3.36× |           12.29× |           16.88× |

- 2^10

```
//for n
#define TWO_TO_TEN 1 << 10
```

![210](/docs/210.png)

- 2^13

```
//for n
#define TWO_TO_THIRTEEN 1 << 13
```

![213](/docs/213.png)

- 2^14

```
//for n
#define TWO_TO_FOURTEEN 1 << 14
```

![214](/docs/214.png)

- how about 1003 (to testify Boundary-check) ?

```
//for n
#define ONE_THOUSAND_THREE 1003
```

![1003](/docs/1003.png)

## Implementations

---

### initial params

```C
	int n = TWO_TO_THIRTEEN;

	size_t bytesA = (size_t)n * (size_t)n * sizeof(float);
	size_t bytesV = (size_t)n * sizeof(float);

	float* A = (float*)malloc(bytesA);
	float* X = (float*)malloc(bytesV);
	float* Y_c = (float*)malloc(bytesV);
	float* Y_x86_64 = (float*)malloc(bytesV);
	float* Y_xmm = (float*)malloc(bytesV);
	float* Y_ymm = (float*)malloc(bytesV);
```

#### Timing Method

```C
QueryPerformanceFrequency
```

### Correct call parameters

```C

// Lambda function call
double time_run(void (*func)(int,float*, float*, float*),int n, float* A, float* X, float* Y) {

	LARGE_INTEGER frequency, start, end;
	double exec_time;

	QueryPerformanceFrequency(&frequency);
	QueryPerformanceCounter(&start);

	func(n, A, X, Y);

	QueryPerformanceCounter(&end);
	exec_time = (double)(end.QuadPart - start.QuadPart)*100000 / frequency.QuadPart;

	return exec_time;
}
```

```C

     // Run and get execution time of each kernel function
     exec_time_c[curr_run] = time_run(matvec_c, n, A, X, Y_c);
     exec_time_x86_64[curr_run] = time_run(matvec_x86_64_asm, n, A, X, Y_x86_64);
	exec_time_xmm[curr_run] = time_run(matvec_xmm_avx2_asm, n, A, X, Y_xmm);
	exec_time_ymm[curr_run] = time_run(matvec_ymm_avx2_asm, n, A, X, Y_ymm);
```

### Assembly XMM & YMM

```assembly languary
; y = A * x (single-precision float)
; x86-64 kernel
; RCX = n (int)
; RDX = A (float*)
; R8 = x (float*)
; R9 = y (float*)
```

### Volatile handle

```
    push r12
    push rdi
    push rsi

done:
    pop rsi
    pop rdi
    pop r12
    ret

```

### Boundary check handle

C reference uses double accumulation for numerical stability.

XMM kernel: 128-bit, processes 4 floats/iter; n4 = n & ~3, scalar tail for n%4.

YMM kernel: 256-bit, processes 8 floats/iter; n8 = n & ~7, scalar tail for n%8.

use `and` to check boundary and save the reminder into another funcion
Example if `ymm` to ` and     r12d, -8`,`xmm` to `and     r12d, -4`

Example: n=1003 (YMM) with only last n%W elements non-zero (W=8) to validate reminder path.

```C
    mov     r12d, ecx
    and     r12d, -8
    xor     edi, edi
loop_inner:
    cmp     edi, r12d
    jge     boundary_check
    vmovups ymm1, [r11 + rdi*4]
    vmovups ymm2, [r8  + rdi*4]
    vmulps  ymm1, ymm1, ymm2
    vaddps  ymm0, ymm0, ymm1
    add     edi, 8
    jmp     loop_inner

 boundary_check:
    vextractf128 xmm1, ymm0, 1
    vhaddps xmm0, xmm0, xmm1
    vhaddps xmm0, xmm0, xmm0
    vhaddps xmm0, xmm0, xmm0
    vmovaps xmm3, xmm0
    mov eax, r12d

 reminder_add:
    cmp eax,ecx
    jge store
    vmovss   xmm1, dword [r11 + rax*4]
    vmulss   xmm1, xmm1, dword [r8 + rax*4]
    vaddss   xmm3, xmm3, xmm1
    inc eax
    jmp reminder_add
```

### Correctness check

```C


boolean compare_results(const char* name_a,  const float* A,
	const char* name_b, const float* B, int n) {
	int tail_start = n > 3 ? n - 3 : 0;

	printf("kernel (%s) first3: ", name_a);
	for (int i = 0; i < 3 && i < n; ++i) printf("%.6f ", A[i]);
	printf(" | last3: ");
	for (int i = tail_start; i < n; ++i) printf("%.6f ", A[i]);
	printf("\n");

	printf("kernel (%s) first3: ", name_b);
	for (int i = 0; i < 3 && i < n; ++i) printf("%.6f ", B[i]);
	printf(" | last3: ");
	for (int i = tail_start; i < n; ++i) printf("%.6f ", B[i]);
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

```

## Comparative Results

---

## Analysis

The SIMD YMM kernel consistently outperformed both the scalar (x86) and SIMD XMM implementations — achieving up to ~28× speedup compared to the baseline C version for small-to-medium matrices.

The SIMD kernels achieve significant acceleration because each instruction processes multiple data elements in parallel.
While the scalar C implementation handles one float at a time, the XMM version processes 4 floats per iteration (128-bit), and the YMM version processes 8 floats per iteration (256-bit).
This reduces loop iterations, branch overhead, and instruction decoding, allowing the CPU’s wide vector pipelines and execution units to operate at full efficiency.
In addition, SIMD instructions perform contiguous memory accesses, improving cache locality and prefetching behavior.

Within the XMM/YMM kernel, two key operations make vector reduction efficient:

`vextractf128` extracts the upper 128-bit lane of a 256-bit YMM register, enabling the combination of the upper and lower halves of the vector directly in registers.

`vhaddps` performs horizontal addition inside each vector, summing elements pairwise without scalar loops.
By chaining two `vhaddps` instructions, the kernel quickly reduces eight partial sums to a single scalar value, avoiding costly memory stores or extra loop passes.

Overall, these operations allow the AVX2 version to fully exploit data-level parallelism while keeping all arithmetic within registers — leading to a 10+× speedup over the scalar baseline, depending on matrix size and memory alignment.

`Cache Efficiency` SIMD kernels reuse contiguous memory regions effectively.Streaming access patterns benefit from hardware prefetchers, minimizing cache misses.

---

## Troubleshooting

### What happened with MEMCMP not equal but EPSILON works?

In this project we first tried to verify correctness by comparing the output vectors byte-by-byte:

```C
	memcmp(Y_c, Y_simd, sizeof(float)*n)
```

This failed even when the printed numbers looked the same. The reason is that `floating-point` results depend on accumulation/order and rounding

Since `memcmp` does a strict byte-for-byte comparison, any tiny rounding difference makes it report `not equal`.

Fix: compare numerically with a `epsilon` instead of bitwise equality. We treat two floats as equal if their difference is within a small absolute/relative bound:

```C
static inline boolean nearly_equal(float a, float b, float eps) {
	float diff = fabsf(a - b);
	float scale = fmaxf(1.0f, fmaxf(fabsf(a), fabsf(b)));
	return diff <= eps * scale;
}
```

With this epsilon check, `C` vs. `XMM` and `C` vs. `YMM` both pass.
