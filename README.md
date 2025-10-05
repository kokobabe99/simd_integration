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
     │   ├── asmfunc2.asm          (SIMD XMM)
     │   ├── asmfunc3.asm          (SIMD YMM)
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
void matvec_c(int n, float* A, float* x, float* y) {
	for (int i = 0; i < n; ++i) {
		double s = 0.0;
		float* row = A + (size_t)i * n;
		for (int j = 0; j < n; ++j)
			s += (double)row[j] * (double)x[j];
		y[i] = (float)s;
	}
}
```

### Results

```
Matrix sizes: n ∈ {1003,2^10,2^13,2^14} RUN in 30 times

#define RUNS 10
```

| n            |      C (ms) |   x86 (ms) |  XMM (ms) |  YMM (ms) | x86 Speedup VS C | XMM Speedup VS C | YMM Speedup VS C |
| ------------ | ----------: | ---------: | --------: | --------: | ---------------: | ---------------: | ---------------: |
| 1003         |     521.193 |    166.493 |    28.962 |    18.537 |            3.13× |           18.00× |           28.12× |
| 2^10 (1024)  |     411.454 |    122.419 |    53.028 |    32.243 |            3.36× |            7.76× |           12.76× |
| 2^13 (8192)  |  25,186.700 |  7,837.390 | 2,255.319 | 1,971.116 |            3.21× |           11.17× |           12.78× |
| 2^14 (16384) | 101,722.746 | 31,554.394 | 9,077.505 | 7,961.411 |            3.22× |           11.21× |           12.78× |

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
		printf("last 3 for Correctness check");
		print_first_last("kernel (C  ):", Y_c, n);
		print_first_last("kernel (x86):", Y_x86_64, n);
		print_first_last("kernel (xmm):", y_xmm_v2, n);
		print_first_last("kernel (ymm)", Y_ymm, n);
```

## Comparative Results

---

## Analysis

Why SIMD is faster: data-level parallelism (4-wide/8-wide), reduced loop overhead, contiguous loads.

XMM vs YMM: wider vectors, reduction cost, memory-bandwidth limits.

Impact of boundary handling: tail fraction shrinks as n grows.

Effects of alignment (vmovups vs vmovaps) and AVX↔SSE transition (mitigated with VEX forms and vzeroupper).

Floating-point summation order → bitwise differences; hence tolerance-based checks.

---

### What happened with MEMCMP ?
