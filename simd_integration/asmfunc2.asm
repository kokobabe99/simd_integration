

; y = A * x (single-precision float) 
; AVX2 / XMM knernel (128-bit) —4 * float per loop 32-bit
; RCX = n (int) 
; RDX = A (float*) 
; R8 = x (float*) 
; R9 = y (float*)

default rel
bits 64
section .text
global matvec_xmm_avx2_asm

matvec_xmm_avx2_asm:
    xor rsi, rsi ;outer loop counter i

loop_outer:
    xorps xmm0, xmm0 ;result per row    
    mov r10, rsi
    imul r10, rcx
    lea r11, [rdx + r10 * 4] ;row (&A[i])
    xor rdi, rdi ;inner loop counter j

loop_inner:
    vmovdqu xmm1, [r11 + rdi * 4]
    vmovdqu xmm2, [r8 + rdi * 4]
    vmulps xmm1, xmm1, xmm2
    vaddps xmm0, xmm0, xmm1
    
    add rdi, 4
    cmp rdi, rcx
    jl loop_inner

store:
    vhaddps xmm0, xmm0
    vhaddps xmm0, xmm0
    vmovss [r9 + rsi * 4], xmm0
    
    inc rsi
    cmp rsi, rcx
    jl loop_outer

done:
    ret     
   
   
