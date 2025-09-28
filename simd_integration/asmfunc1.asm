

; y = A * x (single-precision float) 
; x86-64 kernel
; RCX = n (int) 
; RDX = A (float*) 
; R8 = x (float*) 
; R9 = y (float*)

default rel
bits 64
section .text
global matvec_x86_64_asm

matvec_x86_64_asm:
    xor rsi, rsi ;outer loop counter i

loop_outer:
    xorps xmm0, xmm0 ;result per row    
    mov r10, rsi
    imul r10, rcx
    lea r11, [rdx + r10 * 4] ;row (&A[i])
    xor rdi, rdi ;inner loop counter j

loop_inner:
    movss xmm1, [r11 + rdi * 4]
    mulss xmm1, [r8 + rdi * 4]
    addss xmm0, xmm1
    
    inc rdi
    cmp rdi, rcx
    jl loop_inner

store:
    movss [r9 + rsi * 4], xmm0
    
    inc rsi
    cmp rsi, rcx
    jl loop_outer

done:
    ret     
   
