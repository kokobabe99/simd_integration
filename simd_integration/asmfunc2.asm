

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

     
   
