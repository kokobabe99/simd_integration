

; y = A * x (single-precision float) 
; AVX2 / ymm knernel (256-bit) —8 * float per loop 32-bit
; RCX = n (int) 
; RDX = A (float*) 
; R8 = x (float*) 
; R9 = y (float*)

default rel
bits 64
section .text
global matvec_ymm_avx2_asm

matvec_ymm_avx2_asm:

    push r12
    push rdi
    push rsi
    xor rsi, rsi ;outer loop counter i

loop_outer:
    cmp     esi, ecx
    jge     done
    vxorps  ymm0, ymm0, ymm0 

    mov     r10d, esi
    imul    r10d, ecx
    lea     r11, [rdx + r10*4]           

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
store:
    vmovss dword [r9 + rsi * 4], xmm3
    inc esi
    jmp loop_outer
    
done:
    pop rsi
    pop rdi
    pop r12
    vzeroupper
    ret     
   
   
