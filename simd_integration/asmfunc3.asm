

; y = A * x (single-precision float) 
; AVX2 / YMM knernel (256-bit) —8 * float per loop 32-bit,
; RCX = n (int) 
; RDX = A (float*) 
; R8 = x (float*) 
; R9 = y (float*)

default rel
bits 64
section .text
global matvec_ymm_avx2_asm

matvec_ymm_avx2_asm:

     
    sub     rsp, 32
    xor     r10d, r10d  
loop_outer:
    cmp     r10d, ecx      
    jge     done
    mov     eax, r10d   ; row = &A[i,0] = A + i*n*4
    imul    eax, ecx
    lea     r11, [rdx + rax*4]        

    mov     eax, ecx
    and     eax, -8
    mov     dword [rsp+4], eax ; save the reminder of n * 8 in the stack

    vxorps  ymm0, ymm0, ymm0
    mov     dword [rsp+0], 0     

loop_inner:
    mov     eax, dword [rsp+0]        
    cmp     eax, dword [rsp+4] ; check if there is remaining out of 8
    jge     boundary_check

    vmovups ymm1, [r11 + rax*4]        
    vmovups ymm2, [r8  + rax*4]        
    vmulps  ymm1, ymm1, ymm2
    vaddps  ymm0, ymm0, ymm1

    add     eax, 8
    mov     dword [rsp+0], eax         
    jmp     loop_inner

boundary_check:
    vextractf128 xmm1, ymm0, 1
    vaddps      xmm0, xmm0, xmm1      
    vhaddps     xmm0, xmm0, xmm0       
    vhaddps     xmm0, xmm0, xmm0        
    vmovaps     xmm3, xmm0              

    mov     eax, dword [rsp+4]        
reminder_add:
    cmp     eax, ecx           
    jge     store

    vmovss   xmm1, dword [r11 + rax*4]
    vmulss   xmm1, xmm1, dword [r8 + rax*4]
    vaddss   xmm3, xmm3, xmm1 
    inc     eax
    jmp     reminder_add
store:
    vmovss  dword [r9 + r10*4], xmm3
    inc     r10d
    jmp     loop_outer

done:
    add     rsp, 32
    ret
