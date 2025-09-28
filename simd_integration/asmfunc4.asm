

; y = A * x (single-precision float) 
; AVX2 / XMM knernel (128-bit) —4 * float per loop 32-bit
; RCX = n (int) 
; RDX = A (float*) 
; R8 = x (float*) 
; R9 = y (float*)

default rel
bits 64
section .text
global matvec_xmm_avx2_asm_v2

matvec_xmm_avx2_asm_v2:

    push r12
    push rdi
    push rsi
    xor rsi, rsi ;outer loop counter i

loop_outer:
    cmp rsi, rcx
    jge done
    vxorps xmm0, xmm0,xmm0;reset result per row   

    mov r10, rsi
    imul r10, rcx
    lea r11, [rdx + r10 * 4] ;row (&A[i])

    mov r12d, ecx
    and r12d, -4   ; align to 4
    xor edi, edi ;inner loop counter j

loop_inner:

    cmp edi, r12d
    jge boundary_check

    vmovups xmm1, [r11 + rdi * 4]
    vmovups xmm2, [r8 + rdi * 4]
    vmulps xmm1, xmm1, xmm2
    vaddps xmm0, xmm0, xmm1
    
    add edi, 4
    jmp loop_inner


 boundary_check:
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
    ret     
   
   
