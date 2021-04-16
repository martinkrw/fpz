; Copyright (c) 2021 Martin Wanvik <martin.kr.wanvik@gmail.com>
;
; Permission to use, copy, modify, and distribute this software for any
; purpose with or without fee is hereby granted, provided that the above
; copyright notice and this permission notice appear in all copies.
;
; THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
; WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
; MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
; ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
; WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
; ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
; OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

%include "asmmacros.inc"

%macro MUL_64x64_ITER 1
	MUL_64x64 [DST + %1*8], TMP0, OP2_0, [OP1 + %1*8], TMP0
%endmacro

%macro MUL_64x128_ITER 1
	mov OP1_0, [OP1 + %1*8]

	MUL_64x64 [DST + %1*8], TMP0, OP1_0, OP2_0, TMP0
	MUL_64x64 TMP0, TMP1, OP1_0, OP2_1, TMP0, TMP1
%endmacro

%macro MUL_64x192_ITER 1
	mov OP1_0, [OP1 + %1*8]

	MUL_64x64 [DST + %1*8], TMP0, OP1_0, OP2_0, TMP0
	MUL_64x64 TMP0, TMP1, OP1_0, OP2_1, TMP0, TMP1
	MUL_64x64 TMP1, TMP2, OP1_0, OP2_2, TMP1, TMP2
%endmacro

%macro MUL_64x256_ITER 1
	mov OP1_0, [OP1 + %1*8]

	MUL_64x64 [DST + %1*8], TMP0, OP1_0, OP2_0, TMP0
	MUL_64x64 TMP0, TMP1, OP1_0, OP2_1, TMP0, TMP1
	MUL_64x64 TMP1, TMP2, OP1_0, OP2_2, TMP1, TMP2
	MUL_64x64 TMP2, TMP3, OP1_0, OP2_3, TMP2, TMP3
%endmacro

%macro ADDMUL_64x256_ITER 1
	mov OP1_0, [OP1 + %1*8]

	MUL_64x64 [DST + %1*8], TMP0, OP1_0, OP2_0, TMP0, [DST + %1*8]
	MUL_64x64 TMP0, TMP1, OP1_0, OP2_1, TMP0, TMP1
	MUL_64x64 TMP1, TMP2, OP1_0, OP2_2, TMP1, TMP2
	MUL_64x64 TMP2, TMP3, OP1_0, OP2_3, TMP2, TMP3
%endmacro

%macro COMPUTE_ROW 1
%push
%define %$ITERATION 	%1
	mov rax, N1_TMP		
				
	lea INDEX1, [rax - 1]
	neg INDEX1

%%loop:
	%$ITERATION INDEX1
	
	add INDEX1, 1
	jle %%loop
%%end:
%pop
%endmacro

%define FACTOR1	rdi
%define FACTOR2 rsi
%define TERM1	r8	; passed in rdx, but will be moved
%define TERM2	rcx

; ================================================================================
; u128 fpz_addmul_u64_u64(digit f1, digit f2, digit term1, term2)
; 
; Compute and return f1*f2 + term.
; ================================================================================

global fpz_addmul_u64_u64
fpz_addmul_u64_u64:

	mov TERM1, rdx

	MUL_64x64 rax, rdx, FACTOR1, FACTOR2, TERM1, TERM2

	ret


;================================================================================
; digit fpz_addmul_u64(dstptr dst, digit factor, digit term, digcnt dst_cap)
;
; Compute (*dst) * factor + term and store the result in *dst again. Return the 
; most significant digit in case the result overflows the given capacity. It 
; does not set or check the overflow flag; this is used when setting an integer 
; value from a string, something which requires repeated multiplication by 10/16 
; followed by adding a digit value. In that context, overflow is expected 
; (i.e. it is not an indication of a bug in the library), and using the overflow 
; flag to check for this is cumbersome.
;================================================================================

%define DST	rdi
%define OP1	rdi
%define FACTOR	rsi
%define TERM	rdx
%define DST_CAP	rcx

; Registers used for header

%define HDR		r8

; Registers used to hold digit counts for each operand (same as header registers, but 
; used only after a right shift of 32 bits 

%define N1		r8
%define N1_TMP		r8

%define INDEX1		r9	; overlap with N2 

%define OP2_0		r10
%define TMP0		r11

global fpz_addmul_u64_self
fpz_addmul_u64_self:
	xor eax, eax
	mov HDR, HEADER(DST)

	test FACTOR, FACTOR
	setnz al

	shr HDR, 32
	setnz ah

	test al, ah	; Check if the product is going to be zero
	jz .product_zero

; Adjust pointers so that we may index using negative numbers running up to zero, rather
; than positive integers running up to N1 (this is more convenient when there is a shortage 
; of registers, since it means that we only need to read the digit count at the loop entry
; rather than after every iteration; we do it here since it would be more of a hassle not
; to do it (and consistency is good)).

	lea DST, [DST + N1*8]

	mov TMP0, TERM			; make sure TERM gets added into the result
	
	COMPUTE_ROW MUL_64x64_ITER
	
	mov rax, TMP0	; return value will be the most significant digit if there isn't room to write it

; reset the value of DST so that it points to the header

	mov rdx, N1
	neg rdx
	lea DST, [DST + rdx*8]

	cmp DST_CAP(d), N1(d)
	jbe .done

	mov [DST + N1*8 + 8], rax

	neg rax
	adc N1(d), 1

	xor eax, eax

.done:

	mov [DST + size], N1(d)

	ret

.product_zero:
; Set the value of *dst to TERM

	xor eax, eax

	test TERM, TERM
	setnz al
	shl rax, 32

	mov qword [DST], rax		; set sign to positive, overflow flag to 0 and size to (TERM != 0)
	mov qword DIGIT(DST, 0), TERM

	ret

; ==============================================================================
; void fpz_mul_i64_i64(dstptr dst, int64_t op1, int64_t op2, digcnt dst_cap)
;
; Multiply two 64-bit (built-in, signed two's complement) integers and store 
; the result in *dst. Will write one digit unconditionally (even if op1*op2 = 0), 
; and a second one if required and if dst_cap is big enough (sets overflow flag
; if dst_cap is not sufficient).
; ==============================================================================

section .text

%define DST	rdi
%define OP1	rsi
%define OP2	rdx
%define DST_CAP	rcx

global fpz_mul_i64_i64
fpz_mul_i64_i64:

	SIGN_MAGNITUDE OP1, r9, OP1
	SIGN_MAGNITUDE OP2, r8, OP2

	xor r8d, r9d
	and r8d, 1

	mov rax, OP1
	mul OP2

	jnc .one_digit

	cmp DST_CAP(d), 2
	jb overflow

	mov rcx, 0x0000000200000000
	or rcx, r8

	mov HEADER(DST), rcx
	mov DIGIT(DST, 0), rax
	mov DIGIT(DST, 1), rdx

	ret

.one_digit:
	
	xor ecx, ecx

	test rax, rax
	setnz cl

	shl rcx, 32
	or rcx, r8

	mov HEADER(DST), rcx
	mov DIGIT(DST, 0), rax

	ret




; ==============================================================================
; void fpz_mul_i64(dstptr dst, srcptr op1, int64_t op2, digcnt dst_cap)
; void fpz_mul(dstptr dst, srcptr op1, srcptr op2, digcnt dst_cap)
;
; Compute the product op1 and op2 and store the result in *dst, which is assumed
; to have capacity for dst_cap digits. If this is insufficient, the overflow
; flag is set. If any of the non built-in integers have their overflow flag set,
; print an error message and abort the entire process (by calling 
; fpz_use_after_overflow_handler)
; ==============================================================================

; Registers used for header

%define HDR1		r8

; Registers used to hold digit counts for each operand (same as header registers, but 
; used only after a right shift of 32 bits 

%define N1		r8
%define N2		r9

%define INDEX1		r9	; overlap with N2 
%define INDEX2		r8	; overlap with N1

%define OP2_0		r10
%define TMP0		r11

%define OP1_0		r12

%define OP2_1		r13
%define TMP1		r14

%define OP2_2		r15
%define TMP2		rbx

%define OP2_3		rbp
%define TMP3		rcx	; overlap with DST_CAP


%define N1_TMP	r8	; Same as where it originally is, so no need to move

global fpz_mul_i64
fpz_mul_i64:

	mov HDR1, HEADER(OP1)

	CHECK_OVERFLOW HDR1, use_after_overflow
	
; Magnitude gets stored in OP2_0, freeing OP2 (rdx) for other uses

	SIGN_MAGNITUDE OP2, rax, OP2_0	
	setnz dl			; The zero flag will get set if (and only if) OP2 = 0

; Compute (and write) sign of result

	xor eax, HDR1(d)	
	and eax, 1
	
	mov [DST], eax

	shr HDR1, 32
	setnz dh

	test dl, dh	; Check if either of the operands are zero
	jz .zero

; Since the result is nonzero, it requires at least N1 digits. If we were not
; given enough digits, set overflow flag in DST and return

	cmp DST_CAP(d), N1(d)
	jb overflow

; Adjust pointers so that we may index using negative numbers running up to zero, rather
; than positive integers running up to N1 (this is more convenient when there is a shortage 
; of registers, since it means that we only need to read the digit count at the loop entry
; rather than after every iteration; we do it here since it would be more of a hassle not
; to do it (and consistency is good)).

	lea DST, [DST + N1*8]
	lea OP1, [OP1 + N1*8]

	xor TMP0(d), TMP0(d)
	
	COMPUTE_ROW MUL_64x64_ITER
	
	mov rax, TMP0

	add DST, 8	
	mov N2(d), 1
	
	jmp common_epilogue

.zero:

	mov qword [DST], 0
	mov qword DIGIT(DST, 0), 0

	ret

overflow:

	mov qword HEADER(DST), 0x00000100	; set size to 0, sign to positive and overflow_flag to 1

	ret

use_after_overflow:

	USE_AFTER_OVERFLOW_HANDLER

	ret

; ================================================================================
; void fpz_umul_sp(fpz_dstptr dst, fpz_srcptr op1, fpz_srcptr op2, 
;	fpz_digcnt dst_cap, fpz_digcnt n1, fpz_digcnt n2);
; void fpz_umul_mp(fpz_dstptr dst, fpz_srcptr op1, fpz_srcptr op2,
;	fpz_digcnt dst_cap, fpz_digcnt n1, fpz_digcnt n2);
;
; Perform an unsigned multiplication of *op1 and *op2, pointing to arrays of 
; n1 and n2 digits, respectively, and store result in dst (which is assumed to have 
; space for dst_cap digits). 
; The function requires that n1 >= n2. fpz_mul_sp requires that n2 <= 4, and performs 
; the computation in a single pass (hence the name). This is advantageous if 
; op1 and/or op2 coincides with dst, and allows us to skip copying of such operands to 
; the stack. Note that mul_mp doesn't test for this, this is supposed to be done in the 
; caller (these functions are called from fpz_mul, which is written in C++).
; ================================================================================

align 16
global fpz_umul_sp
global fpz_umul_mp

fpz_umul_sp:

	lea DST, [DST + N1*8]
	lea OP1, [OP1 + N1*8]
	lea OP2, [OP2 + N2*8]

	lea rax, [rel jmp_table_sp]
	jmp [rax + N2*8 - 8]

umul_64_sp:

	mov OP2_0, [OP2]
	xor TMP0(d), TMP0(d)

	COMPUTE_ROW MUL_64x64_ITER

	mov rax, TMP0

	add DST, 8
	mov N2(d), 1

	jmp common_epilogue

umul_128_sp:

	WRITE_TO_REDZONE r12, r13, r14

	mov OP2_0, [OP2 - 8]
	mov OP2_1, [OP2]

	xor TMP0(d), TMP0(d)
	xor TMP1(d), TMP1(d)

	COMPUTE_ROW MUL_64x128_ITER

	mov [DST + 8], TMP0
	
	mov rax, TMP1

	add DST, 16
	mov N2(d), 2

	READ_FROM_REDZONE r12, r13, r14

	jmp common_epilogue

umul_192_sp:

	WRITE_TO_REDZONE r12, r13, r14, r15, rbx

	mov OP2_0, [OP2 - 16]
	mov OP2_1, [OP2 - 8]
	mov OP2_2, [OP2]

	xor TMP0(d), TMP0(d)
	xor TMP1(d), TMP1(d)
	xor TMP2(d), TMP2(d)

	COMPUTE_ROW MUL_64x192_ITER

	mov [DST + 8], TMP0
	mov [DST + 16], TMP1

	mov rax, TMP2

	add DST, 24
	mov N2(d), 3

	READ_FROM_REDZONE r12, r13, r14, r15, rbx

	jmp common_epilogue

umul_256_sp:

	WRITE_TO_REDZONE DST_CAP, r12, r13, r14, r15, rbx, rbp

	mov OP2_0, [OP2 - 24]
	mov OP2_1, [OP2 - 16]
	mov OP2_2, [OP2 - 8]
	mov OP2_3, [OP2]

	xor TMP0(d), TMP0(d)
	xor TMP1(d), TMP1(d)
	xor TMP2(d), TMP2(d)
	xor TMP3(d), TMP3(d)

	COMPUTE_ROW MUL_64x256_ITER

	mov [DST + 8], TMP0
	mov [DST + 16], TMP1
	mov [DST + 24], TMP2

	mov rax, TMP3

	add DST, 32
	mov N2(d), 4

	READ_FROM_REDZONE DST_CAP, r12, r13, r14, r15, rbx, rbp

	jmp common_epilogue

align 16
fpz_umul_mp:

	lea DST, [DST + N1*8]
	lea OP1, [OP1 + N1*8]
	lea OP2, [OP2 + N2*8]

%define OP2_TMP		[rsp - 8]
%define N1_TMP		[rsp - 16]
%define N2_TMP		[rsp - 24]

	WRITE_TO_REDZONE OP2, N1, N2, DST_CAP, r12, r13, r14, r15, rbx, rbp

; Set up the value of INDEX2 (this will destroy the value of N1, but that is saved now)

	lea INDEX2, [N2 - 1]
	neg INDEX2

	lea rax, [rel jmp_table_mp]
	and N2, 3		; destroys N2, but that is saved (and won't be used again for a long time)
	jmp [rax + N2*8]

umul_64_mp:

	mov OP2_0, [OP2 + INDEX2*8]
	xor TMP0(d), TMP0(d)

	COMPUTE_ROW MUL_64x64_ITER

	mov OP2, OP2_TMP	; retrieve second operand from the stack again

	mov rax, TMP0

	add DST, 8
	add INDEX2, 1

	jmp main_loop

umul_128_mp:

	mov OP2_0, [OP2 + INDEX2*8]
	mov OP2_1, [OP2 + INDEX2*8 + 8]

	xor TMP0(d), TMP0(d)
	xor TMP1(d), TMP1(d)

	COMPUTE_ROW MUL_64x128_ITER

	mov OP2, OP2_TMP	; retrieve second operand again

	mov [DST + 8], TMP0
	mov rax, TMP1
	
	add DST, 16
	add INDEX2, 2

	jmp main_loop

umul_192_mp:

	mov OP2_0, [OP2 + INDEX2*8]
	mov OP2_1, [OP2 + INDEX2*8 + 8]
	mov OP2_2, [OP2 + INDEX2*8 + 16]

	xor TMP0(d), TMP0(d)
	xor TMP1(d), TMP1(d)
	xor TMP2(d), TMP2(d)

	COMPUTE_ROW MUL_64x192_ITER

	mov OP2, OP2_TMP	; retrieve second operand again

	mov [DST + 8], TMP0
	mov [DST + 16], TMP1

	mov rax, TMP2

	add DST, 24
	add INDEX2, 3

	jmp main_loop

umul_256_mp:

	mov OP2_0, [OP2 + INDEX2*8]
	mov OP2_1, [OP2 + INDEX2*8 + 8]
	mov OP2_2, [OP2 + INDEX2*8 + 16]
	mov OP2_3, [OP2 + INDEX2*8 + 24]

	xor TMP0(d), TMP0(d)
	xor TMP1(d), TMP1(d)
	xor TMP2(d), TMP2(d)
	xor TMP3(d), TMP3(d)

	COMPUTE_ROW MUL_64x256_ITER

	mov OP2, OP2_TMP	; retrieve second operand again

	mov [DST + 8], TMP0
	mov [DST + 16], TMP1
	mov [DST + 24], TMP2

	mov rax, TMP3

	add DST, 32
	add INDEX2, 4

	jmp main_loop

align 16
main_loop:

	mov [DST], rax	; store most significant digit from previous iteration

	mov OP2_0, [OP2 + INDEX2*8]
	mov OP2_1, [OP2 + INDEX2*8 + 8]
	mov OP2_2, [OP2 + INDEX2*8 + 16]
	mov OP2_3, [OP2 + INDEX2*8 + 24]

	xor TMP0(d), TMP0(d)
	xor TMP1(d), TMP1(d)
	xor TMP2(d), TMP2(d)
	xor TMP3(d), TMP3(d)

	COMPUTE_ROW ADDMUL_64x256_ITER

	mov OP2, OP2_TMP	; retrieve second operand again

	mov [DST + 8], TMP0
	mov [DST + 16], TMP1
	mov [DST + 24], TMP2

	mov rax, TMP3

	add DST, 32

	add INDEX2, 4
	jl main_loop

; Note: first argument empty, so skips reading into OP2 (it will not be needed)

	READ_FROM_REDZONE NOWHERE, N1, N2, DST_CAP, r12, r13, r14, r15, rbx, rbp

common_epilogue:

; Compute the total number of digits written up to this point 

	lea N1(d), [N1(d) + N2(d) - 1]
	
; Reset the destination pointer so that it points to the header again

	mov rdx, N1
	neg rdx
	lea DST, [DST + rdx*8 - 8]

; Figure out if we have room to write the final digit

	cmp DST_CAP(d), N1(d)
	jbe .done

	mov [DST + N1*8 + 8], rax

; Increment digit count if final digit written is nonzero

	neg rax
	adc N1(d), 0
	
	xor eax, eax

.done:
	
	test rax, rax
	setnz [DST + overflow_flag]
	mov [DST + size], N1(d)

	ret

section .data

; Jump tables (single and multi-pass multiplication)

jmp_table_sp:
	dq umul_64_sp
	dq umul_128_sp
	dq umul_192_sp
	dq umul_256_sp

jmp_table_mp:

	dq umul_256_mp	
	dq umul_64_mp	
	dq umul_128_mp
	dq umul_192_mp

section .text

;================================================================================
; fpz_digit fpz_usubmul_lowcap_64_self(fpz_dstptr dst, fpz_srcptr op1, fpz_digit op2, fpz_digcnt n);
; fpz_digit fpz_usubmul_64_self(fpz_dstptr dst, fpz_srcptr op1, fpz_digit op2, fpz_digcnt n);
; 
; Internal functions used with division, working with unsigned parameters. The parameters 
; dst and op1 are pointers to arrays of digits, not to the header of a full integer. 
;
; Both functions compute *dst -= op1*op2, where op1 is a pointer to n digits.
; The difference is that the first function assumes that dst also contains n digits,
; while the latter assumes it contains n + 1 digits. In the second case, we perform
; a final subtraction from the most significant digit of dst, while in the first
; case we do not (and simply return the sum of the propagated borrow and digit
; n+1 of op1*op2).
;================================================================================

%define DST	rdi
%define OP1	rsi
%define OP2	rdx
%define N	rcx
%define N1	rcx

%define N1_TMP	rcx	; Same as where it originally is, so no need to move

%define INDEX	r8
%define INDEX1	r8
%define OP2_0	r9
%define TMP0	r10

%macro SUBMUL_64x64_ITER 1
	MUL_64x64 rax, TMP0, OP2_0, [OP1 + %1*8], TMP0

	sub [DST + %1*8], rax
	adc TMP0, 0
%endmacro

global fpz_usubmul_lowcap_64_self
global fpz_usubmul_64_self

fpz_usubmul_lowcap_64_self:

	mov OP2_0, OP2

	lea DST, [DST + N*8 - 8]
	lea OP1, [OP1 + N*8 - 8]

	xor TMP0(d), TMP0(d)
	
	COMPUTE_ROW SUBMUL_64x64_ITER
	
	mov rax, TMP0

	ret

fpz_usubmul_64_self:

	mov OP2_0, OP2

	lea DST, [DST + N*8 - 8]
	lea OP1, [OP1 + N*8 - 8]

	xor TMP0(d), TMP0(d)
	
	COMPUTE_ROW SUBMUL_64x64_ITER
	
	xor eax, eax
	sub [DST + INDEX*8], TMP0
	setc al

	ret
	
