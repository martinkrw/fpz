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

; ADD INDEX, ACTIONS
;
; Iteration macro for adding digit INDEX of OP1 to the same digit of OP2 and 
; storing the result into digit INDEX of RES. 

%macro ADD 2
%if %2 & RESTORE
	neg al
%endif
	mov TMP, DIGIT(OP1, %1)
%ifidn %1, 0
	add TMP, DIGIT(OP2, %1)
%else
	adc TMP, DIGIT(OP2, %1)
%endif
	mov DIGIT(DST, %1), TMP
%if %2 & SAVE
	setc al
%endif
%endmacro

; ADD_CARRY INDEX, ACTIONS
;
; Iteration macro for adding the carry (actually, whatever happens to be stored
; in rax) to digit INDEX of OP1 and then storing the result to digit INDEX of
; RES. 

%macro ADD_CARRY 2
	mov TMP, DIGIT(OP1, %1)
%if %2 & RESTORE
	add TMP, rax
%else
	adc TMP, 0
%endif
	mov DIGIT(DST, %1), TMP
%if %2 & SAVE
	setc al
%endif
%endmacro

; SUB INDEX, ACTIONS
; 
; Iteration macro for subtracting digit INDEX of OP2 from digit INDEX of OP1
; and storing the result into digit INDEX of RES. 

%macro SUB 2
%if %2 & RESTORE
	neg al
%endif
	mov TMP, DIGIT(OP1, %1)
%ifidn %1, 0
	sub TMP, DIGIT(OP2, %1)
%else
	sbb TMP, DIGIT(OP2, %1)
%endif
	mov DIGIT(DST, %1), TMP
%if %2 & SAVE
	setc al
%endif
%endmacro

; Same as ADD_CARRY, except with subtraction and borrow. 

%macro SUB_BORROW 2
	mov TMP, DIGIT(OP1, %1)
%if %2 & RESTORE
	sub TMP, rax
%else
	sbb TMP, 0
%endif
	mov DIGIT(DST, %1), TMP
%if %2 & SAVE
	setc al
%endif
%endmacro

; Loop for adding/subtracting a single digit to/from a larger integer.
; Assumes that digit count (N) is at least 2

%define VARIANT_ADD	1
%define VARIANT_SUB	2

%macro ADD_SUB_64 1
%push
%define %$VARIANT	%1

	xor eax, eax
	mov INDEX(d), 1

%if %$VARIANT == VARIANT_ADD
	add OP2, DIGIT(OP1, 0)
	mov DIGIT(DST, 0), OP2
%else
	mov TMP, DIGIT(OP1, 0)
	sub TMP, OP2
	mov DIGIT(DST, 0), TMP
%endif

	setc al	

	test N, 1
	jnz .loop	; if N is odd, an even number of digits remain, so jump into the loop

	mov INDEX(d), 2

	mov TMP, DIGIT(OP1, 1)
%if %$VARIANT == VARIANT_ADD
	add TMP, rax
%else
	sub TMP, rax
%endif
	mov DIGIT(DST, 1), TMP

	setc al
	
	cmp N, 2
	jbe .end

.loop:

	mov TMP, DIGIT(OP1, INDEX)
	mov OP2, DIGIT(OP1, INDEX + 1)

%if %$VARIANT == VARIANT_ADD
	add TMP, rax
	adc OP2, 0
%else
	sub TMP, rax
	sbb OP2, 0
%endif

	setc al

	mov DIGIT(DST, INDEX), TMP
	mov DIGIT(DST, INDEX + 1), OP2

	add INDEX, 2
	cmp INDEX, N

	jb .loop

.end:
%pop
%endmacro
	
; COUNT_DIGITS
;
; Correct the digit count for the result (used only for subtraction, where 
; a lot of the result digits may end up being zero). 

%macro COUNT_DIGITS 0
%%count_digits_loop:
	cmp qword DIGIT(DST, N1 - 1), 0
	jne %%done

	sub N1, 1
	ja %%count_digits_loop
%%done:
%endmacro

; ================================================================================
; void	fpz_add_i64_i64(dstptr dst, int64_t op1, int64_t op2, digcnt dst_cap);
; void	fpz_sub_i64_i64(dstptr dst, int64_t op1, int64_t op2);
;
; Add/subtract two built-in 64-bit signed integers (op1 and op2) and write the 
; result into dst, which is assumed to have storage for dst_cap digits
; (subtraction can't overflow in this case, so no need to pass that particular 
; parameter).
; ================================================================================

; Function parameters

%define DST		rdi
%define OP1		rsi
%define OP2		rdx
%define DST_CAP		rcx

%define SIGN		r8
%define TMP		r9

section .text

global fpz_add_i64_i64
global fpz_sub_i64_i64

fpz_add_i64_i64:

	xor eax, eax
	xor r8d, r8d
	xor r9d, r9d

	add OP1, OP2

	setnz al

	jge .ok 

; Decide if the result needs one or two digits. The only case where we
; need two digits is if both operands have the value -2^63 (0x8000000000000000), 
; in which case we have of = 1, sf = 0 and zf = 1. Since we're here, we know that
; of != sf, so it suffices that (of and zf) = 1  or (!of or !zf) = 0. 

	setno TMP(b)
	or TMP(b), al

	jz .two_digits	
	
	mov SIGN(d), 1	; set sign to negative
	neg OP1		; negate result

.ok:

	shl rax, 32	; shift digit count to its correct place in header
	or rax, SIGN	; set the sign bit
	
	mov [DST], rax
	mov DIGIT(DST, 0), OP1
	
	xor eax, eax

	ret

.two_digits:

	cmp DST_CAP, 2
	jb overflow

; Header containing a digit count of 2 and a negative sign 

	mov rax, 0x0000000200000001

	mov HEADER(DST), rax
	mov qword DIGIT(DST, 0), 0
	mov qword DIGIT(DST, 1), 1

	xor eax, eax

	ret

; Subtraction is simpler, since there is no way for the result to require 
; two digits

fpz_sub_i64_i64:

	xor eax, eax
	xor SIGN(d), SIGN(d)

	sub OP1, OP2

	setnz al

	jge .ok

	mov SIGN(d), 1
	neg OP1
.ok:
	shl rax, 32
	or rax, SIGN

	mov [DST], rax
	mov DIGIT(DST, 0), OP1

	ret

; ================================================================================
; void fpz_add_i64(dstptr dst, srcptr op1, int64_t op2, digcnt dst_cap);
; void fpz_sub_i64(dstptr dst, srcptr op1, int64_t op2, digcnt dst_cap);
; void fpz_i64_sub(dstptr dst, int64_t op1, srcptr op2, digcnt dst_cap);
;
; Compute the sum/difference between the integers op1 and op2 and store the result
; in dst, which is assumed to have space for dst_cap digits. Set overflow flag if this 
; is insufficient. If any of non-built in integer operands have the overflow flag 
; set, print error message and abort the entire process (by calling 
; fpz_use_after_overflow_handler).
; ================================================================================

; Function parameters

%define DST	rdi
%define OP1	rsi
%define OP2	rdx	
%define DST_CAP	rcx	; capacity of result

; Variables

%define HDR	r8	; header of first fpz integer operand (deliberately coinciding with N)	
%define HDR1	r8	; alternate name
%define N	r8	; digit count of first fpz integer operand
%define N1	r8	; alternate name
%define HDR2	r9	; header of second fpz integer operand
%define N2	r9	; digit count of second fpz integer operand

%define INDEX		r10
%define TMP		r11

global fpz_add_i64
global fpz_sub_i64
global fpz_i64_sub

align 16
fpz_add_i64:

	mov HDR, HEADER(OP1)
	mov eax, HDR(d)

	CHECK_OVERFLOW HDR, use_after_overflow

	shr HDR, 32	; get digit count into N

	SIGN_MAGNITUDE OP2, TMP, OP2
	and TMP(d), 1		; reduce sign of OP2 to 1 bit

	xor TMP(d), eax		; compare signs

	jnz subtract_64
	jmp add_64

fpz_sub_i64:

	mov HDR, HEADER(OP1)
	mov eax, HDR(d)

	CHECK_OVERFLOW HDR, use_after_overflow

	shr HDR, 32		; get digit count into N

	SIGN_MAGNITUDE OP2, TMP, OP2
	and TMP(d), 1

	xor TMP(d), eax		; result of sign comparison in eax

	jz subtract_64
	jmp add_64

; This is needed because subtraction isn't commutative (so we can't rewrite
; i - x to x - i). Performed like an addition with arguments swapped and the 
; first operand negated, i.e. using the identity i - x = -x + i.

fpz_i64_sub:

	mov HDR(d), 1

; exchange arguments

	mov TMP, OP1
	mov OP1, OP2
	mov OP2, TMP

	xor HDR, HEADER(OP1)	; flip the sign of the first (previously second) operand

	CHECK_OVERFLOW HDR, use_after_overflow

	mov eax, HDR(d)

	shr HDR, 32

	SIGN_MAGNITUDE OP2, TMP, OP2
	and TMP(d), 1
	
	xor TMP(d), eax

	jnz subtract_64

; Fall through

; Add magnitudes

add_64:

	mov [DST + sign], eax		; Write sign of result

	cmp N, 1
	ja .large

; Handle the N <= 1 case separately

	xor eax, eax
	mov TMP(d), 2
	
	add OP2, DIGIT(OP1, 0)		; We assume that if N == 0, then DIGIT(OP1, 0) == 0 
	mov DIGIT(DST, 0), OP2

	setnz N(b)			; Store provisional digit count

	setc al				; Store carry into rax
	cmovc N(d), TMP(d)		; Set digit count to 2 if there was a carry

	mov DIGIT(DST, 1), rax		
	mov [DST + size], N(d)

	ret
	
.large:

	cmp DST_CAP(d), N(d)
	jb overflow

	ADD_SUB_64 VARIANT_ADD

	cmp DST_CAP(d), N(d)
	jbe .done

	mov DIGIT(DST, N), rax	; write carry
	add N(d), eax

	xor eax, eax

.done:

	test eax, eax
	setnz [DST + overflow_flag]

	mov [DST + size], N(d)

	ret

; Subtract magnitudes

align 16
subtract_64:

	cmp N, 1
	ja .large

	mov OP1, DIGIT(OP1, 0)
	sub OP1, OP2

	setnz N(b)

	jae .ok

	neg OP1
	xor al, 1

.ok:
	mov DIGIT(DST, 0), OP1
	mov [DST + sign], eax	
	mov [DST + size], N(d)

	ret

.large:

	mov [DST + sign], eax

	cmp DST_CAP(d), N(d)
	jb overflow

	ADD_SUB_64 VARIANT_SUB

	COUNT_DIGITS

	mov [DST + size], N(d)

	ret

; ================================================================================
; void	fpz_add(dstptr dst, srcptr op1, srcptr op2, digcnt dst_cap);
; void	fpz_sub(dstptr dst, srcptr op1, srcptr op2, digcnt dst_cap);
;
; Add/subtract the integers (pointed to by) op1 and op2 and store the result in *dst, 
; which is assumed to have storage for dst_cap digits. Set overflow flag if this 
; is insufficient. If any of non-built in integer operands have the overflow flag 
; set, print error message and abort the entire process (by calling 
; fpz_use_after_overflow_handler).
; ================================================================================

global fpz_add
global fpz_sub

align 16
fpz_add:

	xor HDR2(d), HDR2(d)
	jmp add_sub_common

; Subtraction is equivalent to addition with the second operand negated. 

fpz_sub:

	mov HDR2(d), 1

add_sub_common:

; Load headers, negating the sign of the second operand if subtracting 

	mov HDR1, HEADER(OP1)
	xor HDR2, HEADER(OP2)

	CHECK_OVERFLOW_2 TMP, HDR1, HDR2, use_after_overflow

	mov eax, HDR1(d)
	mov TMP(d), HDR2(d)

	xor TMP(d), eax		; compare signs, store result in TMP (1 if different, 0 if equal)

; get digit counts into N1 and N2 (aka HDR1 and HDR2)

	shr HDR1, 32	
	shr HDR2, 32

; Swap operands and flip sign if OP1 has fewer digits

	cmp N1(d), N2(d)
	jae .ok

; Exchange N1 and N2 (as well as OP1 and OP2) using INDEX as a temporary

	mov INDEX(d), N1(d)		
	mov N1(d), N2(d)
	mov N2(d), INDEX(d)
	
	mov INDEX, OP1
	mov OP1, OP2
	mov OP2, INDEX

; Set the sign of the result to the sign of the second operand (in other words,
; flip the sign of the result if the operand signs are different)

	xor eax, TMP(d)
.ok:

	test TMP(b), 1

	jne subtract

; Fall through 

; Add the magnitudes of the two integers

add:
	cmp N2, 1
	ja .both_large

	mov OP2, DIGIT(OP2, 0)
	jmp add_64

.both_large:

; "Handle" the possibility that we're passed too few digits to hold 
; the shortest possible result.

	cmp DST_CAP(d), N1(d)
	jb overflow

	mov [DST + sign], eax		; Write sign of result 

	xor eax, eax			; Clear carry register

	EXTENDED_LOOP ADD, ADD_CARRY

	cmp DST_CAP(d), N1(d)		; Check if there is room to write the carry
	jbe .done			; If there isn't, just return it regardless of what it is (this serves as an overflow-return)

	mov DIGIT(DST, N1), rax		; Write carry
	add N1(d), eax			; Increment digit count if carry

	xor eax, eax			; Set return value to 0

.done:

	test eax, eax
	setnz [DST + overflow_flag]

	mov [DST + size], N1(d)

	ret

align 16
subtract:

	cmp N2, 1
	ja .both_large

	mov OP2, DIGIT(OP2, 0)
	jmp subtract_64	

.both_large:

; "Handle" the possibility that we're passed too few digits to hold 
; the shortest possible result.

	cmp DST_CAP(d), N1(d)
	jb overflow

	cmp N1, N2
	je .compare

	mov [DST + sign], eax
	xor eax, eax

	EXTENDED_LOOP SUB, SUB_BORROW

	COUNT_DIGITS
	
	mov [DST + size], N1(d)

	ret

.compare:

; The case where N <= 1 is eliminated, so no need to worry about that
; here
	
	lea INDEX, [N - 1]
	
.cmp_loop:

; Digit-by-digit comparison 

	mov TMP, DIGIT(OP1, INDEX)
	cmp TMP, DIGIT(OP2, INDEX)

	ja .ok
	jb .swap

	sub INDEX, 1
	jae .cmp_loop

.zero:

	mov qword DIGIT(DST, 0), 0
	mov qword HEADER(DST), 0

	xor eax, eax 

	ret
	
.swap:

	mov TMP, OP1
	mov OP1, OP2
	mov OP2, TMP
	
	xor eax, 1

.ok:

	mov [DST + sign], eax

	xor eax, eax

	LOOP SUB

	COUNT_DIGITS

	mov [DST + size], N(d)
	
	ret


; Used by every all functions (in this file) which can overflow.

overflow:

	mov qword HEADER(DST), 0x00000100	; set size to 0, sign to positive and overflow_flag to 1
	ret

use_after_overflow:

	USE_AFTER_OVERFLOW_HANDLER

	ret

;================================================================================
; digit fpz_uadd_self(digit* dst, const digit *op, digcnt n);
; digit fpz_usub_self(digit* dst, const digit *op, digcnt n);
;
; Perform unsigned addition/subtraction of *op from *dst, and store the result in
; *dst again. Both pointer parameters are assumed to point to digits, not to the 
; header of an integer. Both *dst and *op are assumed to contain n digits. These 
; are for internal use. 
;================================================================================

; Iteration macros

%macro ADD_SELF 2
%if %2 & RESTORE
	neg al
%endif
	mov TMP, [OP + (%1)*8]
%ifidn %1, 0
	add [DST + (%1)*8], TMP
%else
	adc [DST + (%1)*8], TMP
%endif
%if %2 & SAVE
	setc al
%endif
%endmacro

%macro SUB_SELF 2
%if %2 & RESTORE
	neg al
%endif
	mov TMP, [OP + (%1)*8]
%ifidn %1, 0
	sub [DST + (%1)*8], TMP
%else
	sbb [DST + (%1)*8], TMP
%endif
%if %2 & SAVE
	setc al
%endif
%endmacro

; Function parameters

%define DST	rdi
%define OP	rsi
%define N	rdx

; Variables

%define INDEX	rcx
%define TMP	r8

global fpz_uadd_self
fpz_uadd_self:

	xor eax, eax

	LOOP ADD_SELF

	ret
	
global fpz_usub_self
fpz_usub_self:

	xor eax, eax

	LOOP SUB_SELF

	ret

