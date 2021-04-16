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

;==============================================================================
; long double fpz_get_ld(srcptr op);
;
; Return the long double (80-bit extended precision) value nearest to the value 
; of the integer op (this probably falls into the category of gratuitous assembly,
; so it may be rewritten in C++ at some point).
;==============================================================================

%define OP		rdi

%define EXPONENT	rsi
%define SIGNIFICAND	rax

%define HDR		r8
%define N		r8

%define SIGN		r10

%define EXP_BIAS	16383

global fpz_get_ld
fpz_get_ld:

	mov HDR, HEADER(OP)
	CHECK_OVERFLOW HDR, use_after_overflow
	mov SIGN(d), HDR(d)

	shr HDR, 32

	test N, N
	jz .zero

; To compute the significand, we need the 64 most significant bits. To do this,
; we first load the two most significant digits (if there is only one digit,
; make the other one zero) 

	xor edx, edx			
	mov SIGNIFICAND, DIGIT(OP, N - 1)

	cmp N, 1
	je .l1
	
	mov rdx, DIGIT(OP, N - 2)

.l1:

; Get number of leading zeroes, which is equal to 63 - index of most significant 1-bit

	bsr rcx, SIGNIFICAND
	mov EXPONENT, rcx
	neg cl
	add cl, 63

; Do a double-precision shift to get the most significant 1-bit to bit index 63
; (the extended precision format does not treat the most significant bit of the 
; significand as implicit)

	shld SIGNIFICAND, rdx, cl

; For the purposes of rounding, we need to keep track of remaining bits too, so 
; shift the second digit similarily

	shl rdx, cl

; Compute the exponent value, i.e. bit index (0-based) of most significant bit + EXP_BIAS

	lea r11d, [N(d) - 1]
	shl r11d, 6

	add EXPONENT, r11
	add EXPONENT, EXP_BIAS

; Decide whether or not to round (the magnitude) up (we use round to even as the tie-breaking rule)

	btr rdx, 63	; test the first bit after the significand and clear it afterwards 
	jnc .done	; if that bit was 0, there is no chance we're rounding up

	test rdx, rdx	
	jnz .round_up	; if there are any additional nonzero bits, we round up


; Loop over the remaining digits, checking for any nonzero bits (we use N as an index,
; as it won't be needed after this)

	sub N, 3	; point N to the index of the next digit 
	jb .tie		

.loop:

	cmp qword DIGIT(OP, N), 0
	jne .round_up

	sub N, 1
	jae .loop

.tie:

	test SIGNIFICAND, 1	; test the least significant bit of the significand
	jz .done		; if zero, the significand is even, so we leave it alone.

.round_up:

	add SIGNIFICAND, 1
	jnc .done

; If there was a carry, the significand is now 2^64. We (formally) shift right by 1 bit,
; giving 2^63, and increment the exponent. 

	mov SIGNIFICAND, 0x8000000000000000
	add EXPONENT, 1

.done:

; The largest possible exponent for a normalized extended precision number
; is 32766. For anything above that, we return infinity.

	cmp EXPONENT, 32767
	jae .overflow

; or in the sign bit

	shl SIGN, 15
	or EXPONENT, SIGN

; store floating point value on the stack and retrieve it into st(0)
; (utilize red zone to avoid adjusting stack pointer)

	mov [rsp - 8], EXPONENT
	mov [rsp - 16], SIGNIFICAND

	fld tword [rsp - 16]

	ret

.zero:

	fldz 

	ret

.overflow:

; Return +/- infinity (depending on the sign of the operand)

	mov EXPONENT, 32767
	mov SIGNIFICAND, 0x8000000000000000

	shl SIGN, 15
	or EXPONENT, SIGN

	mov [rsp - 8], EXPONENT
	mov [rsp - 16], SIGNIFICAND

	fld tword [rsp - 16]

	ret
	
use_after_overflow:

	USE_AFTER_OVERFLOW_HANDLER

	ret
