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

%define DST		rdi
%define NBITS		rcx	; passed in rsi, but moved
%define N		rdx

%define INDEX		r8

%define TMP1		rax
%define TMP0		rsi

; ================================================================================
; void	fpz_ushl_1_self(digit *dst, bitcnt nbits, digcnt n)
;
; Shift the first n digits, starting from dst, left by nbits bits, where nbits < 64. 
; Does not consider the sign, and does not write past n digits of dst. 
; (This is for internal use).
; ================================================================================

global fpz_ushl_1_self
fpz_ushl_1_self:

	mov NBITS, rsi	
	and NBITS, 63	

	sub N(d), 1
	mov TMP0, [DST + N*8]

.loop:
	
	mov TMP1, TMP0 
	mov TMP0, [DST + N*8 - 8]

	shld TMP1, TMP0, NBITS(b)
	mov [DST + N*8], TMP1

	sub N(d), 1
	ja .loop

.last:
	
	shl TMP0, NBITS(b)
	mov [DST], TMP0

	ret
	
; ================================================================================
; void fpz_shl(dstptr dst, srcptr op, bitcnt nbits, bitcnt prec)
;
; Shift the integer pointed to by op left by nbits bits and store the result in dst,
; which is assumed to have space for prec bits. If this is not sufficient, set
; the overflow flag in dst. If the overflow flag is set in op, print an error
; message and abort the process (by calling fpz_use_after_overflow_handler)
; ================================================================================

%define DST		rdi
%define OP		rsi
%define PREC		rdx	; passed in rcx, but moved
%define NBITS		rcx	; passed in rdx, but moved

%define HDR		r8
%define N		r8

%define TMP0		r9
%define TMP1		r10

%define LIMIT_OP	rdx	; overlap with PREC, but used only after we're done with it
%define LIMIT_DST	r8	; overlap with N, used only after we're done using N

global fpz_shl
fpz_shl:

; Ensure that nbits and prec reside in correct registers (shl/shld requires
; that we use cl to store the number of bits)

	mov TMP0, rdx		; save nbits 
	mov PREC, rcx		
	mov NBITS, TMP0		

	mov HDR, HEADER(OP)
	CHECK_OVERFLOW HDR, use_after_overflow

	mov TMP1(d), HDR(d)	
	mov [DST], TMP1(d)	; save sign (no harm even if DST and OP coincide)

	shr HDR, 32
	jz zero

	mov TMP0, [OP + N*8]	; get the most significant nonzero digit
	bsr rax, TMP0

; Compute the number of digits required in the result operand

	sub N, 1
	mov r11, N

	shl r11, 6		; multiply digit count by 64
	add r11, rax		; add in most significant bit index (number of bits in operand)

	add r11, NBITS		; add the shift, to get the number of required bits for the result

; r11 now contains the total number of bits required to hold the result

	cmp PREC, r11
	jb overflow		
	
	shr r11, 6	; get digit index of the most significant digit in res
	add r11, 1

	mov [DST + size], r11(d)	; save digit count 

; We will loop from most to least significant digit of OP and DST, respectively.
; Point OP and DST to their respective most significant digits and set their original
; values as lower limits.

	mov LIMIT_OP, OP
	lea OP, [OP + N*8]	

	mov LIMIT_DST, DST	; Note: this will overwrite N
	lea DST, [DST + r11*8]	

	and NBITS, 63

	add eax, NBITS(d)	
	cmp eax, 64
	jb .pre_loop

	xor TMP1(d), TMP1(d)
	shld TMP1, TMP0, NBITS(b)

	mov [DST], TMP1

	sub DST, 8

.pre_loop:

	cmp OP, LIMIT_OP
	je .last
	jb .zero_loop

.loop:

	mov TMP1, TMP0
	mov TMP0, [OP]

	shld TMP1, TMP0, NBITS(b)

	mov [DST], TMP1

	sub DST, 8
	sub OP, 8
	
	cmp OP, LIMIT_OP
	ja .loop 

.last:

	shl TMP0, NBITS(b)
	mov [DST], TMP0

	sub DST, 8
	cmp DST, LIMIT_DST

	jbe .done

.zero_loop:

	mov qword [DST], 0
	
	sub DST, 8
	cmp DST, LIMIT_DST

	ja .zero_loop

.done:
	xor eax, eax

	ret

overflow:

	mov eax, 1

	mov qword DIGIT(DST, 0), 0
	mov qword HEADER(DST), 0

	ret

zero:
	
	xor eax, eax

	mov qword DIGIT(DST, 0), 0
	mov qword HEADER(DST), 0

	ret
; ================================================================================
; digit fpz_shr(dstptr dst, srcptr op, bitcnt nbits, bitcnt prec)
;
; Shift the integer pointed to by op right by nbits bits and store the result in dst,
; which is assumed to have space for prec bits. If this is not sufficient, set
; the overflow flag in dst. If the overflow flag is set in op, print an error
; message and abort the process (by calling fpz_use_after_overflow_handler)
; ================================================================================

global fpz_shr
fpz_shr:

	mov TMP0, rdx
	mov PREC, rcx
	mov NBITS, TMP0

	mov HDR, HEADER(OP)
	CHECK_OVERFLOW HDR, use_after_overflow
	mov TMP1(d), HDR(d)

	mov [DST], TMP1(d)

	shr HDR, 32
	jz zero

	mov TMP0, [OP + N*8]
	bsr rax, TMP0

; Compute the number of bits required in the result operand

	lea r11, [N - 1]
	
	shl r11, 6
	add r11, rax

	sub r11, NBITS
	jb zero

	cmp r11, PREC
	jae overflow

	shr r11, 6
	add r11, 1

	mov [DST + size], r11(d)

	mov TMP1, NBITS
	and NBITS, 63	

	shr TMP1, 6
	
	lea LIMIT_OP, [OP + N*8]
	lea OP, [OP + TMP1*8 + 16]

	lea LIMIT_DST, [DST + r11*8]
	lea DST, [DST + 8]	

	mov TMP1, [OP - 8]

	cmp DST, LIMIT_DST
	jae .last

.loop:
	mov TMP0, TMP1
	mov TMP1, [OP]

	shrd TMP0, TMP1, cl

	mov [DST], TMP0
	
	add OP, 8
	add DST, 8

	cmp DST, LIMIT_DST

	jb .loop

.last:
	mov TMP0, TMP1
	xor TMP1(d), TMP1(d)

	cmp OP, LIMIT_OP
	ja .label1

	mov TMP1, [OP]
.label1:

	shrd TMP0, TMP1, cl
	test TMP0, TMP0
	jz .done

	mov [DST], TMP0

.done:
	xor eax, eax
	ret

; ================================================================================
; digit fpz_shl_i64(dstptr dst, int64_t op, bitcnt nbits, bitcnt prec)
; ================================================================================

global fpz_shl_i64
fpz_shl_i64:

	mov TMP0, rdx
	mov PREC, rcx
	mov NBITS, TMP0

	SIGN_MAGNITUDE OP, rax, OP
	jz zero

	and rax, 1
	mov [DST], eax		; store sign bit

	bsr rax, OP

	mov r11, rax
	add rax, NBITS

	cmp PREC, rax
	jb overflow 

	shr eax, 6		; compute digit index which the most significant digit goes into
	add eax, 1

	mov [DST + size], eax	; store digit count

	and NBITS, 63	

	xor TMP0(d), TMP0(d)
	
	shld TMP0, OP, cl
	test TMP0, TMP0

	jz .last

	mov [OP + rax*8], TMP0
	sub eax, 1

.last:

	shl OP, cl
	mov [DST + rax*8], OP

	lea LIMIT_DST, [DST + rax*8]
	add DST, 8

.zero_loop:

	mov qword [DST], 0

	add DST, 8
	cmp DST, LIMIT_DST
	jb .zero_loop

	xor eax, eax
	ret

use_after_overflow:

	USE_AFTER_OVERFLOW_HANDLER

	ret
