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

struc int_reciprocal
	m0: 		resq	1	; multiplier, least significant digit
	m1: 		resq	1	; multiplier, most significant digit
	norm_shift:	resq	1	; normalization shift
	d1:		resq 	1	; normalized divisor, most significant digit
	d0:		resq	1	; normalized divisor, least significant digit
endstruc		

%define IR		rdi
%define D2		rsi
%define D1		rdx
%define D0		rax	; passed in rcx, but will be moved

%define DIVISOR 	rsi	; note: same register as D2, just alternate name

%define N0		r8
%define N1		r9

%define M0		r10
%define M1		r11

;==============================================================================
; void fpz_reciprocal_64(int_reciprocal *ir, uint64_t d2, uint64_t d1, uint64_t d0)
;
; Given three digits d2, d1, d0 (where d2 is most significant, etc.), compute 
; the integer reciprocal for dividing by the leading 64 bits of d2:d1:d0. 
; 
;
; This is probably not the fastest method.
;==============================================================================

global fpz_reciprocal_64
fpz_reciprocal_64:

; Compute the normalization shift.

; The shld instruction is very picky about what registers can hold the bit counts, 
; so get the least significant digit out of rcx first. 

	mov D0, rcx

	bsr rcx, D2

	jz .div_by_zero			; trigger a division by zero (will result in a SIGFPE on POSIX-compliant systems)
	neg ecx
	add ecx, 63			; ecx now contains the number of leading zeroes in denominator		
	
	shld D2, D1, cl
	shld D1, D0, cl

	mov [IR + norm_shift], rcx	; store the normalization shift
	mov [IR + d1], D2		; store the most significant digit of the normalized divisor
	mov [IR + d0], D1		; store the least significant digit of the normalized divisor.
	
; Some special handling is required for a normalized divisor which is exactly 2^63.
; Make a copy, shift left by one and see if the result is zero.
	
	mov rax, DIVISOR
	shl rax, 1
	jz .maximum

; The entire multiplier can be computed using two 128 by 64 bit divisions, but that is quite
; costly, so we will try to be slightly more clever (since we're dividing by d, the leading
; bits of the multiplier will help us in computing further bits).
; 
; We begin by computing the leading 32 bits of the multiplier. Those can be expressed as floor(2^(32)*(2^(64) - d)/d),
; i.e. a 96-bit integer divided by a 64-bit integer (with the last 32 bits equal to 0). This is dealt with
; using a standard method, (apparently) found in Knuth's TACP (Volume 2, section 4.3.1); applied to our problem, and using
; the instructions available to us, it involves dividing the leading 64 bits of the numerator (i.e. 2^(64) - d) by
; the leading 32 bits of the denominator, and then correcting the result in at most two steps (the quotient thus obtained
; is always >= the true quotient and the difference is at most 2. 

	mov N0, DIVISOR
	neg N0			; High order digit of numerator is 2^(64) - DIVISOR
	xor N1, N1		; Low order digit of numerator is 0

; shift dividend by 32 bits, making a 96 bit integer with 32 trailing zeroes

	shld N1, N0, 32
	shl N0, 32		

	mov rcx, DIVISOR
	shr rcx, 32		; place leading 32 bits of DIVISOR into rcx

	mov edx, N1(d)
	mov rax, N0		; place the low 32 bits of N1 into edx
	shr rax, 32		; put the high 32 bits of N0 into eax

	div ecx			; divide edx:eax by ecx, quotient is now in eax, remainder in edx

	mov M1(d), eax		; save the quotient

	mul DIVISOR		; multiply the quotient by the (entire) divisor

	sub N0, rax
	sbb N1, rdx		; Compute N1:N0 - M1*DIVISOR

	jae .32_bits_done	; if result is >= 0, the 32 quotient bits are correct

	sub M1, 1

	add N0, DIVISOR
	adc N1, 0		; add DIVISOR to N1:N0

	jc .32_bits_done	; if we're back above (or equal) to zero, we're done

	sub M1, 1		; otherwise, make one more adjustment

	add N0, DIVISOR
	adc N1, 0

.32_bits_done:

; We now know the leading 32 bits of the multiplier, and we now proceed to determine the next 32 bits. 
; This requires dividing the remainder 2^(32)*N0 by DIVISOR. Since we have the initial 32 bits of the 
; multiplier, we can use this to approximate the quotient to within one unit again, so this will be 
; a multiplication and a few adds/subs and branches again.

	mov rax, N0
	mul M1		; compute M1*N0 (96 bits)
 
	xor N1, N1

; Multiply N0 by 2^32

	shld N1, N0, 32	
	shl N0, 32
	
	add rax, N0
	adc rdx, N1

	mov ecx, edx	; save quotient (can't put it in M1 yet, because any adjustments would kill the upper 32 bits)

	mov rax, DIVISOR
	mul rcx		; multiply DIVISOR by the low 32 bits of M1

	sub N0, rax
	sbb N1, rdx

	add ecx, 1	; this may overflow, but if it does, the readjustment will happen (the correct value fits in 32 bits)

	sub N0, DIVISOR	; subtract DIVISOR once more
	sbb N1, 0

	jae .64_bits_done	; no adjustment needed if result >= 0	

; oops, that was too much, decrement the quotient and add the DIVISOR back into remainder 

	sub ecx, 1

	add N0, DIVISOR
	adc N1, 0

.64_bits_done:

	mov N1, N0
	xor N0, N0	; shift the remainder by 64 bits

	shl M1, 32	; shift M1 left by 32 bits to make room for the low 32 bits	
	or M1, rcx	; get the remaining bits into M1, completing the leading digit of the multiplier

	mov [IR + m1], M1

; Compute the low 64 bits of the multiplier. Do this by using the high 64 bits as an approximation to the 
; actual multiplier, and use this to divide. The result may be off by one (one unit too small), so an adjustment 
; may be needed. 

	mov rax, N1	; get the remainder into rax
	mul M1		; compute M1*N1

	add rdx, N1	; rdx now contains a candidate for the last digit of the multiplier, but it may be one unit too small.

	mov M0, rdx

	mov rax, DIVISOR
	mul M0

	sub N0, rax
	sbb N1, rdx

	add M0, 1

	sub N0, DIVISOR
	sbb N1, 0

	sbb M0, 0	; if M0 is too large at this point, this will subtract 1

; No need to adjust remainder at this point

	mov [IR + m0], M0

	ret

.maximum:

	mov qword [IR + m0], -1
	mov qword [IR + m1], -1

	ret
	
.div_by_zero:

	mov qword [IR + m0], -1
	mov qword [IR + m1], -1
	mov qword [IR + norm_shift], -1
	mov qword [IR + d1], 0
	mov qword [IR + d0], 0

; Perform a division by zero and let the OS do what it must do

	xor eax, eax
	div al		

	ret

%define N2		rdi
%define N1		rsi	
%define N0		rdx
%define IR		rcx

%define REM		r8

%define M		r9
%define QUOT		r9

%define TMP0		r10
%define TMP1		r11

;==============================================================================
; uint64_t fpz_div128_64(uint64_t n2, uint64_t n1, uint64_t n0, int_reciprocal *ir)
;
; Shift n2:n1:n0 by the normalization shift of the divisor, cut off the least
; significant digit and divide by the integer corresponding to the reciprocal given
; by *ir (this is done using only multiplications). Return the quotient digit.
;==============================================================================

global fpz_div128_64
fpz_div128_64:

	mov rax, IR			; vacate rcx
	mov rcx, [IR + norm_shift]	; load the normalization shift
	
	shld N2, N1, cl
	shld N1, N0, cl

	mov IR, rax

	test N2, N2
	jz .single_digit_div		; things get easier if N2 is zero

	mov M, [IR + m0]

	MUL_64x64 rax, TMP0, N1, M, N1		; discard low qword
	MUL_64x64 TMP0, TMP1, N2, M, N2, TMP0

	mov M, [IR + m1]
	
	MUL_64x64 rax, TMP0, N1, M, TMP0		; discard low qword 
	MUL_64x64 TMP0, TMP1, N2, M, TMP0, TMP1

; We want to return 2^64 - 1 in case the result overflows. First clear out rax ...

	xor eax, eax

	add TMP0, N1
	adc TMP1, N2

; ... and then subtract the carry, setting all bits to 1 if the result overflowed (0 otherwise)
	sbb rax, 0

; does nothing if all bits of rax are set (overflow), copies the bits of TMP1 otherwise (if rax = 0)

	or rax, TMP1

	ret

.single_digit_div:

	xor eax, eax
	mov edx, 1		; the quotient is either zero or one

	cmp N1, [IR + d1]	; compare N1 with the denominator
	cmovae eax, edx

	ret	


;==============================================================================
; uint64_t fpz_divrem128_64(uint64_t n2, uint64_t n1, uint64_t n0, int_reciprocal *ir, uint64_t *rem)
;
; Same as previous function, except that this also computes the remainder and stores it
; in *rem.
;==============================================================================

global fpz_divrem128_64
fpz_divrem128_64:

	mov rax, IR			; vacate rcx

	mov rcx, [IR + norm_shift]	; load the normalization shift
	
	shld N2, N1, cl
	shld N1, N0, cl

	mov IR, rax

	test N2, N2
	jz .single_digit_div		; things get easier if N2 is zero

	mov M, [IR + m0]

	MUL_64x64 rax, TMP0, N1, M, N1		; discard low qword
	MUL_64x64 TMP0, TMP1, N2, M, N2, TMP0

	mov M, [IR + m1]
	
	MUL_64x64 rax, TMP0, N1, M, TMP0		; discard low qword 
	MUL_64x64 TMP0, TMP1, N2, M, TMP0, TMP1

; We want to return 2^64 - 1 in case the result overflows. First clear out rax ...

	xor eax, eax

	add TMP0, N1
	adc TMP1, N2

; ... and then subtract the carry, setting all bits to 1 if the result overflowed (0 otherwise)
	sbb rax, 0

; Does nothing if all bits of rax are set (overflow), copies the bits of TMP1 otherwise (if rax = 0)

	or rax, TMP1

	mov QUOT, rax		; save the quotient digit

	mul qword [IR + d1]	; multiply quotient by leading digit of denominator

	sub N1, rax	; subtract the low qword of product from N1

	sbb N2, rdx	; subtract the high qword of product from N2 (with borrow)
	mov rax, -1

	test N2, N2
	cmovnz N1, rax	; if N2 is nonzero, set N1 = UINT64_MAX
	
	mov [REM], N1

	mov rax, QUOT

	ret

.single_digit_div:
	
	xor eax, eax
	xor edx, edx

	mov TMP0, [IR + d1]
	cmp N1, TMP0

	setae al
	cmovae rdx, TMP0

	sub N1, rdx

	mov [REM], N1

	ret

;================================================================================
; uint64_t fpz_div10_128(uint64_t n2, uint64_t n1, uint64_t n0)
;
; Special case of division by 10.
;================================================================================

global fpz_div10_128
fpz_div10_128:

	shld N2, N1, 60
	shld N1, rdx, 60

	mov M, 0x9999999999999999

	MUL_64x64 rax, TMP0, N1, M, N1		; discard low qword
	MUL_64x64 TMP0, TMP1, N2, M, N2, TMP0

	MUL_64x64 rax, TMP0, N1, M, TMP0		; discard low qword 
	MUL_64x64 TMP0, TMP1, N2, M, TMP0, TMP1

; We want to return 2^64 - 1 in case the result overflows. First clear out rax ...

	xor eax, eax

	add TMP0, N1
	adc TMP1, N2

; ... and then subtract the carry, setting all bits to 1 if the result overflowed (0 otherwise)
	sbb rax, 0

; does nothing if all bits of rax are set (overflow), copies the bits of TMP1 otherwise (if rax = 0)

	or rax, TMP1	

	ret


;================================================================================
; uint64_t fpz_adjust_quotient(const uint64_t *rd, uint64_t n1, uint64_t n0, uint64_t quot, uint64_t rem)
;
;
;================================================================================

%define IR	rdi
%define N1	rsi
%define N0	rdx
%define QUOT	rcx
%define REM	r8

%define D1	r9
%define D0	r10

global fpz_adjust_quotient
fpz_adjust_quotient:

	mov rax, rcx		; vacate rcx temporarily
	mov rcx, [IR + 16]	; load the normalization shift

	shld N1, N0, cl		; rdx is now free (only N1 is used)
	mov rcx, rax		; get the quotient digit back to where it was

	mul qword [IR + 32]		; multiply quotient digit by last digit of (shifted) denominator

	sub N1, rax
	sbb REM, rdx

	jae .end

	mov D1, [IR + 24]
	mov D0, [IR + 32]

	sub QUOT, 1

	add N1, D0
	adc REM, D1

	jc .end

	sub QUOT, 1

.end:

	mov rax, QUOT
	ret

;================================================================================
; uint64_t fpz_div192_128_direct(uint64_t d2, uint64_t d1, uint64_t d0, uint64_t n3, uint64_t n2, uint64_t n1, uint64_t n0)
;
; Shift the three digit integer d2:d1:d0 into a 128-bit integer d, and shift the four digit integer n3:n2:n1:n0
; into the 192-bit n, and then compute and return floor(n/d).
;================================================================================

%define D2	rdi
%define D1	rsi
%define D0	rdx
%define N3	rcx
%define N2	r8
%define N1	r9
%define N0	[rsp + 8]

%define QUOT	r10
%define REM	r11

global div192_128_direct
div192_128_direct:
	
	mov rax, N3		; get N3 temporarily out of the way, to free up rcx for other stuff

	bsr rcx, D2
	neg ecx
	add ecx, 63

	shld D2, D1, cl
	shld D1, D0, cl

	shld rax, N2, cl
	shld N2, N1, cl
	
	mov r10, N0		; load N0 from the stack
	shld N1, r10, cl	

	mov N3, rax
	
	mov rdx, N3
	mov rax, N2

	div D2			; quot is in rax, remainder in rdx

	mov REM, rdx
	mov QUOT, rax		; get the results out of the way

	mul D1			; multiply quotient by second digit of denominator

	sub N1, rax
	sbb REM, rdx

	jae .end

	sub QUOT, 1

	add N1, D1
	adc REM, D2

	jc .end

	sub QUOT, 1

.end:
	mov rax, QUOT
	ret
	
