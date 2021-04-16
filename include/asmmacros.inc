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

struc fpz_header
	sign:		resb	1
	overflow_flag:	resb	1
	pad:		resb	2
	size:		resd	1
endstruc

; pick out the overflow flag from the low 32 bits of the header
%define OVERFLOW_MASK	0x0000FF00	

%define DIGIT(op, index)	[op + (index)*8 + 8]
%define HEADER(op)		[op]

%macro CHECK_OVERFLOW 2
%push 
%define %$HDR %1
%define %$LABEL %2
	test %$HDR(d), OVERFLOW_MASK
	jnz %$LABEL
%pop
%endmacro

%macro CHECK_OVERFLOW_2 4
%push 
%define %$TMP 	%1
%define %$HDR1	%2
%define %$HDR2	%3
%define %$LABEL %4
	mov %$TMP(d), %$HDR1(d)
	or %$TMP(d), %$HDR2(d)
	test %$TMP(d), OVERFLOW_MASK
	jnz %$LABEL
%pop
%endmacro
	
%macro USE_AFTER_OVERFLOW_HANDLER 0
extern fpz_use_after_overflow_handler

	sub rsp, 8	; make sure stack is aligned properly before function call 
			; (assumes that rsp is unmodified after previous function entry point)

	call fpz_use_after_overflow_handler

; just in case the function returns

	add rsp, 8
%endmacro


; For registers r8 through r15, the partial registers, i.e. lowest byte,
; word and dword, can be accessed by appending 'b', 'w' and 'd', respectively, 
; to the end of the register name. This defines the corresponding names for 
; the remaining general purpose registers, making it convenient to access 
; partial registers for registers that are hidden behind macro names. 

%define raxb	al
%define raxw	ax
%define raxd	eax

%define rbxb	bl
%define rbxw	bx
%define rbxd	ebx

%define rcxb	cl
%define rcxw	cx
%define rcxd	ecx

%define rdxb	dl
%define rdxw	dx
%define rdxd	edx

%define rsib	sil
%define rsiw	si
%define rsid	esi

%define rdib	dil
%define rdiw	di
%define rdid	edi

%define rbpw	bp
%define rbpd	ebp

%macro DEFINE_PARTIAL_REG_MACROS 1-*
%rep %0
%define %1(x)	%1 %+ x
%rotate 1
%endrep
%endmacro

DEFINE_PARTIAL_REG_MACROS rax, rcx, rdx, rdi, rsi, r8, r9, r10, r11, r12, r13, r14, r15, rbx, rbp

; WRITE_TO_REDZONE reg0, ...
; READ_FROM_REDZONE reg0, ...
;
; Save/restore the registers given as arguments to/from the red zone.

%define REDZONE_SIZE	128

%define NOWHERE

%macro WRITE_TO_REDZONE	1-*
%push
%assign %$OFFSET	8
%rep %0
%if %$OFFSET >= REDZONE_SIZE
%error "Red zone exhausted: stored values are likely to be corrupted"
%endif
mov [rsp - %$OFFSET], %1
%rotate 1
%assign %$OFFSET	(%$OFFSET + 8)
%endrep
%pop
%endmacro

%macro READ_FROM_REDZONE 1-*
%push
%assign %$OFFSET 	8
%rep %0
%if %$OFFSET >= REDZONE_SIZE
%error "Red zone exhausted: stored values are likely to be corrupted"
%endif
%ifnempty %1
mov %1, [rsp - %$OFFSET]
%endif
%rotate 1
%assign %$OFFSET (%$OFFSET + 8)
%endrep
%pop
%endmacro 

; SIGN_MAGNITUDE INTEGER, SIGN, MAGNITUDE
;
; Split a 64-bit two's complement signed integer (first parameter) into sign 
; (second parameter) and magnitude (third parameter). Uses no branches. 
; The sign will have the value -1 if negative or 0 if positive. 

%macro SIGN_MAGNITUDE 3
%push
%define %$INTEGER		%1
%define %$SIGN		%2
%define %$MAGNITUDE	%3

%ifidn %$SIGN, %$MAGNITUDE
%error "SIGN_MAGNITUDE: Storage for sign (second argument) and storage\
 for magnitude (third argument) can't coincide"
%endif

%ifnidn %$INTEGER, %$MAGNITUDE
	mov %$MAGNITUDE, %$INTEGER
%endif
	mov %$SIGN, %$MAGNITUDE
	sar %$SIGN, 63		
	
	xor %$MAGNITUDE, %$SIGN
	sub %$MAGNITUDE, %$SIGN
%pop
%endmacro

; ERROR_IF_MATCH SYM, ERRMSG, ...
;
; Macro to check if SYM matches anything on a given list (the arguments passed
; after ERRMSG) and signal an error if it does (with the message given as the 
; second argument). Doing this with just ifidn necessitates a lot of repetition 
; (AFAICT, there is no way to say 'if <something> or <something_else> then ... ') 
; without a preprocessor loop as below. 

%macro ERROR_IF_MATCH 2-*
%push
%define %$SYM 	%1
%define %$ERRMSG	%2
%rotate 1
%rep %0 - 2
%rotate 1
%ifidn %$SYM, %1
%error %$ERRMSG
%endif
%endrep
%pop
%endmacro

; MUL_64x64 STORAGE1, STORAGE2, FACTOR1, FACTOR2, <TERM1>, <TERM2>
;
; Writes code to compute FACTOR1*FACTOR2 + TERM1 + TERM2 and assign the result 
; to STORAGE2:STORAGE1. TERM1 and TERM2 are optional (no code is generated to
; handle them, unless they're actually given). 
;
; Caveats:
; - TERM1 and TERM2 must be distinct from rax and rdx, since the latter 
;   registers are clobbered by the mul instruction (we check for this).
; - The macro does not try to decide if a parameter is a register or a 
;   memory location. Certain parameter combinations may generate
;   instructions with two memory operands, which will simply fail to assemble.

%macro MUL_64x64 4-6
%push
%define %$STORAGE1 	%1
%define %$STORAGE2 	%2
%define %$FACTOR1 	%3
%define %$FACTOR2 	%4
%define %$TERM1		%5
%define %$TERM2		%6

; Check that we're not doing anything too stupid.

ERROR_IF_MATCH %$TERM1, "MUL_64x64: First term can't be either rax or rdx", rax, rdx
ERROR_IF_MATCH %$TERM2, "MUL_64x64: Second term can't be either rax or rdx", rax, rdx

%ifidn %$STORAGE1, rdx
%error "MUL_64x64: rdx can't be passed as first storage location."
%endif

%ifidn %$STORAGE1, %$STORAGE2
%error "MUL_64x64: Identical storage locations specified"
%endif

; Multiply FACTOR1 and FACTOR2 (simplest when either of them is rax)

%ifidn %$FACTOR1, rax
	mul qword %$FACTOR2
%elifidn %$FACTOR2, rax
	mul qword %$FACTOR1
%else
	mov rax, %$FACTOR1
	mul qword %$FACTOR2
%endif

; Add in any optional terms.

%if %0 > 4		
%ifidn %$STORAGE1, %$TERM1
; Use TERM1 as destination to avoid a copy at the end
	add %$TERM1, rax
%else
	add rax, %$TERM1
%endif
	adc rdx, 0
%endif

%if %0 > 5
%ifidn %$STORAGE1, %$TERM1
; Use TERM1 as destination to avoid a copy at the end
	add %$TERM1, %$TERM2
%else
	add rax, %$TERM2
%endif
	adc rdx, 0
%endif

; If the low qword of the result is in rax (which it is, unless STORAGE1 
; coincides with TERM1) and it does not belong there, copy it into STORAGE1.

%ifnidn %$STORAGE1, rax
%ifnidn %$STORAGE1, %$TERM1
	mov %$STORAGE1, rax
%endif
%endif

; Copy the high qword of the result (rdx) into STORAGE2, unless it is already 
; there.

%ifnidn %$STORAGE2, rdx
	mov %$STORAGE2, rdx
%endif
%pop
%endmacro


; LOOP ITERATION
; 
; Simple looping macro, which executes a collection of instructions, supplied
; in the form of an iteration macro, N times (where N is assumed to refer to
; a register or memory location). Presently, it is unrolled by a factor of
; 2, but this could be changed. The present iteration count (ranging from 0 to 
; N - 1), is maintained in the register/memory location that the macro INDEX 
; refers to. 
;
; The iteration macro is supposed to expect two parameters, where the first is 
; the present iteration index and the other is an integer indicating whether
; the macro should save and/or restore some particular state. This is indicated
; using the flags defined immediately below, and is (presently) used exclusively 
; to deal with the saving/restoring of the carry flag in the context of addition 
; or subtraction.
;
; Assumes that N > 0 (so this needs a test beforehand to determine if N == 0).
; This is inconsistent with the EXTENDED_LOOP macro. 

%define DO_NOTHING	0
%define RESTORE		1
%define SAVE		2

%macro LOOP 1
%push
%define %$ITERATION 	%1

	xor INDEX(d), INDEX(d)

	test N(d), 1
	jz %%loop

	%$ITERATION 0, SAVE

	mov INDEX(d), 1

	cmp N(d), 1
	je %%end

%%loop:

	%$ITERATION INDEX, RESTORE
	%$ITERATION (INDEX + 1), SAVE

	add INDEX(d), 2
	cmp INDEX(d), N(d)

	jb %%loop
%%end:
%pop
%endmacro

; EXTENDED_LOOP ITERATION1, ITERATION2
;
; Same as above, but this repeats the instructions given in ITERATION1 N2 times
; and ITERATION2 the remaining N1 - N2 times. It assumes that N1 >= N2 and that
; N1 >= 0 and N2 >= 0.


%macro EXTENDED_LOOP 2
%push
%define %$ITERATION1	%1
%define %$ITERATION2	%2

	xor INDEX(d), INDEX(d)

;	test N2, N2
;	jz %%end_loop1

	test N2, 1
	jz %%loop1

	%$ITERATION1 0, SAVE

	mov INDEX, 1
	
;	cmp N2, 1
;	je %%end_loop1

%%loop1:

	%$ITERATION1 INDEX, RESTORE
	%$ITERATION1 (INDEX + 1), SAVE

	add INDEX, 2
	cmp INDEX, N2

	jb %%loop1

%%end_loop1:

	cmp INDEX, N1
	je %%end

%%loop2:

	%$ITERATION2 INDEX, RESTORE | SAVE

	add INDEX, 1
	cmp INDEX, N1

	jb %%loop2

%%end:
%pop
%endmacro

; Common iteration macros.

%macro COPY 2
	mov TMP, [OP + %1*8 + 8]
	mov [DST + %1*8 + 8], TMP
%endmacro

