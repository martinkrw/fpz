// Copyright (c) 2021 Martin Wanvik <martin.kr.wanvik@gmail.com>
//
// Permission to use, copy, modify, and distribute this software for any
// purpose with or without fee is hereby granted, provided that the above
// copyright notice and this permission notice appear in all copies.
//
// THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
// WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
// ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
// WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
// ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
// OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

#include <numbers>
#include "fpz_core.h"

namespace fpz::core {

// Estimate the number of decimal digits in an integer

static FPZ_FORCE_INLINE size_t
base10_size_bound(srcptr op)
{
	bitcnt nbits = fpz_bit_size(op);
	if (nbits == bitcnt{0}) return 1;

	constexpr long double log10_2 = std::numbers::ln2_v<long double> /
	    std::numbers::ln10_v<long double>;
	
	return static_cast<size_t>(log10_2 * nbits.nb) + 1;
}

static FPZ_FORCE_INLINE size_t
base16_size_bound(srcptr op)
{
	bitcnt nbits = fpz_bit_size(op);
	if (nbits == bitcnt{0}) return 1;

	return (nbits.nb + 3) / 4;
}

// Extract the index'th hexadecimal digit (counted from least to most significant) from
// a 64-bit digit value, and return it as a character.

static FPZ_FORCE_INLINE char 
extract_hex_digit(digit value, int index)
{
	int rshift = index*4;
	digit hex_digit = (value >> rshift) & digit{0xF};
	return (hex_digit > 9 ? 'A' + hex_digit - 10 : '0' + hex_digit);
}

// Low-level, internal, unsafe function for writing an integer into a string (in hexadecimal). 
// If negative, a leading minus sign is added. The string is assumed to contain enough 
// characters to hold the value of the integer (including its sign). 
// (This is made into a separate function to make it easier to support other I/O APIs)

static void
fpz_sprint_hex(char *dst, srcptr src)
{
	header hdr = get_header(src);

	if (hdr.size == digcnt{0}) {
		*(dst++) = '0';
		*dst = '\0';
		return;
	}

	if (hdr.sign == NEGATIVE) *(dst++) = '-';

	digit cur_digit = digits(src)[hdr.size.nd - 1];

	// Skip over leading zero digits

	int i = 15;
	for (; extract_hex_digit(cur_digit, i) == '0'; --i);

	for (; i >= 0; --i) *(dst++) = extract_hex_digit(cur_digit, i);
	
	for (digcnt j = hdr.size - digcnt{2}; j < hdr.size; --j) {
		cur_digit = digits(src)[j.nd];
		for (int i = 15; i >= 0; --i) *(dst++) = extract_hex_digit(cur_digit, i);
	}

	*dst = '\0';
}

// Write the value (in hexadecimal) of an integer to a stream.

extern "C" int
fpz_fprint_hex(FILE *stream, srcptr arg) noexcept
{
	// Allocate space for digits, sign and null terminator 
	size_t sz = base16_size_bound(arg) + 2;
	char *str = static_cast<char*>(alloca(sz));

	fpz_sprint_hex(str, arg);
	return std::fputs(str, stream);
}

// Internal, unsafe function for writing an integer into a string (in decimal). It 
// is inconsistent with fpz_sprint_hex, due to the fact that we extract digits from 
// least to most significant (and print in the opposite order). Thus, we assume that 
// dst points to where the null-terminator goes and write digits backwards from there, 
// ending with a minus if the sign is negative. Return a pointer to the last character 
// written this way, i.e. the beginning of the string (which may or may not be the 
// beginning of the allocated string, so one may need a memmove if that is required).
// (This is made into a separate function to make it easier to support other I/O APIs)

static char *
fpz_sprint_dec(char *dst, srcptr op)
{
	// Nul-terminate 
	*dst = '\0';

	header hdr = get_header(op);

	if (hdr.size == digcnt{0}) {
		*(--dst) = '0';
		return dst;
	}

	// Make a (stack-allocated) copy of (the absolute value of) op 
	// that we can modify

	size_t sz = int_alloc_size(hdr.size);
	dstptr copy = static_cast<dstptr>(alloca(sz));
	fpz_set(copy, op, hdr.size);
	set_sign(copy, POSITIVE);

	while (get_size(copy) != digcnt{0}) {
		int dec_digit = fpz_divrem10(copy, copy, hdr.size);
		*(--dst) = '0' + dec_digit;
	}

	if (hdr.sign == NEGATIVE) *(--dst) = '-';

	return dst;
}

extern "C" int
fpz_fprint_dec(FILE *stream, srcptr arg) noexcept
{
	// Allocate space for digits, sign and nul-terminator
	size_t sz = base10_size_bound(arg) + 2; 
	char *str = static_cast<char*>(alloca(sz));

	return std::fputs(fpz_sprint_dec(str + sz - 1, arg), stream);
}

} // namespace fpz::core
