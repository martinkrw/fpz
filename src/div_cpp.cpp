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

#include <cmath>
#include <cstdlib>	/* For alloca */
#include <cinttypes>
#include "fpz_core.h"

namespace fpz::core {

int64_t
fpz_i64_divrem(int64_t *quot, int64_t num, srcptr denom) noexcept
{
	return i64_divrem(quot, num, denom);
}

struct int_reciprocal {
	digit 		multiplier[2];
	uint8_t 	norm_shift;	// normalization shift
	digit 		norm_divisor[2];// normalized divisor
};

// Make sure no unexpected padding occurs
static_assert(sizeof(int_reciprocal) == 5*sizeof(digit), "unexpected layout of struct reciprocal_data");

extern "C" {

digit	fpz_uadd_self(digit *dst, const digit *op, digcnt n) noexcept;
digit 	fpz_usub_self(digit *dst, const digit *op, digcnt n) noexcept;
void	fpz_ushl_1_self(digit *dst, bitcnt bitcnt, digcnt n) noexcept;

digit	fpz_usubmul_64_self(digit *dst, const digit *op1, digit op2, digcnt n) noexcept;
digit	fpz_usubmul_lowcap_64_self(digit *dst, const digit *op1, digit op2, digcnt n) noexcept;

void 	fpz_reciprocal_64(int_reciprocal *ir, digit d2, digit d1, digit d0) noexcept;

digit 	fpz_div128_64(digit n2, digit n1, digit n0, const int_reciprocal* ir) noexcept;
digit 	fpz_divrem128_64(digit n3, digit n2, digit n1, const int_reciprocal* ir, digit *rem) noexcept;
digit 	fpz_adjust_quotient(const int_reciprocal *ir, digit n1, digit n0, digit quot, digit rem) noexcept;

digit 	fpz_div10_128(uint64_t n2, uint64_t n1, uint64_t n0) noexcept;

}

static digit
div192_128(const int_reciprocal *ir, digit n3, digit n2, digit n1, digit n0) noexcept
{
	digit rem;
	digit q = fpz_divrem128_64(n3, n2, n1, ir, &rem);
	return fpz_adjust_quotient(ir, n1, n0, q, rem);
}

// Divide an n-digit integer by a 1-digit integer. Returns magnitude of remainder. 
// digits(quot) and num may be equal. 

static digit
divrem_64(dstptr quot, const digit *num, digit denom, digcnt num_size, uint8_t quot_sign, digcnt quot_cap) noexcept
{
	int_reciprocal ir;
	fpz_reciprocal_64(&ir, denom, 0, 0);

	digit n0 = num[num_size.nd - 1];
	digit q = fpz_div128_64(0, n0, 0, &ir);

	if (quot != nullptr) {
		header quot_hdr = { .sign = quot_sign, .size = num_size - digcnt{1} + digcnt{q != 0} };

		if (FPZ_UNLIKELY(quot_cap < quot_hdr.size)) {
			set_overflow_flag(quot);
			quot = nullptr;	// Stop writing into quot (we still compute the remainder)
		} else {
			if (q != 0) digits(quot)[quot_hdr.size.nd - 1] = q;
			set_header(quot, quot_hdr);
		}
	}

	digit n1 = n0 - q * denom;

	for (digcnt i = num_size - digcnt{2}; i < num_size; --i) {
		n0 = num[i.nd];
		q = fpz_div128_64(n1, n0, 0, &ir);
		n1 = n0 - q * denom;
		if (quot != nullptr) digits(quot)[i.nd] = q;
	}

	return n1;
}

// Divide num by denom (which is a 64-bit built-in signed integer). Returns remainder.
// num and quot are allowed to coincide.

extern "C" int64_t
fpz_divrem_i64(dstptr quot, srcptr num, int64_t denom, digcnt quot_cap) noexcept
{
	header num_hdr = get_header(num);
	uint64_t abs_denom = abs_i64(denom);

	digit rem;

	if (num_hdr.size == digcnt{0} || (num_hdr.size == digcnt{1} && digits(num)[0] < abs_denom)) {
		if (quot != 0) set_i64(quot, 0);
		rem = digits(num)[0];
	} else {
		rem = divrem_64(quot, digits(num), abs_denom, num_hdr.size, 
	    	    static_cast<uint8_t>(num_hdr.sign ^ (denom < 0)), quot_cap);
	}
	
	// Since the remainder is always strictly smaller than the absolute value of the
	// denominator, it can not be INT64_MIN and hence we can negate and still stay within
	// range of an int64_t

	return num_hdr.sign == POSITIVE ? rem : -rem;
}

extern "C" int
fpz_divrem10(dstptr quot, srcptr num, digcnt quot_cap) noexcept
{
	header num_hdr = get_header(num);
	digcnt num_size = num_hdr.size;

	if (num_size == digcnt{0}) {
		if (quot != nullptr) set_i64(quot, 0);
		return 0;
	}

	digit n0 = digits(num)[num_size.nd - 1];
	digit q = fpz_div10_128(0, n0, 0);

	if (quot != nullptr) {
		header quot_hdr = { .sign = num_hdr.sign, .size = num_size - digcnt{1} + digcnt{(q != 0)} };

		if (FPZ_UNLIKELY(quot_cap < quot_hdr.size)) {
			set_overflow_flag(quot);
			quot = nullptr;	// Stop writing into quot
		} else {
			if (q != 0) digits(quot)[quot_hdr.size.nd - 1] = q;
			set_header(quot, quot_hdr);
		}
	}

	digit n1 = n0 - q * 10;

	for (digcnt i = num_size - digcnt{2}; i < num_size; --i) {
		n0 = digits(num)[i.nd];
		q = fpz_div10_128(n1, n0, 0);
		n1 = n0 - q * 10;
		if (quot != nullptr) digits(quot)[i.nd] = q;
	}
	return num_hdr.sign == POSITIVE ? static_cast<int>(n1) : -static_cast<int>(n1);
}

// Divide num by the 128 bit integer denom. None of the pointer parameters can be equal, so
// the caller has to take care of that.

static void
divrem_128(dstptr quot, digit *num, const digit *denom, header num_hdr, header denom_hdr, digcnt quot_cap) noexcept
{
	assert(digits(quot) != num);
	assert(digits(quot) != denom);
	assert(num != denom);

	int_reciprocal ir;
	fpz_reciprocal_64(&ir, denom[1], denom[0], 0);

	digcnt num_size = num_hdr.size;
	num += num_size.nd - 2;
	
	digit q = div192_128(&ir, 0, num[1], num[0], 0);

	if (quot != nullptr) {
		header quot_hdr = { .sign = static_cast<uint8_t>(num_hdr.sign ^ denom_hdr.sign), 
			.size = num_size - digcnt{2} + digcnt{(q != 0)} };

		if (FPZ_UNLIKELY(quot_cap < quot_hdr.size)) {
			set_overflow_flag(quot);
			quot = nullptr;	// Stop writing into quot (we still compute the remainder)
		} else {
			if (q != 0) {
				fpz_usubmul_lowcap_64_self(num, denom, q, digcnt{2});
				digits(quot)[quot_hdr.size.nd - 1] = q;
			}
			set_header(quot, quot_hdr);
		}

		quot += num_hdr.size.nd - 1;
	}

	for (; num_size > digcnt{2}; --num_size) {
		--num;

		q = div192_128(&ir, num[2], num[1], num[0], 0);
		fpz_usubmul_64_self(num, denom, q, digcnt{2});
		if (quot != nullptr) *(--quot) = q;
	}
}

static void
divrem_ge192(dstptr quot, digit *num, const digit *denom, header num_hdr, header denom_hdr, digcnt quot_cap) noexcept
{
	digcnt num_size = num_hdr.size;
	digcnt denom_size = denom_hdr.size;

	int_reciprocal ir;
	fpz_reciprocal_64(&ir, denom[denom_size.nd - 1], denom[denom_size.nd - 2], denom[denom_size.nd - 3]);

	num += num_size.nd - denom_size.nd;

	digit q = div192_128(&ir, 0, num[denom_size.nd - 1], num[denom_size.nd - 2], num[denom_size.nd - 3]); 

	if (q != 0 && fpz_usubmul_lowcap_64_self(num, denom, q, denom_size) != 0) {
		q--;
		fpz_uadd_self(num, denom, denom_size);
	}

	if (quot != nullptr) {
		header quot_hdr = { .sign = static_cast<uint8_t>(num_hdr.sign ^ denom_hdr.sign), 
			.size = num_size - denom_size + digcnt{q != 0} };

		if (FPZ_UNLIKELY(quot_cap < quot_hdr.size)) {
			set_overflow_flag(quot);
			quot = nullptr;	// Stop writing into quot (we still compute the remainder)
		} else {
			if (q != 0) digits(quot)[quot_hdr.size.nd - 1] = q;
			set_header(quot, quot_hdr);
		}
		
		quot += num_size.nd - denom_size.nd + 1;
	}

	for (; num_size > denom_size; --num_size) {
		--num;

		q = div192_128(&ir, num[denom_size.nd], num[denom_size.nd - 1], 
			num[denom_size.nd - 2], num[denom_size.nd - 3]);

		if (fpz_usubmul_64_self(num, denom, q, denom_size) != 0) {
			q--;
			fpz_uadd_self(num, denom, denom_size);
		}

		if (quot != nullptr) *(--quot) = q;
	}
}

static FPZ_FORCE_INLINE bool
is_less_than(srcptr op1, srcptr op2, digcnt n1, digcnt n2) noexcept
{
	if (n1 < n2) return true;
	else if (n1 == n2) return digits(op1)[n1.nd - 1] < digits(op2)[n2.nd - 1];
	return false;
}

// Given numerator (dividend) and denominator (divisor), compute quotient and remainder (i.e. 
// integers quot and rem satisfying num = quot*denom + rem where 0 <= |rem| < |denom|). Either of 
// quot and rem may be null pointers, in which case the integer in question is not written.

extern "C" void
fpz_divrem(dstptr quot, dstptr rem, srcptr num, srcptr denom, digcnt quot_cap, digcnt rem_cap) noexcept
{
	assert(quot != rem);

	header num_hdr = get_header(num);
	header denom_hdr = get_header(denom);

	check_overflow_flag_2(num_hdr, denom_hdr);

	// Handle the trivial case
	if (FPZ_UNLIKELY(num_hdr.size == digcnt{0} || is_less_than(num, denom, num_hdr.size, denom_hdr.size))) {
		if (quot != nullptr) set_i64(quot, 0);
		if (rem != nullptr && rem != num) set(rem, num, rem_cap);
		return;
	}

	if (denom_hdr.size <= digcnt{1}) {
		digit r = divrem_64(quot, digits(num), digits(denom)[0], num_hdr.size, 
		    static_cast<uint8_t>(num_hdr.sign ^ denom_hdr.sign), quot_cap);

		if (rem != nullptr) {
			header rem_hdr = { .sign = num_hdr.sign, .size = digcnt{r != 0} };
			digits(rem)[0] = r;
			set_header(rem, rem_hdr);
		}
		
		return;
	}

	// Temporary storage which the numerator will be copied into; the core division functions will 
	// change this gradually into the correct remainder. If there is space, use rem, since that
	// avoids a copy in the end.

	dstptr rem_tmp = rem;

	// Check if the remainder has enough room to hold the entirety of the numerator. If not, allocate
	// space for it. 

	if (rem == nullptr || rem_cap < num_hdr.size) {
		size_t sz = int_alloc_size(num_hdr.size);
		rem_tmp = static_cast<dstptr>(alloca(sz));
	}

	// If quot or rem_tmp (often the same as rem) coincides with denom, we need to make a copy of denom
	// before doing any calculations (writing into quot and tmp will happen before we're done using denom).

	srcptr denom_tmp = denom;

	if (denom == quot || denom == rem_tmp) {
		size_t sz = int_alloc_size(denom_hdr.size);
		dstptr denom_tmp_0 = static_cast<dstptr>(alloca(sz));
		set(denom_tmp_0, denom, denom_hdr.size);
		denom_tmp = denom_tmp_0;
	}

	if (rem_tmp != num) set(rem_tmp, num, num_hdr.size);
	
	if (denom_hdr.size == digcnt{2}) divrem_128(quot, digits(rem_tmp), digits(denom_tmp), num_hdr, denom_hdr, quot_cap);
	else divrem_ge192(quot, digits(rem_tmp), digits(denom_tmp), num_hdr, denom_hdr, quot_cap);

	if (rem != nullptr) {
		header rem_hdr = { .sign = num_hdr.sign, .size = denom_hdr.size};
		while (rem_hdr.size != digcnt{0} && digits(rem_tmp)[rem_hdr.size.nd - 1] == 0) --rem_hdr.size;

		set_header(rem_tmp, rem_hdr);
		if (rem != rem_tmp) set(rem, rem_tmp, rem_cap);
	}
}

static inline int
fpz_ucmp_n(const digit* op1, const digit* op2, digcnt n) noexcept
{
	for (digcnt i = n - digcnt{1}; i < n; --i) {
		if (op1[i.nd] > op2[i.nd]) return 1;
		else if (op1[i.nd] < op2[i.nd]) return -1;
	}

	return 0;
}

extern "C" long double
fpz_fdiv(srcptr num, srcptr denom) noexcept
{
	bitcnt msb_num = fpz_bit_size(num);
	bitcnt msb_denom = fpz_bit_size(denom);

	if (msb_num == bitcnt{0} && msb_denom == bitcnt{0}) return NAN;	// 0 / 0 
	else if (msb_denom == bitcnt{0}) return INFINITY;		// num / 0, num != 0 
	else if (msb_num == bitcnt{0}) return 0.0;			// 0 / denom, denom != 0

	digcnt num_size = to_digits(msb_num);
	digcnt denom_size = to_digits(msb_denom);

	digcnt size = max(num_size, denom_size) + digcnt{1};
	dstptr num_sh = static_cast<dstptr>(alloca(int_alloc_size(size)));

	bitcnt s = to_bits(size) - bitcnt{1} - msb_num;
	bitcnt t = BITS_PER_DIGIT - bitcnt{1} - remainder_bits(msb_denom);

	// this ensures that num_sh is also normalized when it is shifted by the normalization
	// shift of denom (but num_sh has at least one more digit than denom)
	fpz_shl(num_sh, num, s - t, to_bits(size));

	digit *num_ptr = digits(num_sh) + size.nd - denom_size.nd - 1;
	const digit *denom_ptr = digits(denom);

	// take care of the most significant 1-bit of the significand
	if (fpz_ucmp_n(num_ptr + 1, denom_ptr, denom_size) < 0) {
		// this will shift the most significant bit out of existence
		fpz_ushl_1_self(num_ptr, bitcnt{1}, denom_size + digcnt{1});
		++s;
	}

	fpz_usub_self(num_ptr + 1, denom_ptr, denom_size);

	uint64_t q;
	int_reciprocal ir;

	switch (denom_size.nd) {
	case 1:
		fpz_reciprocal_64(&ir, denom_ptr[0], 0, 0);
		q = fpz_div128_64(num_ptr[1], num_ptr[0], 0, &ir);

		fpz_usubmul_64_self(num_ptr, denom_ptr, q, denom_size);
		break;
	case 2:
		fpz_reciprocal_64(&ir, denom_ptr[1], denom_ptr[0], 0);
		q = div192_128(&ir, num_ptr[2], num_ptr[1], num_ptr[0], 0);

		fpz_usubmul_64_self(num_ptr, denom_ptr, q, denom_size);
		break;
	default:
		fpz_reciprocal_64(&ir, denom_ptr[denom_size.nd - 1], denom_ptr[denom_size.nd - 2], denom_ptr[denom_size.nd - 3]);
		q = div192_128(&ir, num_ptr[denom_size.nd], num_ptr[denom_size.nd - 1], num_ptr[denom_size.nd - 2], num_ptr[denom_size.nd - 3]);

		if (fpz_usubmul_64_self(num_ptr, denom_ptr, q, denom_size) != 0) {
			q--;
			fpz_uadd_self(num_ptr, denom_ptr, denom_size);
		}
		break;
	}

	/* Or in any remaing bits */
	
	int remainder_nonzero = 0;

	for (digcnt i{0}; i < size; ++i) {
		if (num_ptr[i.nd] != 0) {
			remainder_nonzero = 1;
			break;
		}
	}

	bitcnt exp = t - s + to_bits(size - denom_size - digcnt{1});
	int sign = 1 - 2*(get_sign(num) ^ get_sign(denom));

	// Check if we need to round up

	// Is the last bit 0? 
	if (q & 1) {
		q >>= 1;
		// Round up if either next bit is nonzero (round to even) or if remainder is nonzero 
		if ((q & 1) || remainder_nonzero) {
			q += 1;
			// If bit 63 just got set, we've "overflowed", so shift right and increment exponent
			if (q & 0x8000000000000000) {
				q >>= 1;
				++exp;
			}
		}
	} else {
		q >>= 1;
	}

	q |= 0x8000000000000000;

	return sign*ldexpl(static_cast<long double>(q), (exp + bitcnt{1}).nb);
}

} // namespace fpz::core
