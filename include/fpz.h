#ifndef FPZ_H
#define FPZ_H

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

#include <algorithm>
#include <compare>
#include <initializer_list>
#include <numbers>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include "fpz_core.h"

namespace fpz {
namespace detail {

using core::digit, core::word;
using core::dstptr, core::srcptr;
using core::size_type, core::bitcnt, core::digcnt;

using core::BITS_PER_DIGIT;
using core::POSITIVE, core::NEGATIVE;

// Estimate the number of bits required to hold an integer represented
// as a string of digits, optionally starting with a sign (+/-) and containing
// any number of spaces (only ' ' is allowed, no newlines/tabs) and digits. 
// Only decimal and hexadecimal are supported, and decimal is assumed unless the 
// first characters after the sign are '0' and 'x'. The result is exact in 
// the hexadecimal case. Both upper and lower case hexadecimal digits are allowed.

constexpr bitcnt 
estimate_bit_size(const char *str)
{
	unsigned base = 10;

	if (str[0] == '+' || str[0] == '-') str++;
	if (str[0] == '0' && str[1] == 'x') {
		base = 16;
		str+= 2;
	}

	// Count digit characters and compute the magnitude as we go along. 
	// Overflow is not a concern at this time; the value will only be 
	// used once we know (by other means) that it fits in 64 bits or less.

	uint64_t magnitude = 0;
	uint32_t num_digits = 0;

	while (core::is_digit(*str, base) || *str == ' ') {
		if (*str != ' ') {
			++num_digits;
			magnitude = magnitude*base + core::digit_value(*str, base);
		}
		str++;
	}

	if (*str != '\0') throw core::invalid_argument{};

	constexpr const char max_u63_hex[] = "7FFFFFFFFFFFFFFF";
	constexpr const char max_u63_dec[] = "9223372036854775807";

	// The most important thing here is to determine if the magnitude can be represented
	// in (strictly) less than 64 bits. This is definitely not the case if the string has 
	// more digits than max_u63_hex or max_u63_dec (according to base), so in these cases
	// we estimate (fairly accurately, should be at most 4 bits too high in most cases) 
	// without actually computing the value. If num_digits is less than or equal to the 
	// number of digits in these strings (according to base), we know the value fits in 64 bits, 
	// so we get our answer from the computed magnitude.

	if (base == 16 && num_digits > sizeof(max_u63_hex) - 1) {
		return bitcnt { num_digits * 4 };
	} else if (base == 10 && num_digits > sizeof(max_u63_dec) - 1) {
		constexpr long double log2_10 = std::numbers::ln10_v<long double> / 
			std::numbers::ln2_v<long double>;
		return bitcnt{ static_cast<core::size_type>(num_digits * log2_10) + 1 };
	} else { 
		return core::bit_size(magnitude);
	}
}

// Compute the number of bits required to hold an integer with the given
// digits. 

constexpr bitcnt
bit_size(std::initializer_list<digit> list)
{
	size_t first = 0;
	const digit* digits = list.begin();
	while (first < list.size() && digits[first] == 0) first++;

	if (first == list.size()) return bitcnt{0};

	bitcnt ret{core::bit_size(digits[first])};
	
	while (++first < list.size()) ret += bitcnt{BITS_PER_DIGIT};

	return ret;
}

// Check if an unsigned integer is a power of two.

constexpr bool
is_power_of_2(size_t x)
{
	return (x != 0) && ((x & (x - 1)) == 0);
}

// Integers strictly smaller than 64 bits are treated specially, so
// test for that in a uniform way. 

constexpr bool
is_large(bitcnt prec)
{
	return prec >= BITS_PER_DIGIT;
}

// Compute the ceiling of the binary logarithm of an unsigned integer.

constexpr bitcnt
ceil_log2(size_t x)
{
	if (x == 0) throw std::domain_error("log2(0) is undefined");

	bitcnt ret = core::bit_size(x);
	return is_power_of_2(x) ? ret - bitcnt{1} : ret;
}

// Make sure we know about it if an operation returns an invalid result
// (which is indicative of a bug in the library). The check is always performed
// in the constant evaluated case, otherwise it is only performed in debugging builds.

template <typename T> constexpr void
assert_valid(T op, bitcnt prec, const char *filename, int line)
{
	bool valid = true;
	if (std::is_constant_evaluated()) {
		if constexpr (std::is_same_v<T, int64_t>) valid = core::bit_size(op) <= prec;
		else valid = core::is_valid(op, prec);
		if (!valid) throw std::logic_error("operation returned an invalid result");
	} else {
#ifndef NDEBUG
		if constexpr (std::is_same_v<T, int64_t>) valid = core::bit_size(op) <= prec;
		else valid = core::fpz_is_valid(op, prec);
		if (!valid) core::fpz_assertion_failure("operation returned an invalid result", filename, line);
#endif
	}
}

#define FPZ_ASSERT_VALID(op, prec)	::fpz::detail::assert_valid(op, prec, __FILE__, __LINE__)

// Structure used to keep track of precision bounds of expressions (at compile-time). It expresses 
// a bound on the magnitude of an expression as m*2^(exp), where both m (multiplier) and 
// exp (exponent) are positive integers. In principle, this could have been a simple bit count,
// but that would result in linear rather than logarithmic growth of the resulting integer variable
// precisions. This, to be fair, is probably not that important for reasonable expressions that we 
// can write by hand, but if we're computing a sum of 256 terms of 128 bits each in a loop, the result
// would suddenly require 384 bits, as opposed to 136. Such a calculation could be done differently
// from the rest (it is already a separate function, sum(expr_bound b, size_t count)), but that would
// introduce a strange inconsistency which might cause trouble in certain situations.

struct expr_bound {
	core::bitcnt 	exp;
	size_t 		m;

	consteval
	operator bitcnt() const noexcept
	{
		return exp + ceil_log2(m);
	}
};

// Define arithmetic operators for expression bounds, providing bounds on the
// corresponding arithmetic expressions. In general, no overflow checking (for the 
// expr_bound fields) is performed, since we expect to be working with reasonably
// simple expressions involving variables of no more than a few thousand bits.

// Return a bound on e1 + e2 if e1 is bounded by b1 and e2 is bounded by b2.

consteval expr_bound
operator+(expr_bound b1, expr_bound b2) noexcept
{
	// Make sure b1 has the largest exponent
	if (b1.exp < b2.exp) std::swap(b1, b2);

	bitcnt shift = b1.exp - b2.exp;

	// m1*2^e + m2*2^e = (m1 + m2)*2^e
	if (shift == bitcnt{0}) return { .exp = b1.exp, .m = b1.m + b2.m };

	// Convert the bound with the smallest exponent into one with the same
	// exponent as the other bound.

	// If the shift is larger than the with of a size_t, the converted bound
	// will always have multiplier equal to 1

	constexpr int max_bits = std::numeric_limits<size_t>::digits;
	if (shift >= bitcnt{ static_cast<size_type>(max_bits) })
		return { .exp = b1.exp, .m = b1.m + 1 };
	
	size_t mask = (size_t{1} << shift.nb) - 1;
	size_t m = (b2.m >> shift.nb) + ((b2.m & mask) != 0);

	return { .exp = b1.exp, .m = b1.m + m };
}

// Return a bound on e1 - e2 if e1 is bounded by b1 and e2 is bounded by b2 (same as e1 + e2).

consteval expr_bound 
operator-(expr_bound b1, expr_bound b2) noexcept
{
	return b1 + b2;
}

// Return a bound on e1 * e2 if e1 is bounded by b1 and e2 is bounded by b2.

consteval expr_bound 
operator*(expr_bound b1, expr_bound b2)
{
	return { .exp = b1.exp + b2.exp, .m = b1.m*b2.m };
}

// Return a bound on the left shift e << bitcnt if e is bounded by b.

consteval expr_bound 
operator<<(expr_bound b, bitcnt nbits)
{
	return { .exp = b.exp + nbits, .m = b.m };
}

// Return a bound on the right shift e >> bitcnt if e is bounded by b.

consteval expr_bound 
operator>>(expr_bound b, bitcnt nbits)
{
	if (b.exp > nbits) return { .exp = b.exp - nbits, .m = b.m };
	else return { .exp = bitcnt{0}, .m = 1 };
}

// Return a bound on e1 / e2 if e1 is bounded by n and e2 is bounded by d.

consteval expr_bound 
operator/(expr_bound n, expr_bound d)
{
	return n;
}

// Return a bound on e1 % e2 if e1 is bounded by n and e2 is bounded by d.

consteval expr_bound 
operator%(expr_bound n, expr_bound d)
{
	// This appears to be equivalent to std::min(n, d), but this is clearer
	if (static_cast<bitcnt>(n) <= static_cast<bitcnt>(d)) return n;
	else return d;
}

// Return a bound on count*e if e is bounded by b.

consteval expr_bound 
sum(expr_bound b, size_t count)
{
	return { .exp = b.exp, .m = count*b.m };
}

// Return a bound on e^count if e is bounded by b.

consteval expr_bound 
product(expr_bound b, size_t count)
{
	size_t m = b.m;
	for (size_t i = 0; i < count; i++) m *= b.m;

	return { .exp = bitcnt{ static_cast<size_type>(count*b.exp.nb) }, .m = m };
}

} // namespace detail

// External type for doing precision calculation. Users aren't expected to deal
// with the expr_bound structure directly, except implictly through the operators 
// defined above. The results returned from such operators are supposed to be 
// initialize constants of type precision instead, and those constants give
// rise to integer types when passed to the ztype template (they are implicitly
// convertible to size_type). Such initialization requires that we "collapse"
// the expr_bound into a 1-term bound (i.e. we add ceil_log2(m) to the value 
// of exp and set m to 1), because otherwise, we may end up with inconsistent results.
// This is probably best explained through an example (assume a hypothetical 
// implementation where the _bits literal returns an expr_bound and where expr_bound
// implicitly converts to uint32_t):
//
// constexpr expr_bound base = 128_bits;
// constexpr expr_bound b1 = base*base + base;
// constexpr expr_bound b2 = b1*b1 + (b1 - b1)*b1;
//
// ztype<b1> x = ...;
// ztype<b1> y = ...;
// ztype<b2> z = x*x + (x - y)*x;
//
// We expect that z has enough room to hold the result of the expression we assign
// to it (it has the same form as what we assign to b2), but that might not be the case; 
// x and y has the precision given by collapsing b1, while b2 is computed from the 
// uncollapsed bound b1 and that may end up being smaller. The result is that a static
// assert would fire down the road for seemingly no reason (after all, the user did everything
// right). Using the type precision instead fixes the issue, since the collapse is done on 
// initialization (we don't want to make the user do this explicitly):
// 
// constexpr precision base = 128_bits;
// constexpr precision b1 = base*base + base;
// constexpr precision b2 = b1*b1 + (b1 - b1)*b1;
// 
// Users could still do 
//
// constexpr auto b1 = base*base + base;
//
// and get an expr_bound, but since expr_bound doesn't implicitly convert to
// an unsigned integer, at least the program wouldn't compile. 

struct precision : public detail::expr_bound {
	explicit consteval 
	precision(core::bitcnt e) noexcept : expr_bound{ .exp = e, .m = 1 } {}
	
	consteval
	precision(const detail::expr_bound& b) noexcept : 
	    expr_bound{ .exp = static_cast<core::bitcnt>(b), .m = 1 } {}

	consteval
	operator core::size_type() const noexcept
	{
		return exp.nb;
	}
};

consteval precision
operator""_bits(unsigned long long count)
{
	if (count > std::numeric_limits<core::size_type>::max()) throw std::range_error("bit count is too large");
	return precision{ core::bitcnt{static_cast<uint32_t>(count)} };
}

namespace detail {

/* Operator structures.  */

struct op_add {
	static constexpr int64_t neutral_value = 0;

	static constexpr expr_bound 
	bound(expr_bound b1, expr_bound b2)
	{
		return b1 + b2;
	}

	// bound for operation repeated count times

	static constexpr expr_bound 
	bound_rep(expr_bound b, size_t count)
	{
		return sum(b, count);
	}

	static constexpr bitcnt
	term_prec(bitcnt bp, size_t count)
	{
		return bp > ceil_log2(count) ? bp - ceil_log2(count) : bitcnt{0};
	}

	static constexpr int64_t 
	eval(int64_t op1, int64_t op2, bitcnt bit_bound) 
	{
		op1 += op2;
		FPZ_ASSERT_VALID(op1, bit_bound);
		return op1;
	}

	static constexpr void 
	eval(dstptr dst, int64_t op1, int64_t op2, bitcnt prec)
	{ 
		digcnt cap = core::to_digits(prec);

		if (std::is_constant_evaluated()) core::add_i64_i64(dst, op1, op2, cap);
		else core::fpz_add_i64_i64(dst, op1, op2, cap);
		FPZ_ASSERT_VALID(dst, prec);
	}

	static constexpr void 
	eval(dstptr dst, srcptr op1, int64_t op2, bitcnt prec) 
	{ 
		digcnt cap = core::to_digits(prec);

		if (std::is_constant_evaluated()) core::add_i64(dst, op1, op2, cap);
		else core::fpz_add_i64(dst, op1, op2, cap);
		FPZ_ASSERT_VALID(dst, prec);
	}

	static constexpr void 
	eval(dstptr dst, int64_t op1, srcptr op2, bitcnt prec)
	{
		eval(dst, op2, op1, prec);
	}

	static constexpr void 
	eval(dstptr dst, srcptr op1, srcptr op2, bitcnt prec)
	{
		digcnt cap = core::to_digits(prec);
		if (std::is_constant_evaluated()) core::add(dst, op1, op2, cap);
		else core::fpz_add(dst, op1, op2, cap);
		FPZ_ASSERT_VALID(dst, prec);
	}
};

struct op_sub {

	static constexpr int64_t neutral_value = 0;

	static constexpr expr_bound 
	bound(expr_bound b1, expr_bound b2)
	{
		return b1 - b2;
	}

	// bound for operation repeated count times

	static constexpr expr_bound 
	bound_rep(expr_bound b, size_t count)
	{
		return sum(b, count);
	}

	static constexpr bitcnt
	term_prec(bitcnt bp, size_t count)
	{
		return bp > ceil_log2(count) ? bp - ceil_log2(count) : bitcnt{0};
	}

	static constexpr int64_t 
	eval(int64_t op1, int64_t op2, bitcnt bit_bound) 
	{ 
		op1 -= op2;
		FPZ_ASSERT_VALID(op1, bit_bound);
		return op1;
	}

	static constexpr void 
	eval(dstptr dst, int64_t op1, int64_t op2, bitcnt prec)
	{
		if (std::is_constant_evaluated()) core::sub_i64_i64(dst, op1, op2);
		else core::fpz_sub_i64_i64(dst, op1, op2);
		FPZ_ASSERT_VALID(dst, prec);
	}

	static constexpr void 
	eval(dstptr dst, srcptr op1, int64_t op2, bitcnt prec) 
	{ 
		digcnt cap = core::to_digits(prec);

		if (std::is_constant_evaluated()) core::sub_i64(dst, op1, op2, cap);
		else core::fpz_sub_i64(dst, op1, op2, cap);
		FPZ_ASSERT_VALID(dst, prec);
	}

	static constexpr void 
	eval(dstptr dst, int64_t op1, srcptr op2, bitcnt prec)
	{
		digcnt cap = core::to_digits(prec);

		if (std::is_constant_evaluated()) core::i64_sub(dst, op1, op2, cap);
		else core::fpz_i64_sub(dst, op1, op2, cap);
		FPZ_ASSERT_VALID(dst, prec);
	}

	static constexpr void 
	eval(dstptr dst, srcptr op1, srcptr op2, bitcnt prec)
	{
		digcnt cap = core::to_digits(prec);

		if (std::is_constant_evaluated()) core::sub(dst, op1, op2, cap);
		else core::fpz_sub(dst, op1, op2, cap);
		FPZ_ASSERT_VALID(dst, prec);
	}
};

struct op_mul {
	static constexpr int64_t neutral_value = 1;

	static constexpr expr_bound 
	bound(expr_bound b1, expr_bound b2)
	{
		return b1 * b2;
	}

	static constexpr expr_bound 
	bound_rep(expr_bound b, size_t count)
	{
		return product(b, count);
	}

	static constexpr bitcnt
	term_prec(bitcnt bp, size_t count)
	{
		return bitcnt{ static_cast<size_type>(bp.nb / count) };
	}

	static constexpr int64_t 
	eval(int64_t op1, int64_t op2, bitcnt prec) 
	{ 
		op1 *= op2;
		FPZ_ASSERT_VALID(op1, prec);
		return op1;
	}

	static constexpr void 
	eval(dstptr dst, int64_t op1, int64_t op2, bitcnt prec) 
	{
		digcnt cap = core::to_digits(prec);

		if (std::is_constant_evaluated()) core::mul_i64_i64(dst, op1, op2, cap);
		else core::fpz_mul_i64_i64(dst, op1, op2, cap);
		FPZ_ASSERT_VALID(dst, prec);
	}

	static constexpr void 
	eval(dstptr dst, srcptr op1, int64_t op2, bitcnt prec)
	{
		digcnt cap = core::to_digits(prec);

		if (std::is_constant_evaluated()) core::mul_i64(dst, op1, op2, cap);
		else core::mul_i64(dst, op1, op2, cap);
		FPZ_ASSERT_VALID(dst, prec);
	}

	static constexpr void 
	eval(dstptr dst, int64_t op1, srcptr op2, bitcnt prec)
	{
		eval(dst, op2, op1, prec);
	}

	static constexpr void 
	eval(dstptr dst, srcptr op1, srcptr op2, bitcnt prec)
	{
		digcnt cap = core::to_digits(prec);		

		if (std::is_constant_evaluated()) core::mul(dst, op1, op2, cap);
		else core::fpz_mul(dst, op1, op2, cap);
		FPZ_ASSERT_VALID(dst, prec);
	}
};

struct op_div {
	static constexpr expr_bound 
	bound(expr_bound n, expr_bound d)
	{
		return n / d;
	}

	static constexpr int64_t 
	eval(int64_t op1, int64_t op2, bitcnt prec)
	{
		op1 /= op2;
		FPZ_ASSERT_VALID(op1, prec);
		return op1;
	}

	static constexpr int64_t 
	eval(int64_t op1, srcptr op2, bitcnt prec)
	{
		if (std::is_constant_evaluated()) core::i64_divrem(&op1, op1, op2);
		else core::fpz_i64_divrem(&op1, op1, op2);

		FPZ_ASSERT_VALID(op1, prec);
		return op1;
	}

	static constexpr void 
	eval(dstptr dst, srcptr op1, int64_t op2, bitcnt prec)
	{
		digcnt cap = core::to_digits(prec);

		if (std::is_constant_evaluated()) core::divrem_i64(dst, op1, op2, cap);
		else core::fpz_divrem_i64(dst, op1, op2, cap);
		FPZ_ASSERT_VALID(dst, prec);
	}

	static constexpr void 
	eval(dstptr dst, srcptr op1, srcptr op2, bitcnt prec)
	{
		digcnt cap = core::to_digits(prec);
	
		if (std::is_constant_evaluated()) core::divrem(dst, nullptr, op1, op2, cap, digcnt{0});
		else core::fpz_divrem(dst, nullptr, op1, op2, cap, digcnt{0});
		FPZ_ASSERT_VALID(dst, prec);
	}
};

struct op_mod {
	static constexpr expr_bound 
	bound(expr_bound n, expr_bound d)
	{
		return n % d;
	}

	static constexpr int64_t 
	eval(int64_t op1, int64_t op2, bitcnt prec)
	{
		op1 %= op2;
		FPZ_ASSERT_VALID(op1, prec);
		return op1;
	}

	static constexpr int64_t 
	eval(int64_t op1, srcptr op2, bitcnt prec)
	{
		int64_t q{0};
		if (std::is_constant_evaluated()) op1 = core::i64_divrem(&q, op1, op2);
		else op1 = core::fpz_i64_divrem(&q, op1, op2);

		FPZ_ASSERT_VALID(op1, prec);
		return op1;
	}

	static constexpr int64_t 
	eval(srcptr op1, int64_t op2, bitcnt prec)
	{
		if (std::is_constant_evaluated()) op2 = core::divrem_i64(nullptr, op1, op2, digcnt{0});
		else op2 = core::fpz_divrem_i64(nullptr, op1, op2, digcnt{0});
		FPZ_ASSERT_VALID(op2, prec);
		return op2;
	}

	static constexpr void 
	eval(dstptr dst, srcptr op1, srcptr op2, bitcnt prec)
	{
		digcnt cap = core::to_digits(prec);		

		if (std::is_constant_evaluated()) core::divrem(nullptr, dst, op1, op2, digcnt{0}, cap);
		else core::fpz_divrem(nullptr, dst, op1, op2, digcnt{0}, cap);
		FPZ_ASSERT_VALID(dst, prec);
	}
};

struct op_shift_left {
	static constexpr expr_bound 
	bound(expr_bound b, bitcnt nbits)
	{
		return b << nbits;
	}

	static constexpr int64_t 
	eval(int64_t op, bitcnt nbits, bitcnt prec)
	{
		if (std::is_constant_evaluated()) {
			op *= core::digit{1} << nbits.nb;
		} else {
			register uint8_t nbits_tmp asm ("cl") = nbits.nb;
			__asm__(
				"salq %[nbits], %[op]"
				: [op] "+rm" (op)
				: [nbits] "r" (nbits_tmp)
				: "cc"
			);
		}
		FPZ_ASSERT_VALID(op, prec);
		return op;
	}

	static constexpr void 
	eval(dstptr dst, int64_t op, bitcnt nbits, bitcnt prec)
	{
		if (std::is_constant_evaluated()) core::shl_i64(dst, op, nbits, prec);
		else core::fpz_shl_i64(dst, op, nbits, prec);
		FPZ_ASSERT_VALID(dst, prec);
	}

	static constexpr void 
	eval(dstptr dst, srcptr op, bitcnt nbits, bitcnt prec)
	{
		if (std::is_constant_evaluated()) core::shl(dst, op, nbits, prec);
		else core::fpz_shl(dst, op, nbits, prec);
		FPZ_ASSERT_VALID(dst, prec);
	}
};

struct op_shift_right {

	static constexpr expr_bound 
	bound(expr_bound b, bitcnt nbits)
	{
		return b >> nbits;
	}

	static constexpr int64_t 
	eval(int64_t op, bitcnt nbits, bitcnt prec) 
	{
		if (std::is_constant_evaluated()) op /= uint64_t{1} << nbits.nb;
		else {
			register uint8_t nbits_tmp asm ("cl") = nbits.nb;
			__asm__(
				"sarq %[nbits], %[op]"
				: [op] "+rm" (op)
				: [nbits] "r" (nbits_tmp)
				: "cc"
			);
		}
		FPZ_ASSERT_VALID(op, prec);
		return op;
	}

	static constexpr int64_t 
	eval(srcptr op, bitcnt nbits, bitcnt prec)
	{
		digit tmp[2]{0, 0};
		if (std::is_constant_evaluated()) core::shr(tmp, op, nbits, prec);
		else core::fpz_shr(tmp, op, nbits, prec);
		FPZ_ASSERT_VALID(tmp, prec);

		return (core::get_sign(op) == POSITIVE) ? tmp[1] : -tmp[1];
	}

	static constexpr void 
	eval(dstptr dst, srcptr op, bitcnt nbits, bitcnt prec)
	{
		if (std::is_constant_evaluated()) core::shr(dst, op, nbits, prec);
		else core::fpz_shr(dst, op, nbits, prec);
		FPZ_ASSERT_VALID(dst, prec);
	}	
};

struct op_mod_2exp {

	static constexpr expr_bound 
	bound(expr_bound b, bitcnt nbits)
	{
		return expr_bound{ .exp = std::min(b.exp, nbits), .m = b.m};
	}

	static constexpr int64_t
	eval(int64_t op, bitcnt nbits, bitcnt prec)
	{
		digit mask = (core::digit{1} << nbits.nb) - 1;
		if (op >= 0) op &= mask;
		else op = -((-op) & mask);
		FPZ_ASSERT_VALID(op, prec);
		return op;
	}

	static constexpr int64_t
	eval(srcptr op, bitcnt nbits, bitcnt prec)
	{
		core::digit magnitude = core::digits(op)[0];
		core::digit mask = (core::digit{1} << nbits.nb) - 1;
		magnitude &= mask;

		int64_t r = core::get_sign(op) == POSITIVE ? magnitude : -magnitude;
		FPZ_ASSERT_VALID(r, prec);
		return r;
	}

	static constexpr void
	eval(dstptr dst, srcptr op, bitcnt nbits, bitcnt prec) noexcept
	{
		if (std::is_constant_evaluated()) core::truncate(dst, op, nbits);
		else core::fpz_truncate(dst, op, nbits);
		FPZ_ASSERT_VALID(dst, prec);
	}

};

template <typename T> struct expr_traits_nocvref;

template <typename T>
using expr_traits = expr_traits_nocvref<std::remove_cvref_t<T>>;

template <typename T>
concept Expr = expr_traits<T>::is_expr;

template <typename T>
concept SimpleExpr = expr_traits<T>::is_expr && !expr_traits<T>::is_compound;

template <SimpleExpr E> constexpr const auto 
get_value(const E& e) noexcept
{
	if constexpr (expr_traits<E>::has_storage) return e.value;
	else if constexpr (std::is_same_v<decltype(e.z), const int64_t>) return e.z;
	else return e.z.value;
}

// Internal types used as destinations during expression evaluation. local_var is used
// to pass a writable pointer to the value of a ztype to the evaluation function
// of an expression (this happens inside the assignment operator for ztype, when
// the right hand side is an expression). tmp_var is used to hold a pointer to a  
// temporary integer (using ztype itself would require making the various expression
// types friends with ztype, since ztype keeps its value private). 
//
// We use separate types in order to detect the possibility of aliasing, which may
// impose conditions on evaluation (a tmp_var is assumed to be allocated by the 
// expression evaluator, and thus can't alias with any of the integers that make up
// the expression).

template <core::size_type prec_> 
struct local_var {
	static constexpr bitcnt prec{prec_};
	dstptr ptr;
};

template <core::size_type prec_>
struct tmp_var {
	static constexpr bitcnt prec{prec_};
	dstptr ptr;
};

template <Expr E> 
struct tmp_space {
	static constexpr bitcnt prec{expr_traits<E>::bit_bound};
	word 	value[core::to_digits(prec).nd + 1];

	constexpr tmp_var<prec.nb> 
	pass() 
	{ 
		return { value }; 
	}
};

// Evaluate both simple and compound expressions 

template <Expr E, typename D> FPZ_FORCE_INLINE constexpr const auto
eval_any(const E& e, D dst) noexcept
{
	if constexpr (expr_traits<E>::is_compound) return e.eval(dst);
	else return get_value(e);
}

// Negation

template <Expr E> 
struct negate_expr {
	const E& e;

	template <typename D> FPZ_FORCE_INLINE constexpr auto 
	eval(D dst) const noexcept
	{
		using ETr = expr_traits<E>;

		if constexpr (is_large(ETr::bit_bound)) {
			auto src = eval_any(e, dst);
			if constexpr (!ETr::is_compound) {
				if (std::is_constant_evaluated()) core::set(dst.ptr, src, to_digits(D::prec));
				else core::fpz_set(dst.ptr, src, to_digits(D::prec));
			}
			core::negate(dst.ptr);
			return dst.ptr;
		} else {
			return -eval_any(e, nullptr);
		}
	}
};

// Absolute value 

template <Expr E>
struct abs_expr {
	const E& e;

	template <typename D> FPZ_FORCE_INLINE constexpr auto 
	eval(D dst) const noexcept
	{
		using ETr = expr_traits<E>;

		if constexpr (is_large(ETr::bit_bound)) {
			auto src = eval_any(e, dst);
			if constexpr (!ETr::is_compound) {
				if (std::is_constant_evaluated()) core::set(dst.ptr, src, to_digits(D::prec));
				else core::fpz_set(dst.ptr, src, to_digits(D::prec));
			}
			core::set_sign(dst.ptr, core::POSITIVE);
			return dst.ptr;
		} else {
			return core::abs_i64(eval_any(e, nullptr));
		}
	}
};

template <Expr E1, Expr E2, typename Op>
struct binary_expr {
	const E1& 	e1;
	const E2&	e2;

	template <typename D> FPZ_FORCE_INLINE constexpr auto 
	eval(D dst) const noexcept
	{
		using Tr = expr_traits<binary_expr>;
		using E1Tr = expr_traits<E1>;
		using E2Tr = expr_traits<E2>;

		tmp_space<E1> tmp1;
		tmp_space<E2> tmp2;

		if constexpr (is_large(Tr::bit_bound)) {
			constexpr uint32_t tmp1_bits_used = E1Tr::is_compound ? E1Tr::bit_bound.nb : 0;
			constexpr uint32_t tmp2_bits_used = E2Tr::is_compound ? E2Tr::bit_bound.nb : 0;

			// If both e1 and e2 may contain the integer we're writing into, we can never safely
			// write into dst before we have evaluated both subexpressions, so we have to use 
			// temporaries. The same is the case if dst is too small to hold the result of 
			// evaluating either of the subexpressions. The same call is also performed
			// in the case of two simple subexpressions, but the temporaries won't be used in that case.
			// Otherwise, we're safe as long as our first evaluation goes into a temporary,
			// so choose the order that uses the least memory.

			if constexpr ((E1Tr::template subexpr_may_alias_with<D> && E2Tr::template subexpr_may_alias_with<D>) || 
			    (!E1Tr::is_compound && !E2Tr::is_compound) || 
			    (D::prec < std::min(E1Tr::bit_bound, E2Tr::bit_bound))) {
				Op::eval(dst.ptr, eval_any(e1, tmp1.pass()), eval_any(e2, tmp2.pass()), D::prec);
			} else if constexpr (D::prec < E1Tr::bit_bound || 
			    (D::prec >= E2Tr::bit_bound && tmp1_bits_used <= tmp2_bits_used)) {
				auto src1 = eval_any(e1, tmp1.pass());
				Op::eval(dst.ptr, src1, eval_any(e2, dst), D::prec);
			} else {
				auto src2 = eval_any(e2, tmp2.pass());
				Op::eval(dst.ptr, eval_any(e1, dst), src2, D::prec);
			}
			return dst.ptr;
		} else {
			return Op::eval(eval_any(e1, tmp1.pass()), eval_any(e2, tmp2.pass()), Tr::bit_bound);
		}
	}
};

// Right and left shifts (a.k.a. multiplication and division by powers of two, mul2exp and div2exp)
// as well as mod2exp fall into this expression category (mod2exp doesn't really involve any shifting, 
// but it has the same parameter types as the shifts).

template <Expr E, typename Op, uint32_t nbits>
struct shift_expr {
	const E& e;

	template <typename D> FPZ_FORCE_INLINE constexpr auto 
	eval(D dst) const noexcept
	{
		using Tr = expr_traits<shift_expr>;
		using ETr = expr_traits<E>;

		tmp_space<E> tmp;

		if constexpr (is_large(Tr::bit_bound)) {
			if constexpr (D::prec < ETr::bit_bound) {
				Op::eval(dst.ptr, eval_any(e, tmp.pass()), bitcnt{nbits}, D::prec);
			} else {
				Op::eval(dst.ptr, eval_any(e, dst), bitcnt{nbits}, D::prec);
			}
			return dst.ptr;
		} else {
			return Op::eval(eval_any(e, tmp.pass()), bitcnt{nbits}, Tr::bit_bound);
		} 
	}
};

} // namespace detail

template <core::size_type prec> 
class ztype;

// Constant 64-bit integer

template <int64_t v>
struct const_i64_type {
	static constexpr int64_t z = v;
};

namespace detail {

// Classes supporting a range-based for loop for computing sums or products incrementally

template <size_t count, typename Op, uint32_t dst_prec, uint32_t term_prec>
class accumulate_iter {
	size_t i = 0;
	ztype<dst_prec>& 	dst;
	ztype<term_prec>& 	term;
public:
	constexpr 
	accumulate_iter(ztype<dst_prec>& d, ztype<term_prec>& t) noexcept : dst(d), term(t) {}

	constexpr bool 
	operator!=(size_t cnt) const {
		return i != cnt;
	}

	constexpr accumulate_iter& 
	operator++() noexcept
	{
		// Accumulate current term if there is at least one more iteration 
		// (the last one is handled by the destructor)
		if (FPZ_LIKELY(++i < count)) Op::eval(dst.value, get_value(dst), get_value(term), dst.prec);
		return *this;
	}

	constexpr std::pair<ztype<term_prec>&, const ztype<dst_prec>&>
	operator*() noexcept
	{
		return { term = const_i64_type<Op::neutral_value>{}, dst };
	}
};

struct empty{};

template <core::size_type dst_prec, size_t count, typename Op, Expr E, typename Termtype>
class accumulate_range {
	ztype<dst_prec>& dst;

	// Accumulate into a temporary rather than the destination variable. This prevents
	// the destination from overflowing due to modification during the loop. It also
	// means that the destination variable can be freely used during the loop. 

	ztype<dst_prec> dst_tmp;

	// Allow users to pass their own term variable (which may save some stack memory under 
	// some circumstances, provided the variable is used for other things). If one is not passed,
	// make one.

	std::conditional_t<std::is_same_v<Termtype, empty&>, 
		ztype<Op::term_prec(bitcnt{dst_prec}, count).nb>, Termtype&> term;

	static_assert(count > 0, "iteration count can't be 0");
	static_assert(expr_traits<decltype(term)>::has_storage, "inappropriate term type provided");
	static_assert(core::bitcnt{dst_prec} >= Op::bound_rep(expr_traits<E>::bound, count), "precision of destination is too small");
public:
	constexpr 
	accumulate_range(ztype<dst_prec>& d, const E& init_expr, Termtype&& t) noexcept : 
	    dst(d), dst_tmp(init_expr), term(t) {}

	constexpr 
	~accumulate_range() noexcept
	{
		Op::eval(dst.value, get_value(dst_tmp), get_value(term), dst.prec);
	}
	
	constexpr auto
	begin() noexcept
	{
		return accumulate_iter<count, Op, dst_prec, 
		    std::remove_reference_t<decltype(term)>::prec.nb>{dst_tmp, term};
	}

	constexpr size_t 
	end() const noexcept
	{
		return count;
	}
};

} // namespace detail 


template <core::size_type prec_> class ztype {
public:
	static constexpr core::bitcnt prec{prec_};
	static constexpr uint64_t max_abs_value = detail::is_large(prec) ? 0 : (uint64_t{1} << prec.nb) - 1;
	
private:
	static_assert(prec > core::bitcnt{0}, "precision must be at least 1 bit");

	// Integer value; either an array of sufficiently many digits to hold prec bits + header,
	// or a built-in 64-bit signed integer (if prec is 63 bits or less)
	std::conditional_t<detail::is_large(prec), 
		core::word[core::to_digits(prec).nd + 1], int64_t> value;

	constexpr ztype& 
	set_i64_nothrow(int64_t v) noexcept
	{
		if constexpr (detail::is_large(prec)) {
			core::set_i64(value, v);
		} else {
			value = v;
		}
		return *this;
	}

	ztype(int64_t v) noexcept
	{
		(void)set_i64_nothrow(v);
	}

	template <size_t count, typename Op, detail::Expr E, typename Termtype> FPZ_FORCE_INLINE auto 
	accumulate(const E& init_expr, Termtype&& term) noexcept
	{
		return detail::accumulate_range<prec.nb, count, Op, E, Termtype>(*this, init_expr, term);
	}

public:

	// Default constructor; do nothing in runtime contexts, but zero-initialize in 
	// constexpr contexts
	constexpr 
	ztype() noexcept
	{
		if (std::is_constant_evaluated()) {
			if constexpr (detail::is_large(prec)) {
				for (core::word& w : value) w = 0;
			} else {
				value = 0;
			}
		}
	}

	// Construct from arbitrary (statically bounded) expression
	template <detail::Expr E> explicit FPZ_FORCE_INLINE constexpr
	ztype(const E& e) noexcept
	{
		(void)operator=(e);
	}

	template <detail::Expr E> FPZ_FORCE_INLINE constexpr ztype&
	operator=(const E& rhs) noexcept
	{
		using detail::expr_traits, detail::is_large, detail::local_var;
		using detail::get_value, detail::eval_any, core::to_digits;

		static_assert(prec >= expr_traits<E>::bit_bound, "target precision is too small");

		if constexpr (is_large(expr_traits<E>::bit_bound)) {
			if constexpr (expr_traits<E>::is_compound) {
				// Pass a wrapped, writable pointer to our digits 
				// to the evaluation method for the expression
				local_var<prec.nb> dst{value};
				rhs.eval(dst);
			} else {
				if (std::is_constant_evaluated()) core::set(value, get_value(rhs), to_digits(prec));
				else core::fpz_set(value, get_value(rhs), to_digits(prec));
			}
		} else {
			set_i64_nothrow(eval_any(rhs, nullptr));
		}

		return *this;
	}

	// Named constructor (there is no other way to construct from a "bare" int64_t)
	static constexpr ztype
	from_i64(int64_t v) noexcept(detail::is_large(prec))
	{
		return (ztype{}).set_i64(v);
	}

	// Named assignment 
	constexpr ztype&
	set_i64(int64_t v) noexcept(detail::is_large(prec))
	{
		if constexpr (!detail::is_large(prec)) {
			if (core::abs_i64(v) > max_abs_value) throw core::integer_overflow{};
		}
		return set_i64_nothrow(v);
	}

	constexpr int64_t
	get_i64() const noexcept(!detail::is_large(prec))
	{
		if constexpr (detail::is_large(prec)) {
			if (std::is_constant_evaluated()) return core::get_i64(value);
			else return core::fpz_get_i64(value);
		} else {
			return value;
		}
	}

	constexpr 
	ztype(const char *str)
	{
		(void)operator=(str);
	}

	constexpr ztype&
	operator=(const char *str)
	{
		using detail::is_large;

		if constexpr (is_large(prec)) {
			if (std::is_constant_evaluated()) core::set_str(value, str, prec);
			else core::fpz_set_str(value, str, prec);
		} else {
			if (std::is_constant_evaluated()) value = core::i64_set_str(str, prec);
			else value = core::fpz_i64_set_str(str, prec);
		}

		return *this;
	}

	constexpr
	ztype(std::initializer_list<core::digit> digits) 
	{ 
		(void)operator=(digits);
	}

	constexpr ztype&
	operator=(std::initializer_list<core::digit> digits)
	{
		using detail::is_large;

		core::digcnt num_digits{ static_cast<uint32_t>(digits.size()) };

		if constexpr (is_large(prec)) {
			if (std::is_constant_evaluated()) {
				core::set_digits(value, digits.begin(), num_digits, prec);
			} else {
				core::fpz_set_digits(value, digits.begin(), num_digits, prec);
			}
		} else {
			if (std::is_constant_evaluated()) {
				value = core::i64_set_digits(digits.begin(), num_digits, prec);
			} else {
				value = core::fpz_i64_set_digits(digits.begin(), num_digits, prec);
			}
		}

		return *this;
	}

	constexpr ztype&
	negate() noexcept 
	{ 
		if constexpr (detail::is_large(prec)) core::negate(value);
		else value = -value;
		return *this;
	}

	constexpr int 
	signum() const noexcept 
	{
		if constexpr (detail::is_large(prec)) return core::signum(value);
		else return value > 0 ? 1 : value == 0 ? 0 : -1;
	}

	static constexpr ztype
	make_max() noexcept
	{
		return (ztype{}).set_max();
	}

	constexpr ztype&
	set_max() noexcept 
	{
		if constexpr (detail::is_large(prec)) {
			if (std::is_constant_evaluated()) core::set_max(value, prec);
			else core::fpz_set_max(value, prec);
		} else {
			value = max_abs_value;
		}
		return *this;
	}

	static constexpr ztype
	make_min() noexcept
	{
		return (ztype{}).set_min();
	}

	constexpr ztype&
	set_min() noexcept 
	{ 
		set_max();
		return negate();
	}

	static ztype
	make_random() noexcept
	{
		return (ztype{}).set_random();
	}

	ztype&
	set_random() noexcept 
	{
		if constexpr (detail::is_large(prec)) core::fpz_random(value, prec);
		else value = core::fpz_random_i64(prec);
		return *this;
	}

	explicit 
	operator long double() const noexcept 
	{ 
		if constexpr (detail::is_large(prec)) return core::fpz_get_ld(value);
		else return static_cast<long double>(value);
	}

	void 
	print_hex(FILE *stream = stdout) const noexcept 
	{ 
		if constexpr (detail::is_large(prec)) core::fpz_fprint_hex(stream, value);
		else {
			if (value < 0) std::fputc('-', stream);
			std::fprintf(stream, "%" PRIX64, core::abs_i64(value));
		}
	}

	void 
	print(FILE *stream = stdout) const noexcept 
	{ 
		if constexpr (detail::is_large(prec)) core::fpz_fprint_dec(stream, value);
		else std::fprintf(stream, "%" PRIi64, value);
	}

	void 
	println_hex(FILE *stream = stdout) const noexcept
	{
		print_hex(stream);
		std::fputc('\n', stream);
	}

	void 
	println(FILE *stream = stdout) const noexcept
	{
		print(stream);
		std::fputc('\n', stream);
	}

	constexpr core::bitcnt 
	bit_size() const noexcept
	{
		if constexpr (detail::is_large(prec)) {
			if (std::is_constant_evaluated()) return core::bit_size(value);
			else return core::fpz_bit_size(value);
		} else return core::bit_size(value);
	}

	template <size_t count, detail::Expr E = const_i64_type<0>, typename Termtype = detail::empty> FPZ_FORCE_INLINE auto
	sum(const E& init_expr = const_i64_type<0>{}, Termtype&& term = std::move(detail::empty{})) noexcept
	{
		return accumulate<count, detail::op_add>(init_expr, term);
	}

	template <size_t count, detail::Expr E = const_i64_type<1>, typename Termtype = detail::empty> FPZ_FORCE_INLINE auto
	product(const E& init_expr = const_i64_type<1>{}, Termtype&& term = std::move(detail::empty{})) noexcept
	{
		return accumulate<count, detail::op_mul>(init_expr, term);
	}

	template <size_t count, detail::Expr E = const_i64_type<0>, typename Termtype = detail::empty> FPZ_FORCE_INLINE auto
	difference(const E& init_expr = const_i64_type<0>{}, Termtype&& term = std::move(detail::empty{})) noexcept
	{
		return accumulate<count, detail::op_sub>(init_expr, term);
	}

	template <size_t count, typename Op, uint32_t dst_prec, uint32_t term_prec>
	friend class detail::accumulate_iter;

	template <core::size_type prec2, size_t count, typename Op, detail::Expr E, typename Termtype> 
	friend class detail::accumulate_range;

	template <core::size_type prec1, uint32_t prec2, detail::SimpleExpr E1, detail::SimpleExpr E2>
	friend constexpr void divrem(ztype<prec1>& quot, ztype<prec2>& rem, const E1& num, const E2& denom) noexcept;

	template <detail::SimpleExpr E> 
	friend constexpr const auto detail::get_value(const E&) noexcept;
};

template <detail::SimpleExpr E1, detail::SimpleExpr E2> FPZ_FORCE_INLINE constexpr auto
operator<=>(const E1& lhs, const E2& rhs) 
{
	using core::bitcnt;
	using detail::is_large, detail::expr_traits, detail::get_value;

	if constexpr (!is_large(expr_traits<E1>::bit_bound) && !is_large(expr_traits<E2>::bit_bound)) {
		return get_value(lhs) <=> get_value(rhs);
	} else {
		if (std::is_constant_evaluated()) {
			return core::cmp(get_value(lhs), get_value(rhs)) <=> 0;
		} else {
			return core::fpz_cmp(get_value(lhs), get_value(rhs)) <=> 0;
		}
	}
}

template <detail::SimpleExpr E1, detail::SimpleExpr E2> FPZ_FORCE_INLINE constexpr bool
operator==(const E1& lhs, const E2& rhs) 
{
	return (lhs <=> rhs) == 0;
}

template <detail::SimpleExpr E> FPZ_FORCE_INLINE constexpr auto
operator<=>(const E& lhs, int64_t rhs)
{
	using core::bitcnt;
	using detail::is_large, detail::expr_traits, detail::get_value;

	if constexpr (!is_large(expr_traits<E>::bit_bound)) {
		return get_value(lhs) <=> rhs;
	} else {
		if (std::is_constant_evaluated()) {
			return core::cmp(get_value(lhs), rhs) <=> 0;
		} else {
			return core::fpz_cmp(get_value(lhs), rhs) <=> 0;
		}
	}
}

template <detail::SimpleExpr E> FPZ_FORCE_INLINE constexpr bool
operator==(const E& lhs, int64_t rhs)
{
	return (lhs <=> rhs) == 0;
}

// Compare absolute values

template <detail::SimpleExpr E1, detail::SimpleExpr E2> FPZ_FORCE_INLINE int
cmpabs(const E1& lhs, const E2& rhs)
{
	using detail::get_value;

	if (std::is_constant_evaluated()) {
		return core::cmpabs(get_value(lhs), get_value(rhs));
	} else {
		return core::fpz_cmpabs(get_value(lhs), get_value(rhs));
	}
}

template <detail::SimpleExpr E> FPZ_FORCE_INLINE int
cmpabs(const E& lhs, int64_t rhs)
{
	using detail::get_value;

	if (std::is_constant_evaluated()) {
		return core::cmpabs(get_value(lhs), rhs);
	} else {
		return core::fpz_cmpabs(get_value(lhs), rhs);
	}
}

template <core::size_type quot_prec, uint32_t rem_prec, detail::SimpleExpr E1, detail::SimpleExpr E2>
constexpr void 
divrem(ztype<quot_prec>& quot, ztype<rem_prec>& rem, const E1& num, const E2& denom) noexcept
{
	using core::bitcnt, core::to_digits;
	using detail::is_large, detail::expr_traits, detail::get_value;

	using E1Tr = expr_traits<E1>;
	using E2Tr = expr_traits<E2>;

	static_assert(bitcnt{quot_prec} >= E1Tr::bit_bound, "quotient precision is too small");
	static_assert(bitcnt{rem_prec} >= std::min(E1Tr::bit_bound, E2Tr::bit_bound), "remainder precision is too small");

	if constexpr (is_large(E1Tr::bit_bound) && is_large(E2Tr::bit_bound)) {
		if (std::is_constant_evaluated()) {
			core::divrem(quot.value, rem.value, get_value(num), get_value(denom), 
			    to_digits(quot.prec), to_digits(rem.prec));
		} else {
			core::fpz_divrem(quot.value, rem.value, get_value(num), get_value(denom), 
			    to_digits(quot.prec), to_digits(rem.prec));
		}
	} else if constexpr (is_large(E1Tr::bit_bound)) {
		if (std::is_constant_evaluated()) {
			rem.set_int64_nothrow(core::divrem_i64(quot.value, get_value(num), 
			    get_value(denom), to_digits(quot.prec)));
		} else {
			rem.set_int64_nothrow(core::fpz_divrem_i64(quot.value, get_value(num), 
			    get_value(denom), to_digits(quot.prec)));
		}
	} else if constexpr (is_large(E2Tr::bit_bound)) {
		int64_t q{0}, r{0};
		if (std::is_constant_evaluated()) r = core::i64_divrem(&q, get_value(num), get_value(denom));
		else r = core::fpz_i64_divrem(&q, get_value(num), get_value(denom));
		quot.set_i64_nothrow(q);
		rem.set_i64_nothrow(r);
	} else {
		int64_t n = get_value(num);
		int64_t d = get_value(denom);

		int64_t q = n / d;
		rem.set_int64_nothrow(n - q*d);
		quot.set_int64_nothrow(q);
	}

	FPZ_ASSERT_VALID(get_value(quot), quot.prec);
	FPZ_ASSERT_VALID(get_value(rem), rem.prec);
}

template <detail::SimpleExpr E1, detail::SimpleExpr E2> FPZ_FORCE_INLINE long double 
fdiv(const E1& num, const E2& denom)
{
	using core::bitcnt;
	using detail::is_large, detail::expr_traits, detail::get_value;
	using E1Tr = expr_traits<E1>;
	using E2Tr = expr_traits<E2>;

	if constexpr (is_large(E1Tr::bit_bound) && is_large(E2Tr::bit_bound)) 
		return core::fpz_fdiv(get_value(num), get_value(denom));
	else if constexpr (is_large(E1Tr::bit_bound)) {
		auto tmp = ztype<64>::from_i64(get_value(denom));
		return core::fpz_fdiv(get_value(num), get_value(tmp));
	} else if constexpr (is_large(E2Tr::bit_bound)) {
		auto tmp = ztype<64>::from_i64(get_value(num));
		return core::fpz_fdiv(get_value(tmp), get_value(denom));
	} else {
		return static_cast<long double>(get_value(num)) / 
			static_cast<long double>(get_value(denom));
	}
}

template <detail::Expr E> constexpr auto
ztype_from_expr(const E& e) noexcept
{
	return ztype<detail::expr_traits<E>::bit_bound.nb>{e};
}

// Constant from string

template <typename S>
struct const_str_type {
	static constexpr ztype<detail::estimate_bit_size(S::str()).nb> z{S::str()};
};

// Unsigned constant from list of digits (specified from most to least significant)

template <core::digit... ds>
struct const_digits_type {
	static constexpr ztype<detail::bit_size({ ds... }).nb> z{{ ds... }};
};

namespace detail {

// Separate out expression-related properties into their own traits classes,
// rather than cluttering up the integer types with them (this also gets us traits
// defined for all types, not just the stuff we define here).

template <typename T>
struct expr_traits_nocvref {
	static constexpr bool is_expr = false;
	static constexpr bool is_compound = false;
	static constexpr bool has_storage = false;
};

template <core::size_type prec>
struct expr_traits_nocvref<ztype<prec>> {
	static constexpr bool 	is_expr = true;
	static constexpr bool 	is_compound = false;
	static constexpr bool 	has_storage = true;

	static constexpr bitcnt bit_bound{prec};
	static constexpr expr_bound bound{ bit_bound, 1 };

	template <typename D> static constexpr bool subexpr_may_alias_with = 
	    bit_bound == D::prec && std::is_same_v<D, local_var<D::prec.nb>>;
};

// Common stuff for all compound expressions 

struct expr_traits_compound_base {
	static constexpr bool is_expr = true;
	static constexpr bool is_compound = true;
	static constexpr bool has_storage = false;
};

template <typename E>
struct expr_traits_nocvref<negate_expr<E>> : public expr_traits_compound_base {
	static constexpr expr_bound bound = expr_traits<E>::bound;
	static constexpr bitcnt bit_bound = expr_traits<E>::bit_bound;

	template <typename D> static constexpr bool subexpr_may_alias_with = 
		expr_traits<E>::template subexpr_may_alias_with<D>;
};

template <typename E>
struct expr_traits_nocvref<abs_expr<E>> : 
	public expr_traits_nocvref<negate_expr<E>> {};

template <typename E1, typename E2, typename Op>
struct expr_traits_nocvref<binary_expr<E1, E2, Op>> : public expr_traits_compound_base {
	static constexpr expr_bound bound = Op::bound(expr_traits<E1>::bound, expr_traits<E2>::bound);
	static constexpr bitcnt bit_bound = static_cast<bitcnt>(bound);

	template <typename D> static constexpr bool subexpr_may_alias_with = 
	    expr_traits<E1>::template subexpr_may_alias_with<D> ||
	    expr_traits<E2>::template subexpr_may_alias_with<D>;
};

template <typename E, typename Op, uint32_t nbits>
struct expr_traits_nocvref<shift_expr<E, Op, nbits>> : public expr_traits_compound_base {
	
	static constexpr expr_bound bound = Op::bound(expr_traits<E>::bound, bitcnt{nbits});
	static constexpr bitcnt bit_bound = static_cast<bitcnt>(bound);

	template <typename D> static constexpr bool subexpr_may_alias_with = 
	    expr_traits<E>::template subexpr_may_alias_with<D>;
};

template <core::size_type bit_size>
struct expr_traits_const_base {
	static constexpr bool is_expr = true;
	static constexpr bool is_compound = false;
	static constexpr bool has_storage = false;

	static constexpr bitcnt bit_bound { bit_size };
	static constexpr expr_bound bound { .exp = bit_bound, .m = 1 };

	template <typename D> static constexpr bool subexpr_may_alias_with = false;
};

template <int64_t v>
struct expr_traits_nocvref<const_i64_type<v>> : 
	public expr_traits_const_base<core::bit_size(v).nb> {};

template <digit... ds>
struct expr_traits_nocvref<const_digits_type<ds...>>: 
	public expr_traits_const_base<const_digits_type<ds...>::z.bit_size().nb> {};

template <typename S>
struct expr_traits_nocvref<const_str_type<S>>: 
	public expr_traits_const_base<const_str_type<S>::z.bit_size().nb> {};

} // namespace detail

template <int64_t v>
constexpr auto const_i64 = const_i64_type<v>{};

template <core::digit... ds>
constexpr auto const_digits = const_digits_type<ds...>{};

// While we wait for (more widely available) support for CNTTP ...

#define FPZ_CONST_STR(s)	[]() {							\
	struct strprovider { static constexpr const char *str() { return s; } };	\
	return fpz::const_str_type<strprovider>{};					\
}()

// Operators that produce expressions. Compound expressions are dealt with as temporary values 
// (and the structures representing them are expected to be optimized out), while 
// ztype variables are passed by const reference. 

template <detail::Expr T1, detail::Expr T2> constexpr auto
operator+(const T1& e1, const T2& e2)
{
	return detail::binary_expr<T1, T2, detail::op_add> { e1, e2 };
}

template <detail::Expr T1, detail::Expr T2> constexpr auto
operator-(const T1& e1, const T2& e2)
{
	return detail::binary_expr<T1, T2, detail::op_sub> { e1, e2 };
}

template <detail::Expr T1, detail::Expr T2> constexpr auto
operator*(const T1& e1, const T2& e2)
{
	return detail::binary_expr<T1, T2, detail::op_mul> { e1, e2 };
}

template <detail::Expr T1, detail::Expr T2> constexpr auto
operator/(const T1& e1, const T2& e2)
{
	return detail::binary_expr<T1, T2, detail::op_div> { e1, e2 };
}

template <detail::Expr T1, detail::Expr T2> constexpr auto
operator%(const T1& e1, const T2& e2)
{
	return detail::binary_expr<T1, T2, detail::op_mod> { e1, e2 };
}

template <detail::Expr T> constexpr auto
operator-(const T& e)
{
	return detail::negate_expr<T> { e };
}

template <detail::Expr T, int64_t nbits> constexpr auto
operator<<(const T& e, const const_i64_type<nbits>&)
{
	static_assert(nbits >= 0 && nbits <= std::numeric_limits<uint32_t>::max(), "bit count out of range");

	if constexpr (nbits == 0) return e;
	else return detail::shift_expr<T, detail::op_shift_left, nbits> { e };
}

template <detail::Expr T, int64_t nbits> constexpr auto
operator>>(const T& e, const const_i64_type<nbits>&)
{
	static_assert(nbits >= 0 && nbits <= std::numeric_limits<uint32_t>::max(), "bit count out of range");

	if constexpr (nbits == 0) return e;
	else return detail::shift_expr<T, detail::op_shift_right, nbits> { e };
}

template <detail::Expr T, int64_t nbits> constexpr auto
mod2exp(const T& e, const const_i64_type<nbits>&)
{
	if constexpr (nbits >= detail::expr_traits<T>::prec) return e;
	else return detail::shift_expr<T, detail::op_mod_2exp, nbits> { e };
}

template <detail::Expr T> constexpr auto
abs(const T& e)
{
	return detail::abs_expr<T> { e };
}

} // namespace fpz

using fpz::operator""_bits;

#endif	// FPZ_H
