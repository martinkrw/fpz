#ifndef FPZ_CORE_H
#define FPZ_CORE_H

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

#include <bit>
#include <cassert>
#include <cinttypes>
#include <compare>
#include <concepts>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

// Header for the low-level "core" functions of fpz, written
// in C++ (20) and assembly. Deals with fixed precision integers represented
// as an array of (64-bit) words. 

namespace fpz::core {

// Force the compiler to do the right thing in cases where we think we know better

#define FPZ_FORCE_INLINE	__attribute__((always_inline)) inline

// Make branches that are never supposed to be taken unless the library or 
// program has bugs (such as the ones that check the overflow flag on a 
// variable we want to use and then aborts the program if that is set)
// correctly predicted by static prediction (should be replaced by
// the attributes [[likely]] and [[unlikely]] eventually).

#define FPZ_UNLIKELY(cond)	__builtin_expect(cond, false)
#define FPZ_LIKELY(cond)	__builtin_expect(cond, true)

// Some typedefs (only intended as annotations)

using size_type = uint32_t; // 4 billion bits should be enough for everyone ...

using word = uint64_t;
using digit = uint64_t;
using dstptr = word *;
using srcptr = const word *;



// Sign values

constexpr uint8_t POSITIVE = 0;
constexpr uint8_t NEGATIVE = 1;

struct u128 {
	digit low;
	digit high;
};

// Exception structures

struct error : public std::exception { };

struct integer_overflow : public error { 
	virtual const char *
	what() const noexcept 
	{ 
		return "integer overflow";
	}
};

struct invalid_argument : public error { 
	virtual const char *
	what() const noexcept 
	{
		return "invalid argument";
	}
};

// Split a two's complement signed integer (the parameter i) into sign 
// (0 if POSITIVE, 1 if NEGATIVE, returned in *sign) and magnitude (return value) 
// without branching.

constexpr uint64_t
sign_magnitude(int64_t i, uint8_t *sign) noexcept
{
	uint64_t s = -(i < 0);
	if (sign != nullptr) *sign = s & 1;
	return (s ^ static_cast<uint64_t>(i)) - s;
}

// Compute absolute value (magnitude) of an int64_t as a uint64_t (-INT64_MIN is 
// outside the range of an int64_t, so we can't return a value of the same type)

constexpr uint64_t
abs_i64(int64_t i) noexcept
{
	return sign_magnitude(i, nullptr);
}

// Shorthand

constexpr digit MAX_DIGIT = std::numeric_limits<digit>::max();
constexpr uint64_t I64_MAXABS_NEG = abs_i64(std::numeric_limits<int64_t>::min());
constexpr uint64_t I64_MAXABS_POS = abs_i64(std::numeric_limits<int64_t>::max());

// Structures representing bit and digit counts, respectively. We make them
// incompatible types, rather than simple typedefs, so that they can't be
// used interchangeably by mistake.

struct bitcnt {
	size_type nb;
	//explicit constexpr operator bool() { return nb != 0; }
	constexpr auto operator<=>(const bitcnt&) const = default;
	constexpr bitcnt operator+=(bitcnt b) { nb += b.nb; return *this; }
	constexpr bitcnt operator-=(bitcnt b) { nb -= b.nb; return *this; }

	constexpr bitcnt operator++() { ++nb; return *this; }
	constexpr bitcnt operator--() { --nb; return *this; }
};

// Digit/word size in bits. Note that this can't be changed without comprehensive
// changes to the code.

constexpr bitcnt BITS_PER_DIGIT{64};

struct digcnt {
	size_type nd;
	//explicit constexpr operator bool() { return nd != 0; }
	constexpr auto operator<=>(const digcnt&) const = default;
	constexpr digcnt operator+=(digcnt d) { nd += d.nd; return *this; }
	constexpr digcnt operator-=(digcnt d) { nd -= d.nd; return *this; }

	constexpr digcnt operator++() { ++nd; return *this; }
	constexpr digcnt operator--() { --nd; return *this; }
};

constexpr bitcnt operator+(bitcnt b1, bitcnt b2) { return bitcnt{ b1.nb + b2.nb }; }
constexpr bitcnt operator-(bitcnt b1, bitcnt b2) { return bitcnt{ b1.nb - b2.nb }; }

constexpr digcnt operator+(digcnt d1, digcnt d2) { return digcnt{ d1.nd + d2.nd }; }
constexpr digcnt operator-(digcnt d1, digcnt d2) { return digcnt{ d1.nd - d2.nd }; }

constexpr digcnt 
to_digits(bitcnt b) 
{
	return digcnt{ (b.nb + BITS_PER_DIGIT.nb - 1)/ BITS_PER_DIGIT.nb }; 
}

constexpr digcnt
whole_digits(bitcnt b)
{
	return digcnt{ b.nb / BITS_PER_DIGIT.nb };
}

constexpr bitcnt
to_bits(digcnt d)
{
	return bitcnt{ d.nd * BITS_PER_DIGIT.nb };	
}

constexpr bitcnt
remainder_bits(bitcnt b)
{
	return bitcnt{ b.nb % BITS_PER_DIGIT.nb };
}

FPZ_FORCE_INLINE constexpr bitcnt
bit_size(uint64_t x)
{
	return bitcnt{ static_cast<size_type>(std::bit_width(x)) };
}

FPZ_FORCE_INLINE constexpr bitcnt
bit_size_i64(int64_t x)
{
	return bit_size(abs_i64(x));
}

FPZ_FORCE_INLINE constexpr digcnt
max(digcnt i, digcnt j)
{
	return i >= j ? i : j; 
}

FPZ_FORCE_INLINE constexpr digcnt
min(digcnt i, digcnt j)
{
	return i > j ? j : i;
}


// Separately compiled library functions (either assembly or C++).
// Use a prefix to prevent library names from conflicting with names in 
// user code (the functions written in C++ don't strictly need this,
// but we do it anyway for consistency)

extern "C" {

bool 		fpz_is_valid(srcptr, bitcnt) noexcept;
int		fpz_signum(srcptr) noexcept;

void		fpz_set_max(dstptr, bitcnt) noexcept;

void		fpz_set(dstptr, srcptr, digcnt) noexcept;

void		fpz_set_digits(dstptr, const digit *, digcnt, bitcnt);
void		fpz_set_str(dstptr, const char *, bitcnt);
int64_t		fpz_i64_set_str(const char *, bitcnt);
int64_t		fpz_i64_set_digits(const digit *, digcnt, bitcnt);

int64_t		fpz_get_i64(srcptr);

long double	fpz_get_ld(srcptr) noexcept;
digit		fpz_get_bit(srcptr, bitcnt) noexcept;

void		fpz_random(dstptr, bitcnt) noexcept;
int64_t		fpz_random_i64(bitcnt) noexcept;
void		fpz_random_sparse(dstptr, bitcnt) noexcept;

bitcnt		fpz_bit_size(srcptr op) noexcept;

int		fpz_fprint_hex(FILE *, srcptr) noexcept;
int		fpz_fprint_dec(FILE *, srcptr) noexcept;

void		fpz_fprint_debug(FILE*, srcptr, bitcnt) noexcept;

void		fpz_add_i64_i64(dstptr, int64_t, int64_t, digcnt) noexcept;
void		fpz_sub_i64_i64(dstptr, int64_t, int64_t) noexcept;

void	 	fpz_add_i64(dstptr, srcptr, int64_t, digcnt) noexcept;
void		fpz_sub_i64(dstptr, srcptr, int64_t, digcnt) noexcept;
void		fpz_i64_sub(dstptr, int64_t, srcptr, digcnt) noexcept;

void		fpz_add(dstptr, srcptr, srcptr, digcnt) noexcept;
void		fpz_sub(dstptr, srcptr, srcptr, digcnt) noexcept;

void		fpz_mul_i64_i64(dstptr, int64_t, int64_t, digcnt) noexcept;
void		fpz_mul_i64(dstptr, srcptr, int64_t, digcnt) noexcept;

u128		fpz_addmul_u64_u64(digit, digit, digit, digit) noexcept;
digit		fpz_addmul_u64_self(dstptr, digit, digit, digcnt) noexcept;
void		fpz_mul(dstptr, srcptr, srcptr, digcnt) noexcept;

int64_t		fpz_i64_divrem(int64_t*, int64_t, srcptr) noexcept;

int64_t		fpz_divrem_i64(dstptr, srcptr, int64_t, digcnt) noexcept;
void		fpz_divrem(dstptr, dstptr, srcptr, srcptr, digcnt, digcnt) noexcept;
int		fpz_divrem10(dstptr, srcptr, digcnt) noexcept;

long double	fpz_fdiv(srcptr num, srcptr denom) noexcept;

void 		fpz_shl(dstptr, srcptr, bitcnt, bitcnt) noexcept;
void 		fpz_shl_i64(dstptr, int64_t, bitcnt, bitcnt) noexcept;
void 		fpz_shr(dstptr, srcptr, bitcnt, bitcnt) noexcept;

void 		fpz_truncate(dstptr, srcptr, bitcnt) noexcept;

void		fpz_assertion_failure(const char *, const char *, int) noexcept;
void 		fpz_use_after_overflow_handler(void) noexcept;

}	// extern "C"

int		fpz_cmp(srcptr op1, srcptr op2) noexcept;
int		fpz_cmp(srcptr op1, int64_t op2) noexcept;
int		fpz_cmpabs(srcptr op1, srcptr op2) noexcept;
int		fpz_cmpabs(srcptr, int64_t) noexcept;
int		fpz_cmpabs(int64_t op1, int64_t op2) noexcept;

// Integer header containing sign and size, which is the number of digits
// used to represent the current value. An fpz integer is simply an array
// of 64-bit words, the first of which is the header and the remaining ones
// are digits, stored from least to most significant. 

typedef struct {
	uint8_t 	sign;
	uint8_t 	overflow_flag;
	uint16_t	pad;
	digcnt 		size;
} header;

// Just in case something is very wrong ... 
// (if any of these fail, our assembly code needs to change) 

static_assert(sizeof(header) == sizeof(word), "header has unexpected size");
static_assert(offsetof(header, overflow_flag) == 1 && 
	offsetof(header, size) == 4, "header has unexpected layout");

// Number of bytes to allocate for a for an integer variable with space
// for num_digits digits (used only with alloca, as things get less verbose with new)

constexpr size_t 
int_alloc_size(digcnt num_digits)
{
	return (static_cast<size_t>(num_digits.nd) + 1)*sizeof(word);
}

// Read and write headers.

// Had it not been for the requirement that this be constexpr, 
// this would be a simple memcpy of 64 bits from src[0] to the return
// value. So, while we wait for std::bit_cast, this is it.

constexpr header
get_header(srcptr src) noexcept
{
	// Will be replaced by std::bit_cast eventually
	return __builtin_bit_cast(header, *src);
}

constexpr void
set_header(dstptr dst, header hdr) noexcept
{
	// Will be replaced by std::bit_cast eventually
	dst[0] = __builtin_bit_cast(word, hdr);
}

constexpr digcnt
get_size(srcptr src) noexcept
{
	return digcnt{ static_cast<size_type>(src[0] >> 32) };
}

constexpr void 
set_size(dstptr dst, digcnt size) noexcept
{
	dst[0] = (word{size.nd} << 32) | static_cast<size_type>(dst[0]);
}

constexpr uint8_t
get_sign(srcptr src) noexcept
{
	return static_cast<uint8_t>(src[0]);
}

constexpr void
set_sign(dstptr dst, uint8_t sign) noexcept
{
	dst[0] &= 0xFFFFFFFFFFFFFF00;
	dst[0] |= sign;
}

constexpr void
negate(dstptr dst) noexcept
{
	dst[0] ^= word{1};
}

// Signum function (not to be confused with sign)
constexpr int
signum(srcptr op) noexcept
{
	header hdr = get_header(op);
	if (hdr.size != digcnt{0}) return hdr.sign == POSITIVE ? 1 : -1;
	else return 0;
}

// This is intentionally not marked constexpr, even though it could be. This
// ensures that compile-time computations simply stop on overflow. 

FPZ_FORCE_INLINE void
set_overflow_flag(dstptr dst) noexcept
{
	header hdr = { .sign = POSITIVE, .overflow_flag = 1, .size = digcnt{0} };
	set_header(dst, hdr);
}

// Check if the overflow flag is set and take appropriate action.

constexpr void
check_overflow_flag(header hdr) noexcept
{
	if (FPZ_UNLIKELY(hdr.overflow_flag)) 
		fpz_use_after_overflow_handler();
}

constexpr void
check_overflow_flag_2(header hdr1, header hdr2) noexcept
{
	if (FPZ_UNLIKELY(hdr1.overflow_flag | hdr2.overflow_flag)) 
		fpz_use_after_overflow_handler();
}

template <typename T> 
requires std::same_as<T, dstptr> || std::same_as<T, srcptr> constexpr T
digits(T ptr) noexcept
{
	return ptr + 1;
}

FPZ_FORCE_INLINE constexpr bitcnt
bit_size(srcptr op) noexcept
{
	header hdr = get_header(op);
	check_overflow_flag(hdr);

	if (hdr.size == digcnt{0}) return bitcnt{0};
	else return to_bits(hdr.size - digcnt{1}) + 
		bit_size(digits(op)[hdr.size.nd - 1]);
}

// Check if op contains a valid integer of precision prec.

FPZ_FORCE_INLINE constexpr bool
is_valid(srcptr op, bitcnt prec) noexcept
{
	header hdr = get_header(op);
	digcnt digprec = to_digits(prec);

	if (hdr.overflow_flag) return false;

	if (hdr.sign != POSITIVE && hdr.sign != NEGATIVE) return false;

	if (hdr.size == digcnt{0}) return digits(op)[0] == 0;

	// bit_width() must access the size'th digit, which should 
	// be avoided if size > digprec. 
	if (hdr.size == digprec) return prec >= bit_size(op);
	return hdr.size < digprec;
}

// Set the value of an fpz integer to that of a built-in signed integer 

FPZ_FORCE_INLINE constexpr void
set_i64(dstptr dst, int64_t value) noexcept
{
	header hdr{ .sign = value < 0 ? NEGATIVE : POSITIVE, .size = digcnt{value != 0} };
	digits(dst)[0] = sign_magnitude(value, &hdr.sign);
	set_header(dst, hdr);
}

// Set the value of dst to the value of src. This assumes that dst is 
// large enough to contain the value of src. 

FPZ_FORCE_INLINE constexpr void
set(dstptr dst, srcptr src, digcnt digprec) noexcept
{
	header hdr = get_header(src);
	check_overflow_flag(hdr);

	if (FPZ_UNLIKELY(digprec < hdr.size)) {
		set_overflow_flag(dst);
		return;
	}

	set_header(dst, hdr);

	// This ensures that the first digit is copied even if hdr.size == 0
	digits(dst)[0] = digits(src)[0];

	for (digcnt i{1}; i < hdr.size; ++i) 
		digits(dst)[i.nd] = digits(src)[i.nd];
}

FPZ_FORCE_INLINE constexpr void
set_max(dstptr dst, bitcnt prec) noexcept
{
	digcnt n = whole_digits(prec);
	bitcnt rem = remainder_bits(prec);

	for (digcnt i{0}; i < n; ++i) digits(dst)[i.nd] = MAX_DIGIT;

	if (rem != bitcnt{0}) {
		digits(dst)[n.nd] = (digit{1} << rem.nb) - 1;
		++n;
	}

	header hdr = { .sign = POSITIVE, .size = n };
	set_header(dst, hdr);
}

// Convert an fpz integer to an int64_t, if it is small enough to fit. If so, return it.
// Otherwise, throw an overflow_error.

FPZ_FORCE_INLINE constexpr int64_t
get_i64(srcptr src)
{
	header hdr = get_header(src);
	check_overflow_flag(hdr);

	if (hdr.size > digcnt{1}) throw integer_overflow{};

	digit d = digits(src)[0];
	if (hdr.sign == NEGATIVE && d <= I64_MAXABS_NEG) return -d;
	else if (hdr.sign == POSITIVE && d <= I64_MAXABS_POS) return d;

	throw integer_overflow{};
}

FPZ_FORCE_INLINE constexpr digit
get_bit(srcptr src, bitcnt index) noexcept
{
	header hdr = get_header(src);
	check_overflow_flag(hdr);

	if (index >= to_bits(hdr.size)) return 0;

	digit d = digits(src)[whole_digits(index).nd];
	bitcnt rem = remainder_bits(index);
	digit mask = digit{1} << rem.nb;
	return (d & mask) >> rem.nb;
}

FPZ_FORCE_INLINE constexpr void
set_bit(dstptr dst, bitcnt index, bitcnt prec)
{
	if (FPZ_UNLIKELY(index >= prec)) {
		set_overflow_flag(dst);
		throw integer_overflow{};
	}

	digcnt size = get_size(dst);

	// Check if the size needs to be adjusted
	if (index >= to_bits(size)) {
		digcnt newsize = to_digits(index + bitcnt{1});

		// Zero out any digits in between current size and new size
		for (digcnt i = size - digcnt{1}; i < newsize; ++i) 
			digits(dst)[i.nd] = 0;
		set_size(dst, newsize);
	}

	digit d = digit{1} << remainder_bits(index).nb;
	digits(dst)[whole_digits(index).nd] |= d;
}

FPZ_FORCE_INLINE constexpr void
truncate(dstptr dst, srcptr src, bitcnt prec) noexcept
{
	header hdr = get_header(src);
	check_overflow_flag(hdr);
	
	header dst_hdr = { .sign = hdr.sign, 
		.size = min(to_digits(prec), hdr.size) };

	set_header(dst, dst_hdr);

	for (digcnt i{0}; i + digcnt{1} < dst_hdr.size; ++i)
		digits(dst)[i.nd] = digits(src)[i.nd];

	bitcnt rem = remainder_bits(prec);
	digit mask = (rem == bitcnt{0}) ? MAX_DIGIT : (digit{1} << rem.nb) - 1;
	digits(dst)[dst_hdr.size.nd - 1] = digits(src)[dst_hdr.size.nd - 1] & mask;
}

FPZ_FORCE_INLINE constexpr int
cmp_common(const digit *op1, const digit *op2, digcnt n1, digcnt n2) noexcept
{
	if (n1 != n2) return n1 > n2 ? 1 : -1;

	for (digcnt i = n1 - digcnt{1}; i < n1; --i) {
		if (op1[i.nd] > op2[i.nd]) return 1;
		else if (op1[i.nd] < op2[i.nd]) return -1;
	}
	return 0;
}

FPZ_FORCE_INLINE constexpr int
cmpabs_u64(const digit *op1, uint64_t op2, digcnt n) noexcept
{
	if (n > digcnt{1}) return 1;

	digit d = op1[0];
	return d < op2 ? -1 : d > op2;
}

FPZ_FORCE_INLINE constexpr int
cmp(srcptr op1, srcptr op2) noexcept
{
	header hdr1 = get_header(op1);
	header hdr2 = get_header(op2);

	check_overflow_flag_2(hdr1, hdr2);

	int factor = 1 - 2*static_cast<int>(hdr1.sign);

	if (hdr1.sign == hdr2.sign) return factor*cmp_common(digits(op1), digits(op2), hdr1.size, hdr2.size);
	else return factor*(hdr1.size != digcnt{0} || hdr2.size != digcnt{0});
}

FPZ_FORCE_INLINE constexpr int
cmp(srcptr op1, int64_t op2) noexcept
{
	uint8_t sign_op2 = POSITIVE;
	digit abs_op2 = sign_magnitude(op2, &sign_op2);

	header hdr = get_header(op1);
	check_overflow_flag(hdr);

	int factor = 1 - 2*static_cast<int>(hdr.sign);

	if (hdr.sign == sign_op2) 
		return factor*cmpabs_u64(digits(op1), abs_op2, hdr.size);
	else return factor*(hdr.size != digcnt{0} || op2 != 0);
}

FPZ_FORCE_INLINE constexpr int
cmp(int64_t op1, srcptr op2) noexcept
{
	return -cmp(op2, op1);
}

FPZ_FORCE_INLINE constexpr int
cmpabs(srcptr op1, srcptr op2) noexcept
{
	header hdr1 = get_header(op1);
	header hdr2 = get_header(op2);

	check_overflow_flag_2(hdr1, hdr2);

	if (hdr1.size != hdr2.size) 
		return hdr1.size < hdr2.size ? -1 : hdr1.size > hdr2.size;
	else return cmp_common(digits(op1), digits(op2), hdr1.size, hdr2.size);
}

FPZ_FORCE_INLINE constexpr int
cmpabs(srcptr op1, int64_t op2) noexcept
{
	digit abs_op2 = abs_i64(op2);
	header hdr = get_header(op1);

	check_overflow_flag(hdr);

	return cmpabs_u64(digits(op1), abs_op2, hdr.size);
}

FPZ_FORCE_INLINE constexpr int
cmpabs(int64_t op1, int64_t op2) noexcept
{
	digit abs_op1 = abs_i64(op1);
	digit abs_op2 = abs_i64(op2);

	return abs_op1 < abs_op2 ? -1 : abs_op1 > abs_op2;
}

// Addition/subtraction


// Add with carry

constexpr u128
adc(digit term1, digit term2, digit carry_in) noexcept
{
	digit sum = term1 + term2 + carry_in;
	return u128{ .low = sum, .high = (term1 + term2 < term1 || sum < carry_in) }; 
}

// Subtract with borrow

constexpr u128
sbb(digit term1, digit term2, digit borrow_in) noexcept
{
	digit diff = term1 - term2 - borrow_in;
	return u128{ .low = diff, .high = (term1 < term2 || term1 - term2 < borrow_in) };
}

constexpr digit
uadd(digit *dst, const digit *op1, const digit *op2, digcnt n1, digcnt n2) noexcept
{
	u128 r = { 0, 0 };

	for (digcnt i{0}; i < n1 || i < n2; ++i) {
		digit d1 = i < n1 ? op1[i.nd] : 0;
		digit d2 = i < n2 ? op2[i.nd] : 0;

		r = adc(d1, d2, r.high);

		dst[i.nd] = r.low;
	}
	
	return r.high;
}

constexpr void
usub(digit *dst, const digit *op1, const digit *op2, digcnt n1, digcnt n2) noexcept
{
	u128 r = { 0, 0 };

	for (digcnt i{0}; i < n1 || i < n2; ++i) {
		digit d1 = i < n1 ? op1[i.nd] : 0;
		digit d2 = i < n2 ? op2[i.nd] : 0;

		r = sbb(d1, d2, r.high);
		dst[i.nd] = r.low;
	}
}

constexpr void
add(dstptr dst, srcptr op1, srcptr op2, digcnt dst_cap) noexcept
{
	header hdr1 = get_header(op1);
	header hdr2 = get_header(op2);

	check_overflow_flag_2(hdr1, hdr2);

	header hdr = { .sign = hdr1.sign, .size = max(hdr1.size, hdr2.size) };

	if (dst_cap < hdr.size) {
		set_overflow_flag(dst);
		return;
	}

	digit carry = 0;

	if (hdr1.sign == hdr2.sign) {
		carry = uadd(digits(dst), digits(op1), digits(op2), hdr1.size, hdr2.size);

		if (hdr.size + digcnt{1} <= dst_cap) {
			digits(dst)[hdr.size.nd] = carry;
			hdr.size += digcnt{carry != 0};
		} else if (carry != 0) {
			set_overflow_flag(dst);
			return;
		}

	} else {
		if (cmpabs(op1, op2) >= 0) {
			usub(digits(dst), digits(op1), digits(op2), hdr1.size, hdr2.size);
		} else {
			hdr.sign ^= 1;
			usub(digits(dst), digits(op2), digits(op1), hdr2.size, hdr1.size);
		}

		while (hdr.size != digcnt{0} && digits(dst)[hdr.size.nd - 1] == 0) --hdr.size;
	}

	set_header(dst, hdr);
}

constexpr void
sub(dstptr dst, srcptr op1, srcptr op2, digcnt dst_cap) noexcept
{
	header hdr1 = get_header(op1);
	header hdr2 = get_header(op2);

	check_overflow_flag_2(hdr1, hdr2);

	header hdr = { .sign = hdr1.sign, .size = max(hdr1.size, hdr2.size) };

	if (dst_cap < hdr.size) {
		set_overflow_flag(dst);
		return;
	}

	digit carry = 0;

	if (hdr1.sign != hdr2.sign) {
		carry = uadd(digits(dst), digits(op1), digits(op2), hdr1.size, hdr2.size);

		if (hdr.size + digcnt{1} <= dst_cap) {
			digits(dst)[hdr.size.nd] = carry;
			hdr.size += digcnt{carry != 0};
		} else if (carry != 0) {
			set_overflow_flag(dst);
			return;
		}

	} else {
		if (cmpabs(op1, op2) >= 0) {
			usub(digits(dst), digits(op1), digits(op2), hdr1.size, hdr2.size);
		} else {
			hdr.sign ^= 1;
			usub(digits(dst), digits(op2), digits(op1), hdr2.size, hdr1.size);
		}

		while (hdr.size != digcnt{0} && digits(dst)[hdr.size.nd - 1] == 0) --hdr.size;
	}
	
	set_header(dst, hdr);
}

constexpr void
add_i64_i64(dstptr dst, int64_t op1, int64_t op2, digcnt digprec) noexcept
{	
	word op1_tmp[2];
	word op2_tmp[2];

	set_i64(op1_tmp, op1);
	set_i64(op2_tmp, op2);

	add(dst, op1_tmp, op2_tmp, digprec);
}

constexpr void
sub_i64_i64(dstptr dst, int64_t op1, int64_t op2) noexcept
{
	word op1_tmp[2];
	word op2_tmp[2];

	set_i64(op1_tmp, op1);
	set_i64(op2_tmp, op2);

	sub(dst, op1_tmp, op2_tmp, digcnt{2});
}

constexpr void
add_i64(dstptr dst, srcptr op1, int64_t op2, digcnt digprec) noexcept
{
	word op2_tmp[2];
	set_i64(op2_tmp, op2);

	return add(dst, op1, op2_tmp, digprec);
}

constexpr void
sub_i64(dstptr dst, srcptr op1, int64_t op2, digcnt digprec) noexcept
{
	word op2_tmp[2];
	set_i64(op2_tmp, op2);

	sub(dst, op1, op2_tmp, digprec);
}

constexpr void
i64_sub(dstptr dst, int64_t op1, srcptr op2, digcnt digprec) noexcept
{
	word op1_tmp[2];
	set_i64(op1_tmp, op1);

	sub(dst, op1_tmp, op2, digprec);
}

// Multiplication

// Compute high:low = f1*f2 + term1 + term2

constexpr u128
addmul_u64_u64(digit f1, digit f2, digit term1, digit term2) noexcept
{
	constexpr word mask = 0x00000000FFFFFFFF;
	digit lh = (f1 & mask)*(f2 >> 32);
	digit hl = (f1 >> 32)*(f2 & mask);
	
	u128 res = adc((f1 & mask)*(f2 & mask), lh << 32, 0);
	res.high += (f1 >> 32) * (f2 >> 32);

	res.low += hl << 32;
	res.high += (res.low < (hl << 32));
	res.high += (lh >> 32) + (hl >> 32);
	
	res.low += term1;
	res.high += res.low < term1;

	res.low += term2;
	res.high += res.low < term2;

	return res;
}

constexpr void
mul(dstptr dst, srcptr op1, srcptr op2, digcnt digprec) noexcept
{
	header hdr1 = get_header(op1);
	header hdr2 = get_header(op2);

	check_overflow_flag_2(hdr1, hdr2);

	header hdr{ .sign = static_cast<uint8_t>(hdr1.sign ^ hdr2.sign), 
		.size = hdr1.size + hdr2.size - digcnt{1} };

	/* If either operand is zero, the result is zero */
	if (hdr1.size == digcnt{0} || hdr2.size == digcnt{0}) {
		set_i64(dst, 0);
		return;
	} 

	if (digprec < hdr.size) {
		set_overflow_flag(dst);
		return;
	}

	// Writable pointers to temporary copies of operands, allowing us to 
	// free them later

	dstptr op1_tmp = 0;
	dstptr op2_tmp = 0;

	if (dst == op1) {
		op1_tmp = new digit[hdr1.size.nd + 1];
		set(op1_tmp, op1, hdr1.size);
		op1 = op1_tmp;
	}

	if (dst == op2) {
		op2_tmp = new digit[hdr2.size.nd + 1];
		set(op2_tmp, op2, hdr2.size);
		op2 = op2_tmp;
	}

	u128 r = { 0, 0 };
	digit d1 = digits(op1)[0];
	for (digcnt j{0}; j < hdr2.size; ++j) {
		digit d2 = digits(op2)[j.nd];
		r = addmul_u64_u64(d1, d2, r.high, 0);
		digits(dst)[j.nd] = r.low;
	}

	for (digcnt i{1}; i < hdr1.size; ++i) {
		digits(dst)[(i + hdr2.size).nd - 1] = r.high;
		r.high = 0;
		digit d1 = digits(op1)[i.nd];

		for (digcnt j{0}; j < hdr2.size; ++j) {
			digit d2 = digits(op2)[j.nd];
			r = addmul_u64_u64(d1, d2, r.high, digits(dst)[(i + j).nd]);
			digits(dst)[(i + j).nd] = r.low;
		}
	}

	if (hdr.size + digcnt{1} <= digprec) {
		digits(dst)[hdr.size.nd] = r.high;
		hdr.size += digcnt{(r.high != 0)};
	} else if (r.high != 0) {
		set_overflow_flag(dst);
		return;
	}	

	if (op1_tmp != 0) delete[] op1_tmp;
	if (op2_tmp != 0) delete[] op2_tmp;

	set_header(dst, hdr);
}

constexpr void
mul_i64_i64(dstptr dst, int64_t op1, int64_t op2, digcnt digprec) noexcept
{
	word op1_tmp[2];
	word op2_tmp[2];

	set_i64(op1_tmp, op1);
	set_i64(op2_tmp, op2);

	mul(dst, op1_tmp, op2_tmp, digprec);
}

constexpr void
mul_i64(dstptr dst, srcptr op1, int64_t op2, digcnt digprec) noexcept
{
	word op2_tmp[2];
	set_i64(op2_tmp, op2);

	mul(dst, op1, op2_tmp, digprec);
}

// Shifts

// Double-precision left shift (shift the value of digit bitcnt bits to the left,
// store resulting digit into the low part of the result, and pick up whatever gets 
// shifted out of the digit into the high part of the result). Named after the 
// x86 instruction, even though it behaves somewhat differently. 

constexpr u128
shld(digit digit, bitcnt nbits) noexcept
{
	// Make sure we don't shift right by 64 bits in one go (which would happen if
	// bitcnt.nb == 0)
	bitcnt rshift = BITS_PER_DIGIT - nbits - bitcnt{1};	
	return u128{ .low = digit << nbits.nb, .high = (digit >> rshift.nb) >> 1 };
}

constexpr void
shl(dstptr dst, srcptr op, bitcnt nbits, bitcnt prec) noexcept
{
	header hdr = get_header(op);
	check_overflow_flag(hdr);

	if (hdr.size == digcnt{0}) {
		set_i64(dst, 0);
		return;
	}

	bitcnt bitsize = bit_size(op) + nbits;
	
	if (prec < bitsize) {
		set_overflow_flag(dst);
		return;
	}

	header dst_hdr = { .sign = hdr.sign, .size = to_digits(bitsize) };
	set_header(dst, dst_hdr);

	// remaining bits
	nbits = remainder_bits(nbits);

	digcnt i = dst_hdr.size - digcnt{1};
	digcnt j = hdr.size - digcnt{1};

	u128 r = shld(digits(op)[j.nd--], nbits);

	if (r.high != 0) digits(dst)[i.nd--] = r.high;

	digit high = r.low;

	while (j < hdr.size) {
		r = shld(digits(op)[j.nd--], nbits);
		digits(dst)[i.nd--] = r.high | high;
		high = r.low;
	}

	digits(dst)[i.nd--] = high;

	while (i < dst_hdr.size) digits(dst)[i.nd--] = 0;
}

constexpr void
shl_i64(dstptr dst, int64_t op, bitcnt nbits, struct bitcnt prec) noexcept
{
	digit op_tmp[2]{0,0};
	set_i64(op_tmp, op);
	shl(dst, op_tmp, nbits, prec);
}

constexpr u128
shrd(digit d, bitcnt nbits) noexcept
{
	bitcnt lshift = BITS_PER_DIGIT - nbits - bitcnt{1};
	return u128{ .low = (d << lshift.nb) << 1, .high = d >> nbits.nb };
}

constexpr void
shr(dstptr dst, srcptr op, bitcnt nbits, bitcnt prec) noexcept
{
	header hdr = get_header(op);
	check_overflow_flag(hdr);

	bitcnt bitsize = bit_size(op) - nbits;
	
	if (prec < bitsize) {
		set_overflow_flag(dst);
		return;
	}

	header dst_hdr = { .sign = hdr.sign, .size = to_digits(bitsize) };
	set_header(dst, dst_hdr);

	// The number of remaining bits
	bitcnt rem = remainder_bits(nbits);

	digcnt i{0};
	digcnt j = whole_digits(nbits);

	u128 r = shrd(digits(op)[j.nd++], rem);
	digit low = r.high;	// The low bits of r get shifted out of existence

	while (i < dst_hdr.size - digcnt{1}) {
		r = shrd(digits(op)[j.nd++], rem);
		digits(dst)[i.nd++] = r.low | low;
		low = r.high;
	}

	r.low = 0;
	if (j < hdr.size) r = shrd(digits(op)[j.nd], rem);
	digits(dst)[i.nd] = r.low | low;
}

// Division

constexpr FPZ_FORCE_INLINE int64_t
i64_divrem(int64_t *quot, int64_t num, srcptr denom)
{
	header denom_hdr = get_header(denom);
	check_overflow_flag(denom_hdr);

	if (denom_hdr.size > digcnt{1}) {
		*quot = 0;
		return num;
	}

	// s = 0 if signs are equal, -1 otherwise
	uint64_t s = -(denom_hdr.sign != (num < 0));

	uint64_t abs_num = abs_i64(num);
	uint64_t abs_denom = digits(denom)[0];
	uint64_t q = abs_num / abs_denom;

	*quot = __builtin_bit_cast(int64_t, (q ^ s) - s);

	s = -(denom_hdr.sign);
	int64_t d = __builtin_bit_cast(int64_t, (abs_denom ^ s) - s);

	return num - *quot * d;
}

constexpr void
divrem(dstptr quot, dstptr rem, srcptr num, srcptr denom, digcnt quot_prec, digcnt rem_prec) noexcept
{
	header num_hdr = get_header(num);
	header denom_hdr = get_header(denom);
	digcnt denom_size = denom_hdr.size;
	digcnt num_size = num_hdr.size;

	check_overflow_flag_2(num_hdr, denom_hdr);

	if (cmp(num, denom) < 0) {
		if (quot != nullptr) set_i64(quot, 0);
		if (rem != nullptr) set(rem, num, rem_prec);
		return;
	}

	// It is still possible that quot overflows, but that will be detected 
	// later anyway 
	if (quot != nullptr && quot_prec + digcnt{1} < num_size - denom_size) {
		set_overflow_flag(quot);
		return;	
	}

	dstptr rem_tmp = rem;
	digcnt rem_tmp_prec = rem_prec;

	// We need the remainder to have at least one more bit than the denominator,
	// so allocate additional digits if required
	if (rem_tmp == nullptr || rem_prec < denom_size + digcnt{1}) {
		rem_tmp = new digit[denom_size.nd + 2]{};
		rem_tmp_prec = denom_size + digcnt{1};
	}

	srcptr denom_tmp = denom;
	dstptr denom_tmp_0 = nullptr;

	// denom is needed throughout the computation, so if it coincides with
	// any of the operands we want to write into, we have to copy it elsewhere
	if (denom == quot || denom == rem_tmp) {
		denom_tmp_0 = new digit[denom_size.nd + 1];
		set(denom_tmp_0, denom, denom_size);
		denom_tmp = denom_tmp_0;
	}

	if (quot != nullptr) {
		set_i64(quot, 0);
		for (size_type i = 0; digcnt{i} < quot_prec; ++i) digits(quot)[i] = 0;
	}

	set_i64(rem_tmp, 0);
	for (size_type i = 0; digcnt{i} < rem_tmp_prec; ++i) digits(rem_tmp)[i] = 0;

	bitcnt bitsize = bit_size(num);

	for (bitcnt i = bitsize - bitcnt{1}; i < bitsize; --i) {
		shl(rem_tmp, rem_tmp, bitcnt{1}, to_bits(rem_tmp_prec));
		if (get_bit(num, i)) set_bit(rem_tmp, bitcnt{0}, to_bits(rem_tmp_prec));
		
		if (cmp(rem_tmp, denom_tmp) >= 0) {
			sub(rem_tmp, rem_tmp, denom_tmp, rem_tmp_prec);
			if (quot != nullptr) set_bit(quot, i, to_bits(quot_prec));
		}
	}

	if (rem_tmp != rem) {
		if (rem != nullptr) set(rem, rem_tmp, rem_prec);
		delete[] rem_tmp;
	}
	
	if (denom_tmp != denom) delete[] denom_tmp;
} 

constexpr int64_t
divrem_i64(dstptr quot, srcptr num, int64_t denom, digcnt quot_prec) noexcept
{
	digit rem_tmp[2]{0,0};
	digit denom_tmp[2]{0,0};
	set_i64(denom_tmp, denom);

	divrem(quot, rem_tmp, num, denom_tmp, quot_prec, digcnt{1});
	return get_sign(rem_tmp) == POSITIVE ? digits(rem_tmp)[0] : -digits(rem_tmp)[0];
}

// Set the digit values of dst to that of the given digit array, which is assumed to contain
// size digits. These are assumed to be stored from most to least significant, i.e.
// the opposite of the usual order (this is to make ordering of digits in an initialization
// list of digits consistent with specification of strings of base 10 / base 16 digits).
// Return ERROR_OVERFLOW if the array of digits represent an integer greater than 2^prec - 1,
// otherwise set the value and return SUCCESS.

FPZ_FORCE_INLINE constexpr void
set_digits(dstptr dst, const digit *digit_values, digcnt size, bitcnt prec)
{
	digcnt first{0};
	while (first < size && digit_values[first.nd] == 0) ++first;

	size -= first;

	if (size == digcnt{0}) {
		set_i64(dst, 0);
		return;
	}

	digcnt digprec = to_digits(prec);

	bitcnt rem = remainder_bits(prec);
	digit limit = (rem == bitcnt{0}) ? MAX_DIGIT : MAX_DIGIT >> (BITS_PER_DIGIT - rem).nb;

	if (size > digprec || (size == digprec && digit_values[first.nd] > limit)) 
		throw integer_overflow{};

	header hdr = { .sign = POSITIVE, .size = size };
	set_header(dst, hdr);

	// Write digits in reverse order into dst
	for (digcnt i{0}; i < size; ++i) digits(dst)[i.nd] = digit_values[(size - i).nd - 1];
}

FPZ_FORCE_INLINE constexpr int64_t
i64_set_digits(const digit *digit_values, digcnt size, bitcnt prec)
{
	digcnt first{0};
	while (first < size && digit_values[first.nd] == 0) ++first;

	size -= first;

	if (size == digcnt{0}) return 0;

	digit limit = MAX_DIGIT >> (BITS_PER_DIGIT - prec).nb;
	if (size > digcnt{1} || (digit_values[first.nd] > limit)) throw integer_overflow{};

	return digit_values[first.nd];
}

FPZ_FORCE_INLINE constexpr int 
digit_value(int c, unsigned base)
{
	int value = c - '0';
	if (base == 10 || (value >= 0 && value < 10)) return value;

	value = c - 'A';
	if (base == 16 && value >= 0 && value < 6) return 10 + value;

	value = c - 'a';
	if (base == 16 && value >= 0 && value < 6) return 10 + value;

	return -1;
}

FPZ_FORCE_INLINE constexpr bool 
is_digit(int c, unsigned base)
{
	int v = digit_value(c, base);
	return v >= 0 && v < base;
}

FPZ_FORCE_INLINE constexpr int64_t
i64_set_str(const char *str, bitcnt prec)
{
	unsigned char sign = str[0] == '-' ? NEGATIVE : POSITIVE;

	unsigned base = 10;
	if (str[0] == '+' || str[0] == '-') str++;
	if (str[0] == '0' && str[1] == 'x') {
		base = 16;
		str+= 2;
	}

	u128 res{ 0, 0 };

	while (*str != '\0') {
		// Allow blanks anywhere in the string
		if (*str == ' ') {
			str++;
			continue;
		}

		int digit = digit_value(*str, base);
		if (digit == -1) throw invalid_argument{};

		if (std::is_constant_evaluated()) res = addmul_u64_u64(res.low, base, digit, 0);
		else res = fpz_addmul_u64_u64(res.low, base, digit, 0);

		if (res.high != 0) throw integer_overflow{};
		
		str++;
	}

	if (prec < bit_size(res.low)) throw integer_overflow{};
	
	return sign == POSITIVE ? res.low : -res.low;
}

FPZ_FORCE_INLINE constexpr void
set_str(dstptr dst, const char *str, bitcnt prec)
{
	set_i64(dst, 0);
	digcnt digprec = to_digits(prec);

	unsigned char sign = str[0] == '-' ? NEGATIVE : POSITIVE;

	unsigned base = 10;
	if (str[0] == '+' || str[0] == '-') str++;
	if (str[0] == '0' && str[1] == 'x') {
		base = 16;
		str+= 2;
	}

	digit overflow = 0;

	while (*str != '\0') {
		// Allow blanks anywhere in the string
		if (*str == ' ') {
			str++;
			continue;
		}

		int digit = digit_value(*str, base);

		if (FPZ_UNLIKELY(digit == -1)) {
			set_i64(dst, 0);
			throw invalid_argument{};
		}

		if (std::is_constant_evaluated()) {
			// If any of these overflow, the computation simply stops
			mul_i64(dst, dst, base, digprec);
			add_i64(dst, dst, digit, digprec);
		} else {
			overflow = fpz_addmul_u64_self(dst, base, digit, digprec);
		}

		if (FPZ_UNLIKELY(overflow != 0)) {
			set_i64(dst, 0);
			throw integer_overflow{};
		}
		
		str++;
	}

	if (std::is_constant_evaluated()) overflow = prec < bit_size(dst);
	else overflow = prec < fpz_bit_size(dst);

	if (FPZ_UNLIKELY(overflow != 0)) {
		set_i64(dst, 0);
		throw integer_overflow{};
	}
	
	set_sign(dst, sign);
}


} // namespace fpz::core

#endif
