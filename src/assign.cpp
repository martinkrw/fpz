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

#include "fpz_core.h"

namespace fpz::core {

extern "C" void
fpz_set_max(dstptr dst, bitcnt bitprec) noexcept
{
	set_max(dst, bitprec);
}

extern "C" void
fpz_set(dstptr dst, srcptr src, digcnt digprec) noexcept
{
	set(dst, src, digprec);
}

// Set the value of dst to the value of src, truncated to the given
// target precision if it is too large.

extern "C" void
fpz_truncate(dstptr dst, srcptr src, bitcnt bitprec) noexcept
{
	truncate(dst, src, bitprec);
}

extern "C" void
fpz_set_digits(dstptr dst, const digit *digits, digcnt size, bitcnt bitprec)
{
	set_digits(dst, digits, size, bitprec);
}

extern "C" int64_t
fpz_i64_set_digits(const digit *digits, digcnt size, bitcnt bitprec)
{
	return i64_set_digits(digits, size, bitprec);
}

extern "C" void
fpz_set_str(dstptr dst, const char *str, bitcnt bitprec)
{
	return set_str(dst, str, bitprec);
}

extern "C" int64_t
fpz_i64_set_str(const char *str, bitcnt bitprec) 
{
	return i64_set_str(str, bitprec);
}

// Write a uniformly distributed random integer of the 
// given bit precision into dst. 

extern "C" void
fpz_random(dstptr dst, bitcnt bitprec) noexcept
{
	digcnt size = to_digits(bitprec);
	bitcnt rem = remainder_bits(bitprec);

	arc4random_buf(digits(dst), size.nd * sizeof *dst);

	if (rem != bitcnt{0}) digits(dst)[size.nd - 1] &= (digit{1} << rem.nb) - 1;

	while (size != digcnt{0} && digits(dst)[size.nd - 1] == 0) --size;

	header hdr = { .sign = static_cast<uint8_t>(arc4random() & 1), .size = size };
	set_header(dst, hdr);
}


// Return a uniformly distributed random 
// integer of precision bitprec < 64.

extern "C" int64_t
fpz_random_i64(bitcnt bitprec) noexcept
{
	uint64_t v;
	arc4random_buf(&v, sizeof v);

	uint64_t s = -(v & 1);	// all 1's (-1) if first bit is set, all 0's otherwise
	v >>= (BITS_PER_DIGIT - bitprec).nb;

	return __builtin_bit_cast(int64_t, (v ^ s) - s);
}

// Similar to fpz_random, but introduces more contiguous blocks of zero bits than would otherwise 
// typically occur. Probably only useful for testing.

extern "C" void
fpz_random_sparse(dstptr dst, bitcnt bitprec)  noexcept
{
	header hdr = { .sign = static_cast<uint8_t>(arc4random() & 1), .size = digcnt{0} };

	digcnt max_size = to_digits(bitprec);

	for (digcnt i{0}; i < max_size; ++i) {
		digit d, mask;
		arc4random_buf(&d, sizeof d);

		uint32_t r = arc4random();
	
		mask = (r % 4) ? MAX_DIGIT : 0;	// zero out every 4th digit
		mask >>= ((r >> 2) % 4) * 16;
		mask <<= ((r >> 4) % 4) * 16;

		d &= mask;
		digits(dst)[i.nd] = d;

		if (d != 0) hdr.size = i + digcnt{1};
	}

	set_header(dst, hdr);
}

} // namespace fpz::core
