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

#include "fpz.h"

namespace fpz::test {

using core::bitcnt;
using core::to_digits;
using core::size_type;
using fpz::detail::get_value;

// For testing, we want to call low-level functions directly, which 
// requires write access to the actual header/digits. 
 
template <size_type prec> auto
get_mutable_value(ztype<prec>& e) {
	return const_cast<core::dstptr>(get_value(e));
}

template <size_type prec> ztype<prec>
make_random_sparse()
{
	ztype<prec> z;
	fpz_random_sparse(get_mutable_value(z), z.prec);
	return z;
}

} // namespace fpz::test

#define EXPECT(cond)	do {								\
	if (!(cond)) {									\
		printf("FAIL\n\t(Condition %s failed on line %d)\n", #cond, __LINE__);	\
		exit(1);								\
	}										\
} while (false)

#define EXPECT_OVERFLOW(srcptr)		do {					\
	fpz::core::header hdr = fpz::core::get_header(srcptr);			\
	EXPECT(hdr.overflow_flag != 0);						\
} while (false)

#define EXPECT_NO_OVERFLOW(srcptr) 	do {					\
	fpz::core::header hdr = fpz::core::get_header(srcptr);			\
	EXPECT(hdr.overflow_flag == 0);						\
} while (false)

#define BEGIN_TEST(funcname)	printf("Testing %s ... ", #funcname)
#define END_TEST		printf("PASS\n")
