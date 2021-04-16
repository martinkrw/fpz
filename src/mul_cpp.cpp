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

#include <algorithm>	// for std::swap
#include <cstdlib>

#include "fpz_core.h"

namespace fpz::core {

extern "C" {

void fpz_umul_sp(dstptr, srcptr, srcptr, digcnt, digcnt, digcnt) noexcept;
void fpz_umul_mp(dstptr, srcptr, srcptr, digcnt, digcnt, digcnt) noexcept;

}

extern "C" void
fpz_mul(dstptr dst, srcptr op1, srcptr op2, digcnt digprec) noexcept
{
	header hdr1 = get_header(op1);
	header hdr2 = get_header(op2);

	check_overflow_flag_2(hdr1, hdr2);

	set_sign(dst, static_cast<uint8_t>(hdr1.sign ^ hdr2.sign));

	digcnt n1 = hdr1.size;
	digcnt n2 = hdr2.size;

	if (n1 < n2) {
		std::swap(n1, n2);
		std::swap(op1, op2);
	}

	if (n2 == digcnt{0}) {
		set_i64(dst, 0);
		return;
	} else if (digprec < n1 + n2 - digcnt{1}) {
		set_overflow_flag(dst);
		return;
	}

	if (n2 <= digcnt{4}) {
		fpz_umul_sp(dst, op1, op2, digprec, n1, n2);
		return;
	}

	// If the computation is going to take multiple passes, the destination
	// operand can't overlap with any of the source operands, so copy 
	// to stack if this is the case.

	if (dst == op1) {
		dstptr op1_tmp;
		size_t sz = int_alloc_size(n1);
		op1_tmp = static_cast<dstptr>(alloca(sz));
		op1 = static_cast<srcptr>(std::memcpy(op1_tmp, op1, sz));
	}

	if (dst == op2) {
		dstptr op2_tmp;
		size_t sz = int_alloc_size(n2);
		op2_tmp = static_cast<dstptr>(alloca(sz));
		op2 = static_cast<srcptr>(std::memcpy(op2_tmp, op2, sz));
	}

	fpz_umul_mp(dst, op1, op2, digprec, n1, n2);
	return;
}

} // namespace fpz::core
