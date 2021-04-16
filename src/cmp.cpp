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

extern "C" bitcnt
fpz_bit_size(srcptr op) noexcept
{
	return bit_size(op);
}

extern "C" int
fpz_signum(srcptr op) noexcept 
{
	return signum(op);
}

int
fpz_cmp(srcptr op1, srcptr op2) noexcept
{
	return cmp(op1, op2);
}


int
fpz_cmp(srcptr op1, int64_t op2) noexcept
{
	return cmp(op1, op2);
}

int
fpz_cmpabs(srcptr op1, srcptr op2) noexcept
{
	return cmpabs(op1, op2);
}

int
fpz_cmpabs(srcptr op1, int64_t op2) noexcept
{
	return cmpabs(op1, op2);
}

int
fpz_cmpabs(int64_t op1, int64_t op2) noexcept
{
	return cmpabs(op1, op2);
}

} // namespace fpz::core
