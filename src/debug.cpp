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

#include <cinttypes>

#include "fpz_core.h"

namespace fpz::core {

extern "C" bool
fpz_is_valid(srcptr op, bitcnt bitprec) noexcept
{
	return is_valid(op, bitprec);
}

extern "C" void
fpz_fprint_debug(FILE* stream, srcptr src, bitcnt bitprec) noexcept
{
	header hdr = get_header(src);
	std::fprintf(stream, "{ .sign = %" PRIu8 ", .of = %" PRIu8 ", .size = %" PRIu32 "}, capacity = %" PRIu32 " bits", 
	    hdr.sign, hdr.overflow_flag, hdr.size.nd, bitprec.nb);

	digcnt digprec = to_digits(bitprec);
	if (hdr.size > digprec || bitprec < bit_size(src)) std::fputs(" (overflowed)", stream);
	std::fputc('\n', stream);

	if (hdr.size == digcnt{0}) std::fprintf(stream, "%.16" PRIX64, digit{0});

	digcnt sz = min(hdr.size, digprec);
	for (digcnt i{0}; i < sz; ++i) std::fprintf(stream, "%.16" PRIX64 " ", digits(src)[i.nd]);
	std::fputc('\n', stream);
}

extern "C" void
fpz_assertion_failure(const char *msg, const char *file, int line) noexcept
{
	std::fprintf(stderr, "%s: file %s, line %d\n", msg, file, line);
	std::abort();
}

extern "C" void
fpz_use_after_overflow_handler() noexcept
{
	std::fprintf(stderr, "libfpz: use of overflowed integer variable detected ... aborting\n");
	std::abort();
}

} // namespace fpz::core
