#include "test.h"

namespace fpz::test {

constexpr size_t num_tests = 4096;

using core::divrem, core::divrem_i64;
using core::fpz_divrem, core::fpz_divrem_i64;

void
test_divrem()
{
	BEGIN_TEST(divrem);

	for (size_t i = 0; i < num_tests; ++i) {
		auto n = make_random_sparse<512>();
		auto d = make_random_sparse<256>();

		if (d == 0) {
			--i;
			continue;
		}

		ztype<512> q;
		ztype<256> r;

		divrem(q, r, n, d);

		auto v = ztype_from_expr(q*d + r);

		EXPECT(n == v);
	}

	END_TEST;
}

} // namespace fpz::test

int
main(void)
{
	fpz::test::test_divrem();
}
