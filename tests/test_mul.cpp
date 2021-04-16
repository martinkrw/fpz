#include "test.h"

namespace fpz::test {

constexpr size_t num_tests = 8192;

using core::mul_i64_i64, core::mul_i64, core::mul;
using core::fpz_mul_i64_i64, core::fpz_mul_i64, core::fpz_mul;

constexpr int64_t int64_min = std::numeric_limits<int64_t>::min();

void
test_mul_i64_i64()
{
	BEGIN_TEST(mul_i64_i64);

	// Make sure overflow flag is correctly set

	ztype<64> z1, z2;
	mul_i64_i64(get_mutable_value(z1), int64_min, int64_min, to_digits(z1.prec));
	fpz_mul_i64_i64(get_mutable_value(z2), int64_min, int64_min, to_digits(z2.prec));

	EXPECT_OVERFLOW(get_value(z1));
	EXPECT_OVERFLOW(get_value(z2));

	for (size_t i = 0; i < 64; ++i) {
		int64_t x = core::fpz_random_i64(bitcnt{63});
		int64_t y = core::fpz_random_i64(bitcnt{63});

		ztype<128> z1, z2;

		mul_i64_i64(get_mutable_value(z1), x, y, to_digits(z1.prec));
		fpz_mul_i64_i64(get_mutable_value(z2), x, y, to_digits(z2.prec));

		EXPECT_NO_OVERFLOW(get_value(z1));
		EXPECT_NO_OVERFLOW(get_value(z2));

		EXPECT(z1 == z2);

		// Make sure x * y doesn't overflow (with built-in arithmetic)
		x /= (uint64_t{1} << 33);
		y /= (uint64_t{1} << 33);
	
		fpz_mul_i64_i64(get_mutable_value(z2), x, y, to_digits(z2.prec));
		EXPECT_NO_OVERFLOW(get_value(z2));

		EXPECT(z2 == x*y);
	}

	END_TEST;
}

void
test_mul_i64()
{
	BEGIN_TEST(mul_i64);
	
	auto m = ztype<256>::make_max();

	ztype<256> z1, z2;

	mul_i64(get_mutable_value(z1), get_value(m), 999999, to_digits(z1.prec));
	fpz_mul_i64(get_mutable_value(z2), get_value(m), 999999, to_digits(z2.prec));

	EXPECT_OVERFLOW(get_value(z1));
	EXPECT_OVERFLOW(get_value(z2));

	for (size_t i = 0; i < num_tests; ++i) {
		auto x = make_random_sparse<256>();
		int64_t y = fpz_random_i64(bitcnt{63});

		ztype<320> z1, z2;

		mul_i64(get_mutable_value(z1), get_value(x), y, to_digits(z1.prec));
		fpz_mul_i64(get_mutable_value(z2), get_value(x), y, to_digits(z2.prec));

		EXPECT(z1 == z2);
	}

	END_TEST;
}

void
test_mul()
{
	BEGIN_TEST(mul);

	auto m1 = ztype<256>::make_max();
	auto m2 = ztype<256>::make_max();

	mul(get_mutable_value(m1), get_value(m1), get_value(m1), to_digits(m1.prec));
	fpz_mul(get_mutable_value(m2), get_value(m2), get_value(m2), to_digits(m2.prec));

	EXPECT_OVERFLOW(get_value(m1));
	EXPECT_OVERFLOW(get_value(m2));

	for (size_t i = 0; i < num_tests; ++i) {
		auto x = make_random_sparse<512>();
		auto y = make_random_sparse<192>();
		auto z = make_random_sparse<512>();

		ztype<1024> z1, z2;

		mul(get_mutable_value(z1), get_value(x), get_value(y), to_digits(z1.prec));
		fpz_mul(get_mutable_value(z2), get_value(x), get_value(y), to_digits(z2.prec));

		EXPECT_NO_OVERFLOW(get_value(z1));
		EXPECT_NO_OVERFLOW(get_value(z2));
		EXPECT(z1 == z2);

		mul(get_mutable_value(z1), get_value(x), get_value(z), to_digits(z1.prec));
		fpz_mul(get_mutable_value(z2), get_value(x), get_value(z), to_digits(z2.prec));

		EXPECT_NO_OVERFLOW(get_value(z1));
		EXPECT_NO_OVERFLOW(get_value(z2));
		EXPECT(z1 == z2);
	}

	END_TEST;
}

} // namespace fpz::test

int
main(void)
{
	fpz::test::test_mul_i64_i64();
	fpz::test::test_mul_i64();
	fpz::test::test_mul();
}

