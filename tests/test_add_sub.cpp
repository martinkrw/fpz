#include "fpz.h"
#include "test.h"

namespace fpz::test {

constexpr size_t num_tests = 8192;

using core::add_i64_i64, core::sub_i64_i64, core::add_i64, core::sub_i64, core::i64_sub;
using core::add, core::sub, core::fpz_add, core::fpz_sub;
using core::fpz_add_i64_i64, core::fpz_sub_i64_i64, core::fpz_add_i64, core::fpz_sub_i64, core::fpz_i64_sub;

constexpr int64_t int64_min = std::numeric_limits<int64_t>::min();

void
test_add_sub_i64_i64()
{
	BEGIN_TEST(add_i64_i64/sub_i64_i64);

	// Make sure int64_min is dealt with correctly
	constexpr auto sum1 = ztype_from_expr(const_i64<int64_min> + const_i64<int64_min>);
	auto sum2 = ztype_from_expr(const_i64<int64_min> + const_i64<int64_min>);

	EXPECT(sum1 == sum2);
	EXPECT(sum1 == ztype_from_expr(-const_digits<1, 0>));

	// Make sure overflow flag is correctly set
	ztype<64> z1, z2;
	add_i64_i64(get_mutable_value(z1), int64_min, int64_min, to_digits(z1.prec));
	fpz_add_i64_i64(get_mutable_value(z2), int64_min, int64_min, to_digits(z2.prec));

	EXPECT_OVERFLOW(get_value(z1));
	EXPECT_OVERFLOW(get_value(z2));

	for (size_t i = 0; i < 64; ++i) {
		int64_t x = core::fpz_random_i64(bitcnt{63});
		int64_t y = core::fpz_random_i64(bitcnt{63});

		add_i64_i64(get_mutable_value(z1), x, y, to_digits(z1.prec));
		fpz_add_i64_i64(get_mutable_value(z2), x, y, to_digits(z2.prec));
		EXPECT_NO_OVERFLOW(get_value(z1));
		EXPECT_NO_OVERFLOW(get_value(z2));

		EXPECT(z1 == z2);

		sub_i64_i64(get_mutable_value(z1), x, y);
		fpz_sub_i64_i64(get_mutable_value(z2), x, y);

		EXPECT(z1 == z2);

		// Make sure x +/- y doesn't overflow (with built-in arithmetic)
		x /= 2;
		y /= 2;
	
		fpz_add_i64_i64(get_mutable_value(z1), x, y, to_digits(z2.prec));
		fpz_sub_i64_i64(get_mutable_value(z2), x, y);

		EXPECT_NO_OVERFLOW(get_value(z1));
		EXPECT_NO_OVERFLOW(get_value(z2));

		EXPECT(z1 == (x + y));
		EXPECT(z2 == (x - y));
	}

	END_TEST;
}

void
test_add_sub_i64()
{
	BEGIN_TEST(add_i64/sub_i64/i64_sub);
	
	auto m = ztype<256>::make_max();

	ztype<256> z1, z2;

	add_i64(get_mutable_value(z1), get_value(m), 1, to_digits(z1.prec));
	fpz_add_i64(get_mutable_value(z2), get_value(m), 1, to_digits(z2.prec));

	EXPECT_OVERFLOW(get_value(z1));
	EXPECT_OVERFLOW(get_value(z2));

	sub_i64(get_mutable_value(z1), get_value(m), -1, to_digits(z1.prec));
	fpz_sub_i64(get_mutable_value(z2), get_value(m), -1, to_digits(z2.prec));
	
	EXPECT_OVERFLOW(get_value(z1));
	EXPECT_OVERFLOW(get_value(z2));
	
	m.negate();

	i64_sub(get_mutable_value(z1), 1, get_value(m), to_digits(z1.prec));
	fpz_i64_sub(get_mutable_value(z1), 1, get_value(m), to_digits(z2.prec));

	EXPECT_OVERFLOW(get_value(z1));
	EXPECT_OVERFLOW(get_value(z2));

	for (size_t i = 0; i < num_tests; ++i) {
		auto x = make_random_sparse<256>();
		int64_t y = fpz_random_i64(bitcnt{63});

		ztype<257> z1, z2;

		add_i64(get_mutable_value(z1), get_value(x), y, to_digits(z1.prec));
		fpz_add_i64(get_mutable_value(z2), get_value(x), y, to_digits(z2.prec));

		EXPECT(z1 == z2);

		sub_i64(get_mutable_value(z1), get_value(x), y, to_digits(z1.prec));
		fpz_sub_i64(get_mutable_value(z2), get_value(x), y, to_digits(z1.prec));

		EXPECT(z1 == z2);

		i64_sub(get_mutable_value(z1), y, get_value(x), to_digits(z1.prec));
		fpz_i64_sub(get_mutable_value(z2), y, get_value(x), to_digits(z1.prec));
	
		EXPECT(z1 == z2);
	}

	END_TEST;
}

void
test_add_sub()
{
	BEGIN_TEST(add/sub);

	auto m1 = ztype<256>::make_max();
	auto m2 = ztype<256>::make_max();

	add(get_mutable_value(m1), get_value(m1), get_value(m1), to_digits(m1.prec));
	fpz_add(get_mutable_value(m2), get_value(m2), get_value(m2), to_digits(m2.prec));

	EXPECT_OVERFLOW(get_value(m1));
	EXPECT_OVERFLOW(get_value(m2));

	m1.set_max();
	m2.set_max();
	m2.negate();

	sub(get_mutable_value(m1), get_value(m1), get_value(m2), to_digits(m1.prec));
	EXPECT_OVERFLOW(get_value(m1));

	m1.set_max();
	fpz_sub(get_mutable_value(m1), get_value(m2), get_value(m1), to_digits(m1.prec));
	EXPECT_OVERFLOW(get_value(m1));

	for (size_t i = 0; i < num_tests; ++i) {
		auto x = make_random_sparse<512>();
		auto y = make_random_sparse<384>();

		ztype<513> z1, z2;

		add(get_mutable_value(z1), get_value(x), get_value(y), to_digits(z1.prec));
		fpz_add(get_mutable_value(z2), get_value(x), get_value(y), to_digits(z2.prec));

		EXPECT_NO_OVERFLOW(get_value(z1));
		EXPECT_NO_OVERFLOW(get_value(z2));
		EXPECT(z1 == z2);

		sub(get_mutable_value(z1), get_value(y), get_value(x), to_digits(z1.prec));
		fpz_sub(get_mutable_value(z2), get_value(y), get_value(x), to_digits(z2.prec));

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
	fpz::test::test_add_sub_i64_i64();
	fpz::test::test_add_sub_i64();
	fpz::test::test_add_sub();
}


