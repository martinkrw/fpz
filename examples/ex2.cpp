#include "fpz.h"

using fpz::ztype_from_expr;
using fpz::ztype;
using fpz::const_i64;
using fpz::const_digits;

// Compile-time arithmetic

constexpr auto z1 = ztype_from_expr(FPZ_CONST_STR("1 000 000 000 000 000 000 001 000"));

constexpr auto z2 = ztype_from_expr(z1 + const_i64<4>*z1 + const_i64<2>*z1);
constexpr auto z3 = ztype_from_expr(z2*z2 + const_i64<42>);
constexpr auto z4 = ztype_from_expr(const_digits<1, 0>); // set the value to 2^64

constexpr auto z5 = ztype_from_expr(FPZ_CONST_STR("1 000 000 000 000 000 000 000 000 000 000") / 
	FPZ_CONST_STR("1 000 000 000 000 000 000 000 000"));

constexpr auto z6 = ztype_from_expr(const_i64<1> << const_i64<64>);

static_assert(z4 == ztype_from_expr(const_i64<1> << const_i64<64>), "incorrect!");

static_assert(z1 > const_i64<-5>, "incorrect!");

int
main(void)
{
	z1.println();
	z2.println();
	z3.println();
	z4.println();
	z6.println();
	z5.println();
}
