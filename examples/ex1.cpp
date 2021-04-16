#include "fpz.h"

using fpz::ztype;
using fpz::const_i64;
using fpz::const_digits;

// Runtime arithmetic

int
main(void)
{
	// Get two random, 128 bit integers
	const auto z1 = ztype<128>::make_random();
	const auto z2 = ztype<128>::make_random();

	// Multiply them (the precision of product is automatically deduced)
	auto product = ztype_from_expr(z1*z2);

	std::printf("z1*z2 = ");
	product.println();

	// Conversion to floating point
	std::printf("(long double)(z1*z2) = %Lf\n", static_cast<long double>(product));
	
	// Some simple checks to see if arithmetic works as expected
	auto z3 = ztype_from_expr(z1*z2 + z1);
	auto z4 = ztype_from_expr((z3 - z1) / z1);

	auto z5 = ztype_from_expr(z1);
	z5.println();

	if (z4 == z2) std::printf("correct!\n");
	else std::printf("incorrect\n");

	auto z6 = ztype_from_expr(const_digits<1,0> + const_digits<5>);
	z6.println();

	const auto x = ztype<512>::make_random();
	const auto y = ztype<256>::make_random();

	auto q = ztype_from_expr(x / y);
	auto r = ztype_from_expr(x % y);

	if (ztype_from_expr(q*y + r) == x) std::printf("correct!\n");
	else std::printf("incorrect!\n");

	if (ztype_from_expr(z1 << const_i64<65>) == ztype_from_expr(z1*const_i64<0x4000000000000000>*const_i64<8>)) 
		std::printf("correct!\n");
	else std::printf("incorrect!\n");

}
