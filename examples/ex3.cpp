#include "fpz.h"

// Compute a 4x4 determinant using a range-based for-loop.

using fpz::precision;
using fpz::const_i64;
using fpz::ztype;

// Separate namespace for precision calculations

namespace precisions {
	constexpr precision entry = 128_bits;	// row entries are 128 bits

	// The total determinant is a sum of terms of the form below 

	constexpr precision det4x4_term = (entry*entry - entry*entry)*(entry*entry - entry*entry);
	constexpr precision det4x4 = sum(det4x4_term, 6);
}

// precision implicitly converts to uint32_t

using entry_type = ztype<precisions::entry>;
using det4x4_term = ztype<precisions::det4x4_term>;
using det4x4_type = ztype<precisions::det4x4>;

struct row4 {
	entry_type c[4];
	row4() : c { entry_type::make_random(), entry_type::make_random(), 
		entry_type::make_random(), entry_type::make_random() } {}
};

static void
compute_det4x4(det4x4_type& out, const row4& row0, const row4& row1, const row4& row2, const row4& row3)
{
	int perm[6][4] = { {0,1,2,3}, {0,2,1,3},
			  {0,3,1,2}, {1,2,0,3},
			  {1,3,0,2}, {2,3,0,1} };

	int signatures[6] = { 1, -1, 1, 1, -1, 1 };

	int i = 0;

	// Add together the values assigned to the term variable for each iteration
	// and store the result into out. 

	// Pass our own term variable to the loop
	det4x4_term t;
	for (auto [term, value] : out.sum<6>(const_i64<0>, t)) {
		int s0 = perm[i][0];
		int s1 = perm[i][1];
		int s2 = perm[i][2];
		int s3 = perm[i][3];

		assert(&term == &t);

		term = (row0.c[s0]*row1.c[s1] - row0.c[s1]*row1.c[s0]) *
		    (row2.c[s2]*row3.c[s3] - row2.c[s3]*row3.c[s2]);
		if (signatures[i] == -1) term.negate();

		printf("iteration %d, current value = ", i);
		value.println(); // Note: value is read-only

		++i;
	}
}

int
main(void)
{
	row4 *rows = new row4[3];
	
	det4x4_type out;
	// repeat the second row, to make the determinant come out to zero
	compute_det4x4(out, rows[0], rows[1], rows[1], rows[2]);

	printf("value of determinant: ");
	out.println();	// should be zero 

	delete[] rows;
}

