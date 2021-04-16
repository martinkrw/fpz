#include <cinttypes>
#include "fpz.h"

// Example showing the possibility of defining custom expression types
// (this is really no different from the built-in ones)

using fpz::ztype;
using fpz::ztype_from_expr;
using fpz::detail::Expr;
using fpz::detail::expr_traits_nocvref;
using fpz::detail::expr_traits;
using fpz::detail::eval_any;
using fpz::detail::is_large;

template <Expr E>
struct example_expr {
	const E& e;
	
	template <typename D> FPZ_FORCE_INLINE constexpr auto
	eval(D dst) const noexcept 
	{
		using Tr = expr_traits<example_expr>;
		puts("Evaluating example_expr ... ");
		if constexpr (is_large(Tr::bit_bound)) {
			eval_any(e, dst);
			printf("If this was not a toy example, we would\n");
			printf("\t- evaluate e,\n");
			printf("\t- transform the result into another integer (requiring at most %" PRIu32 " bits), and\n", D::prec.nb);
			printf("\t- write it into %p.\n", static_cast<void*>(dst.ptr));
			return dst.ptr;
		} else {
			printf("If this was not a toy example, we would\n");
			printf("\t- evaluate e,\n");
			printf("\t- transform the result into an int64_t (requiring at most %" PRIu32 " bits), and\n", D::prec.nb);
			printf("\t- return it.\n");
			return eval_any(e, nullptr);
		}

		eval_any(e, dst);

	}
};

// Specialize the expression traits template for our new expression type

template <typename E>
struct expr_traits_nocvref<example_expr<E>> : public fpz::detail::expr_traits_compound_base {

	// Precision bounds are the same as for E, since we're not doing anything interesting
	static constexpr bitcnt bit_bound = expr_traits<E>::bit_bound;
	static constexpr expr_bound bound = expr_traits<E>::bound;

	template <typename D> static constexpr bool subexpr_may_alias_with = 
		expr_traits<E>::template subexpr_may_alias_with<D>;
};

template <Expr E> constexpr auto
example_operator(const E& e)
{
	return example_expr<E> { e };
}

int
main()
{
	auto x = ztype<128>::from_i64(42);
	auto y = ztype_from_expr(abs(x - example_operator(x + x)));
	y.println();
}
