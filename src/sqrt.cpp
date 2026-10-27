// the final step: pick the square-root method and turn a dependency into a factor of n

#include "sqrt.hpp"
#include "sqrt_newton.hpp"
#include "sqrt_couveignes.hpp"

bool try_dependency(const vector<nfs_relation>& rels, const vector<size_t>& dep, const poly& f, const factor_bases& fb, const bigint& m, const bigint& n, bigint& factor)
{
	// the sign trick of Couveignes needs odd degree
	if ((f.size() - 1) % 2 == 1) return try_dependency_couveignes(rels, dep, f, fb, m, n, factor);
	return try_dependency_newton(rels, dep, f, m, n, factor);
}
