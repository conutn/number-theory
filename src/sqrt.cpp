// final step: pick square-root method + turn dependency into factor of n

#include "sqrt.hpp"
#include "sqrt_newton.hpp"

// tries one dependency: return true and set factor if it gives a nontrivial factor of n
bool try_dependency(const vector<nfs_relation>& rels, const vector<size_t>& dep, const poly& f, const factor_bases& fb, const bigint& m, const bigint& n, bigint& factor)
{
	return try_dependency_newton(rels, dep, f, fb, m, n, factor);
}
