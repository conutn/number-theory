#pragma once
// final step: pick square-root method + turn dependency into factor of n

#include "common.hpp"
#include "factor_base.hpp"
#include "sieve.hpp"

// try one dependency
// uses Newton square root for large factor bases, Couveigne's is O(n^3/2)
// return true if nontrivial factor found & set factor
bool try_dependency(const vector<nfs_relation>& rels, const vector<size_t>& dep, const poly& f, const factor_bases& fb, const bigint& m, const bigint& n, bigint& factor);
