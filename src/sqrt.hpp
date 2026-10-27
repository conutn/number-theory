#pragma once
// the final step: pick the square-root method and turn a dependency into a factor of n

#include "common.hpp"
#include "factor_base.hpp"
#include "sieve.hpp"

// try one dependency
// odd degree: Couveignes, even degree: p-adic Newton fallback
// returns true if nontrivial factor found, sets factor
bool try_dependency(const vector<nfs_relation>& rels, const vector<size_t>& dep, const poly& f, const factor_bases& fb, const bigint& m, const bigint& n, bigint& factor);
