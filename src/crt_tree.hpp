#pragma once
// product and remainder trees + CRT reconstruction

#include "common.hpp"

bigint balanced_product(vector<bigint> v);

// value mod p for each prime
vector<u64> reduce_all(const bigint& value, const vector<u64>& primes);

// returns vector X, |X[c]|<P/2, s.t. X[c] = res[i][c] (mod primes[i])
// for all i and c < d
vector<bigint> crt_tree(const vector<u64>& primes, const vector<vector<u64>>& res, size_t d);
