#pragma once
// product and remainder trees, and the Chinese remainder reconstruction built on them
// product / remainder trees.  With k primes the CRT used to cost O(k^2) word operations (Garner);
// with balanced trees it is dominated by a few big multiplications and divisions

#include "common.hpp"

bigint balanced_product(vector<bigint> v);

// value mod p for each prime
vector<u64> reduce_all(const bigint& value, const vector<u64>& primes);

// the integer x, |x| < P/2, with x = res[i][c] (mod primes[i]) for each coefficient c
//   x = sum_i c_i * (P/p_i)  (mod P),  c_i = res_i * ((P/p_i)^-1 mod p_i)
// (P/p_i) mod p_i comes from  P mod p_i^2 = p_i * ((P/p_i) mod p_i)  (remainder tree),
// the sum from the tree recurrence S = S_L * P_R + S_R * P_L
vector<bigint> crt_tree(const vector<u64>& primes, const vector<vector<u64>>& res, size_t d);
