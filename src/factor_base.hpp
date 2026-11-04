#pragma once
// factor bases (rational, algebraic w/ projective roots, quadratic chars), bound heuristics, norms

#include "common.hpp"

inline constexpr uint64_t FB_PROJECTIVE = UINT64_MAX; // for p | lc(f)

struct fb_pair
{
	uint64_t p;
	uint64_t r; // rational: m mod p, algebraic/quadratic: f(r)=0 mod p (or FB_PROJECTIVE)
};

struct factor_bases
{
	vector<fb_pair> rational;
	vector<fb_pair> algebraic; // w/ projective entries
	vector<fb_pair> quadratic;
	// sign is a separate column not stored here
	bool parity = false; // extra column w/ all ones, all dependencies then have even number of relations
	                     // necessary when lc(f)!=1: rational side of square root gets a factor lc^k, square only for even k (every element of algebraic relation contains lc)
};

// smoothness bound for rational + algebraic
uint64_t choose_bound(const bigint& n);

// # of quadratic char pairs
size_t choose_num_qc(const bigint& n);

// build all three bases
// rat_bound, alg_bound: primes up to these
// num_qc: # of (q, s) pairs, w/ q > max(rat_bound, alg_bound)
factor_bases make_factor_bases(const bigint& n, const poly& f, const bigint& m, uint64_t rat_bound, uint64_t alg_bound, size_t num_qc);

// N(a+b*theta)
bigint norm(poly f, bigint a, bigint b);

// # of matrix columns: sign + all bases (+ parity)
size_t num_columns(const factor_bases& fb);
