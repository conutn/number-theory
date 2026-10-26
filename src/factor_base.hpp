#pragma once
// factor bases (rational, algebraic incl. projective roots, quadratic characters), bound heuristics, norms
// factor bases
// rational base  r=m mod p
// algebraic base f(r)=0 mod p
// projective     r=FB_PROJECTIVE for p | lc(f)
// quadratic base f(s)=0 mod q
// sign (-1) is a separate column not stored here

#include "common.hpp"

inline constexpr uint64_t FB_PROJECTIVE = UINT64_MAX;

struct fb_pair
{
	uint64_t p;
	uint64_t r;
};

struct factor_bases
{
	vector<fb_pair> rational;
	vector<fb_pair> algebraic; // includes projective entries
	vector<fb_pair> quadratic;
	bool parity = false; // extra last column of all ones: every dependency then has an even number of relations.
	                     // needed when lc(f) != 1: the rational side of the square root gets a factor lc^k, a square only for even k
};

uint64_t choose_bound(const bigint& n);

size_t choose_num_qc(const bigint& n);

// build all three bases.
//  rat_bound, alg_bound: primes up to these
//  num_qc: # of (q, s) pairs, with q > max(rat_bound, alg_bound)
factor_bases make_factor_bases(const poly& f, const bigint& m, uint64_t rat_bound, uint64_t alg_bound, size_t num_qc);

// N(a+b*theta)
bigint norm(poly f, bigint a, bigint b);

size_t num_columns(const factor_bases& fb);
