#pragma once
// relation finding: line sieve on both sides, candidate factoring (trial or re-sieve), quadratic characters
//
// column format of the mod 2 matrix (one row per relation (a,b)):
//  0               sign of a+b*m
//  1 ... R         rational base
//  R+1 ... R+A     algebraic base (projective primes p | lc(f) included)
//  R+A+1 ... R+A+Q quadratic characters
//  R+A+Q+1         parity: 1 in every row (only if lc(f) != 1, makes every dependency even)

#include "common.hpp"
#include "factor_base.hpp"

// a prime of a factorization: index into fb.rational (rfac) or fb.algebraic (afac) and its FULL exponent
// (cols only keeps the parities). the square root needs the real exponents, prod q^(e/2),
// so keeping them here means it never has to factor anything again
struct fac_entry
{
	uint32_t idx;
	uint32_t e;
};

struct nfs_relation
{
	int64_t a;
	int64_t b; // b>0, gcd(a,b)=1
	vector<uint32_t> cols; // columns where entry is 1 (column 0 first if a+b*m < 0)
	vector<fac_entry> rfac; // rational factorization of |a+b*m|
	vector<fac_entry> afac; // algebraic: v_p(N(a+b*theta)) for each (p,r) base entry dividing it
	uint64_t lp_r = 0; // large prime on the rational side, 0 = none
	uint64_t lp_a = 0; // large prime on the algebraic side, 0 = none
	uint64_t lp_a_root = 0; // root r of the algebraic ideal (lp_a, r)
	vector<uint32_t> members; // combined relations: indices of the partial relations it is made of
};

struct sieve_params
{
	int64_t A = 0; // try a in [-A, A], 0 -> auto
	uint64_t b_max = 0; // 0 -> auto (the caller retries with bigger parameters if this is not enough)
	size_t target = 0; // # of relations desired, 0 -> #columns+10
	float slack = 1.0f; // fudge factor (natural log) for the rational side
	float aslack = 3.0f; // same for the algebraic side (prime powers are not sieved there)
	int cand_mode = 0; // how candidates are factored: 0 = choose per line, 1 = always test every base prime, 2 = always re-sieve
	uint64_t lp_bound = 0; // large prime bound on each side, 0 -> no large primes
};

inline constexpr uint64_t DEFAULT_B_MAX = 5000;

// appends relations to out until out.size() >= sp.target or b exceeds b_max
// relations with a large prime go to partials
// next_b is the first b still to be processed (in/out), so the sieve can be resumed with a larger target
void find_relations(const poly& f, const bigint& m, const factor_bases& fb, sieve_params sp, vector<nfs_relation>& out, vector<nfs_relation>& partials, uint64_t& next_b);

// number of independent cycles = number of full relations obtainable by combining the partials
size_t count_lp_cycles(const vector<nfs_relation>& partials);

// append one combined relation per fundamental cycle of the large-prime graph to out
void combine_partials(const vector<nfs_relation>& partials, vector<nfs_relation>& out);
