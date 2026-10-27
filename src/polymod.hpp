#pragma once
// polynomials over Z and F_p: gcd, powers mod f, roots and factorization degrees mod p, irreducibility test

#include "common.hpp"

// put coefficients in [0, p), remove leading zeros
poly reduce(poly a, const bigint& p);

// gcd(a,b) in F_p[x], returns monic polynomial
poly poly_gcd(poly a, poly b, const bigint& p);

// a*b mod (f, p)
// requires:
//  a reduced
//  b reduced
//  f reduced & monic & deg>=1
poly mul_mod(const poly& a, const poly& b, const poly& f, const bigint& p);

// distinct solutions of f(x)=0 (mod p)
vector<bigint> roots(const poly& f, const bigint& p);

// b^e mod (f, p)
// requires:
//  b reduced
//  f reduced & monic & deg>=1
poly pow_mod(const poly& b, const bigint& e, const poly& f, const bigint& p);

// derivative over F_p
// requires:
//  f reduced
poly derivative(const poly& f, const bigint& p);

// degrees of irreducible factors of f over F_p
// requires:
//  f reduced & monic & squarefree & deg>=1
vector<size_t> factor_degrees(const poly& f, const bigint& p);

enum class irred
{
	proven, // irreducible over Q
	undetermined // reducible or irreducible (x^4+1)
};

// irreducibility test
// requires:
//  f0 primitive & small degree
// can only prove irreducibility, does not prove reducibility
irred is_irreducible(const poly& f0, size_t num_primes = 30);

// f(x) mod p
bigint eval_mod(const poly& f, const bigint& x, const bigint& p);

// monic version of f(x) = sum c_j x^j, lc = c_d
// alpha = lc*theta is an algebraic integer with the monic minimal polynomial
//   fhat(x) = lc^(d-1) f(x/lc) = x^d + sum_{j<d} c_j lc^(d-1-j) x^j
// (lc = 1: fhat = f)
// relations a+b*theta become lc*a + b*alpha in Z[alpha]; theta -> m corresponds to alpha -> lc*m
poly monic_transform(const poly& f);
