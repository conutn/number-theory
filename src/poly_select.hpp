#pragma once
// GNFS polynomial selection: plain base-m polynomials + lattice (LLL) search for better, non-monic ones
//
// the rational polynomial is always g(x) = x - m. any f with f(m) = 0 (mod n) works, f need not be monic:
// the lattice { c : sum c_i m^i = 0 (mod n) } has dimension d+1 and determinant n, so LLL finds f with
// coefficients ~ n^(1/(d+1)) against ~ n^(1/d) for base-m. weighting coordinate i by skew^i makes the high
// coefficients small and the low ones large, which suits a sieve region with |a| >> b
//
// the leading coefficient is restricted to squarefree numbers with small prime factors (they are in the factor base as
// projective primes, other lc would put large primes into every relation with p | b). candidates are ranked by a cheap
// proxy (expected relations in a few lines, with Murphy's alpha) and the best are scored with the planner's predicted
// total time. the base-m polynomial is always among them, so the result is never worse than the default by the model

#include "common.hpp"

struct gnfs_poly
{
	poly f; // base-m: monic, f(m) = n. lattice search: any lc, f(m) = 0 (mod n)
	bigint m;
};

// base-m polynomial selection
// requires:
//  n odd & composite & not perfect power
// d=0 uses default degree for polynomial
gnfs_poly make_gnfs_poly(const bigint& n, unsigned d = 0, bool balanced = true, unsigned max_tries = 50);

// returns the polynomial the sieve will use: best of the lattice candidates and base-m
// forced_d = 0: try degrees d0-1 .. d0+1 around choose_degree(n) and keep the best by predicted time
// always = false: if the default (base-m) polynomial is predicted to finish in under ~25 s the search is skipped
gnfs_poly select_polynomial(const bigint& n, unsigned forced_d = 0, bool verbose = true, bool always = false);
