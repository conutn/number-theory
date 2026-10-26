#pragma once
// GNFS polynomial selection: plain base-m polynomials and the lattice (LLL) search for better, non-monic ones
//
// The rational polynomial stays g(x) = x - m.  Any f with f(m) = 0 (mod n) works, f need not be monic:
// the lattice L = { (c_0..c_d) : sum c_i m^i = 0 (mod n) } has dimension d+1 and determinant n, so LLL finds
// f with coefficients ~ n^(1/(d+1)), against ~ n^(1/d) for base-m (f(m) = n, monic).  Weighting coordinate i
// by skew^i makes the high coefficients small and the low ones large, which suits a sieve region with |a| >> b.
//
// the leading coefficient lc is restricted to squarefree numbers with small prime factors (they are in the factor base as
// projective primes; other lc would put large primes into every relation with p | b).  Candidates (several m, several skews,
// small combinations of the reduced vectors) are ranked by a cheap proxy (expected relations in a few lines, with
// Murphy's alpha) and the best few are scored with the planner's predicted total time; the plain base-m polynomial is always
// among them, so the result is never worse than the default according to the model.

#include "common.hpp"

struct gnfs_poly
{
	poly f; // base-m: monic, f(m) = n.  lattice search: any lc, f(m) = 0 (mod n).  the rational polynomial is always g(x) = x - m
	bigint m;
};

// base-m polynomial selection
// requires:
//  n odd & composite & not perfect power
// d=0 uses default degree for polynomial
gnfs_poly make_gnfs_poly(const bigint& n, unsigned d = 0, bool balanced = true, unsigned max_tries = 50);

// the polynomial the sieve will use: searches the lattice candidates (see above) and falls back to base-m.
// forced_d = 0: try degrees d0-1 .. d0+1 around choose_degree(n) and keep the best by predicted time
// always = false: if the default (base-m) polynomial is already predicted to finish in under ~25 s the search is skipped
gnfs_poly select_polynomial(const bigint& n, unsigned forced_d = 0, bool verbose = true, bool always = false);
