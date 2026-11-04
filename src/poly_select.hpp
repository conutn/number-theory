#pragma once
// GNFS polynomial selection: plain base-m polynomials + time-budgeted lattice (LLL) search for better, non-monic ones

// search runs until time budget is used up, budget is a fraction of predicted best time

// stages:
//   1. fixed grid of lattices (m = m0*k/2, a few skews), 8 best go to planner
//   2. random global search aimed at target lc L (10^3..10^7), skew follows from L
//     m is drawn just above smallest usable value (n/L)^(1/d)
//     followed by search around 'good' polynomials
//     translations f(x+t) & rotations f + w(x-m) & neighbouring skews
//   3. final planner evaluation of best candidates

#include "common.hpp"

struct gnfs_poly
{
	poly f; // base-m: monic, f(m) = n. lattice search: any lc, f(m) = 0 (mod n)
	bigint m;
};

// base-m polynomial selection
// requires:
//   n odd & composite & not perfect power
// d = 0 uses default degree for polynomial
gnfs_poly make_gnfs_poly(const bigint& n, unsigned d = 0, bool balanced = true, unsigned max_tries = 50);

// return polynomial sieve will use, best of lattice candidates and base-m
// always = false: if default polynomial is predicted under ~25 s search is skipped
gnfs_poly select_polynomial(const bigint& n, unsigned forced_d = 0, bool verbose = true, bool always = false);
