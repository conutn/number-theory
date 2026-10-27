#pragma once
// square root in Z[alpha] by p-adic Newton lifting (any degree; fallback for even degree)
//
// f may be non-monic: the work is done in Z[alpha] with alpha = lc*theta and the monic polynomial fhat (see monic_transform)
// given dependency S, gamma = fhat'(alpha)^2 prod(lc*a+b*alpha) is a square beta^2 in Z[alpha]   (lc = 1: fhat = f, alpha = theta)
// phi is a homomorphism theta -> m (mod n), phi(gamma) = f'(m)^2 prod(a+b*m), so
//  phi(beta)^2 = (f'(m) * sqrt(prod(a+b*m)))^2 (mod n)
// and gcd(phi(beta) +- f'(m) sqrt(prod(a+b*m)), n) is hopefully a nontrivial factor
//
// finding beta p-adically: sqrt of gamma in F_p[x]/f, p prime & f irred mod p
// then lift gamma^(-1/2) by Newton iteration
//  z <- z (3 - gamma z^2) / 2
// until beta = gamma*z satisfies beta^2 = gamma (in Z). the last step is CRT/Hensel reconstruction of the coefficients

#include "common.hpp"
#include "sieve.hpp"

// try one dependency
// returns true if nontrivial factor found, sets factor
bool try_dependency_newton(const vector<nfs_relation>& rels, const vector<size_t>& dep, const poly& f0, const bigint& m, const bigint& n, bigint& factor);
