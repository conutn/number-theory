#pragma once
// Couveignes' square root (odd degree): word-size primes, sign fixed by the norm, CRT by product tree
//
// beta^2 = gamma = f'(theta)^2 prod(a+b*theta) in Z[theta], f monic of ODD degree d
// instead of lifting gamma^(-1/2) p-adically (every step works with the full-size gamma, whose coefficients have
// millions of bits for big n), work modulo many word-size primes p_i with f irreducible mod p_i:
//   1. gamma mod (f, p_i) is built relation by relation in O(d) per relation, 64-bit arithmetic
//   2. r_i = sqrt(gamma) in F_{p_i}[x]/f. with p_i = 3 (mod 4) and d odd, p_i^d = 3 (mod 4), so
//      r_i = gamma^((p_i^d+1)/4): a single exponentiation, no Tonelli-Shanks
//   3. sign: d odd => Norm(-r) = -Norm(r), so exactly one of +-r has Norm = T (mod p_i), where
//      T = |disc f| * prod q^(e_q/2) is a fixed integer with T^2 = Norm(gamma) (e_q = exponent of q in
//      prod Norm(a+b*theta), all even for a dependency). choosing the sign by the norm makes all the r_i
//      residues of the SAME beta (rather than a random +-beta each)
//   4. coefficients of beta from CRT over the p_i. the bound on their size comes from the embeddings
//      sigma_j(beta) = f'(theta_j) sqrt(prod(a+b*theta_j)), so no big arithmetic is needed to know how many
//      primes to use. two spare primes check beta^2 = gamma, more primes are added if the check fails
// the rational side sqrt(prod(a+b*m)) mod n is taken from prime exponents (no huge product)
//
// for even d the sign trick does not work (Norm(-r) = Norm(r)), try_dependency falls back to the
// p-adic Newton method (try_dependency_newton) in that case

#include "common.hpp"
#include "factor_base.hpp"
#include "sieve.hpp"

// try one dependency (Couveignes version), on success factor is a nontrivial factor of n
bool try_dependency_couveignes(const vector<nfs_relation>& rels, const vector<size_t>& dep, const poly& f, const factor_bases& fb, const bigint& m, const bigint& n, bigint& factor);
