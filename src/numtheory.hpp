#pragma once
// small number theory helpers: modular arithmetic, 64-bit mulmod/powmod, Jacobi symbol, integer roots, prime generation

#include "common.hpp"

// generate primes up to lim (inclusive)
template <typename T = uint64_t> // using bigint is much slower and unnecessary most of the time
static vector<T> gen_prime(uint64_t lim)
{
	vector<T> primes;
	if (lim < 2) return primes;
	primes.push_back(T(2));
 
	const uint64_t half = (lim - 1) / 2; // only consider odds
	vector<bool> composite(half + 1, false);
 
	for (uint64_t j = 1; j <= half; ++j)
	{
		if (composite[j]) continue;
 
		const uint64_t p = 2 * j + 1;
		primes.push_back(T(p));
 
		// mask out from p^2 to lim
		if (p <= lim / p)
			for (uint64_t k = (p * p) / 2; k <= half; k += p) composite[k] = true;
	}
 
	return primes;
}

// x mod p -> [0, p)
bigint modp(const bigint& x, const bigint& p);

// a^-1 mod p by extended euclidean 
bigint inv_mod(const bigint& a, const bigint& p);

// floor(n^(1/d))
bigint iroot(const bigint& n, unsigned d);

// log of x>0
double log_cpp(const bigint& x);

// a+b*r mod p
// requires:
//  b>=0
//  p<2^32
uint64_t lin_mod(int64_t a, uint64_t b, uint64_t p, uint64_t r);

// Jacobi symbol (a/n)
// requires:
//  n odd & positive
// note: equal to Legendre symbol if n is prime
int jacobi(uint64_t a, uint64_t n);

inline int ctz64(uint64_t x)
{
	return __builtin_ctzll(x);
}

inline u64 mulmod64(u64 a, u64 b, u64 p) { return (u64)((u128)a * b % p); }

u64 powmod64(u64 a, u64 e, u64 p);

u64 mod_u64(const bigint& x, u64 p);
