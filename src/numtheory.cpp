// small number theory helpers: modular arithmetic, 64-bit mulmod/powmod, Jacobi symbol, integer roots, prime generation

#include "numtheory.hpp"

// x mod p -> [0, p)
bigint modp(const bigint& x, const bigint& p)
{
	bigint r = x % p;
	if (r < 0) r += p; // % keeps the sign of x
	return r;
}

// a^-1 mod p by extended euclidean
// throws domain_error if gcd(a,p) != 1
bigint inv_mod(const bigint& a, const bigint& p)
{
	bigint r0 = p, r1 = modp(a, p), t0 = 0, t1 = 1; // invariant: r_i = t_i * a (mod p)

	while (r1 != 0)
	{
		bigint q = r0 / r1;
		bigint r2 = r0 - q * r1;
		bigint t2 = t0 - q * t1;
		r0 = std::move(r1); r1 = std::move(r2);
		t0 = std::move(t1); t1 = std::move(t2);
	}

	if (r0 != 1) throw domain_error("inv_mod: element not invertible"); // r0 = gcd(a,p)
	return modp(t0, p);
}

// floor(n^(1/d))
bigint iroot(const bigint& n, unsigned d)
{
	if (n < 2) return n;

	bigint lo = 1, hi = bigint(1) << (msb(n) / d + 1); // hi^d > n
	while (lo < hi)
	{
		bigint mid = (lo + hi + 1) >> 1; // round up so lo = mid always makes progress
		if (pow(mid, d) <= n) lo = mid;
		else hi = mid - 1;
	}

	return lo;
}

// log of x>0
// works for x too large to fit in a double
double log_cpp(const bigint& x)
{
	const size_t bits = msb(x) + 1;
	if (bits <= 900) return std::log(x.convert_to<double>()); // double overflows near 2^1024
	const size_t sh = bits - 800;
	return std::log((x >> sh).convert_to<double>()) + double(sh) * 0.6931471805599453; // log(x) ~ log(x >> sh) + sh*log(2)
}

// a+b*r mod p
// requires:
//  b>=0
//  p<2^32
uint64_t lin_mod(int64_t a, uint64_t b, uint64_t p, uint64_t r)
{
	int64_t am = a % (int64_t)p;
	if (am < 0) am += (int64_t)p;
	return ((uint64_t)am + (b % p) * r) % p; // p<2^32 so (b%p)*r fits in 64 bits (r<p)
}

// Jacobi symbol (a/n)
// requires:
//  n odd & positive
// note: equal to Legendre symbol if n is prime
int jacobi(uint64_t a, uint64_t n)
{
	int res = 1;
	a %= n;
	while (a)
	{
		while (!(a & 1))
		{
			a >>= 1;
			const uint64_t r = n & 7;
			if (r == 3 || r == 5) res = -res; // (2/n) = -1 iff n = 3,5 mod 8
		}
		std::swap(a, n);
		if ((a & 3) == 3 && (n & 3) == 3) res = -res; // quadratic reciprocity
		a %= n;
	}
	return n == 1 ? res : 0; // n != 1 means gcd(a,n) > 1
}

// a^e mod p
u64 powmod64(u64 a, u64 e, u64 p)
{
	u64 r = 1 % p; // handles p = 1
	a %= p;
	for (; e; e >>= 1, a = mulmod64(a, a, p)) if (e & 1) r = mulmod64(r, a, p);
	return r;
}

// x mod p -> [0, p)
u64 mod_u64(const bigint& x, u64 p)
{
	bigint r = x % bigint(p);
	if (r < 0) r += p;
	return r.convert_to<u64>();
}
