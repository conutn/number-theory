// small number theory helpers: modular arithmetic, 64-bit mulmod/powmod, Jacobi symbol, integer roots, prime generation

#include "numtheory.hpp"

bigint modp(const bigint& x, const bigint& p)
{
	bigint r = x % p;
	if (r < 0) r += p;
	return r;
}

bigint inv_mod(const bigint& a, const bigint& p)
{
	bigint r0 = p, r1 = modp(a, p), t0 = 0, t1 = 1;

	while (r1 != 0)
	{
		bigint q = r0 / r1;
		bigint r2 = r0 - q * r1;
		bigint t2 = t0 - q * t1;
		r0 = std::move(r1); r1 = std::move(r2);
		t0 = std::move(t1); t1 = std::move(t2);
	}

	if (r0 != 1) throw domain_error("inv_mod: element not invertible");
	return modp(t0, p);
}

bigint iroot(const bigint& n, unsigned d)
{
	if (n < 2) return n;
 
	bigint lo = 1, hi = bigint(1) << (msb(n) / d + 1); // hi^d > n
	while (lo < hi)
	{
		bigint mid = (lo + hi + 1) >> 1;
		if (pow(mid, d) <= n) lo = mid;
		else hi = mid - 1;
	}
 
	return lo;
}

double log_cpp(const bigint& x)
{
	const size_t bits = msb(x) + 1;
	if (bits <= 900) return std::log(x.convert_to<double>());
	const size_t sh = bits - 800;
	return std::log((x >> sh).convert_to<double>()) + double(sh) * 0.6931471805599453;
}

uint64_t lin_mod(int64_t a, uint64_t b, uint64_t p, uint64_t r)
{
	int64_t am = a % (int64_t)p;
	if (am < 0) am += (int64_t)p;
	return ((uint64_t)am + (b % p) * r) % p;
}

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
			if (r == 3 || r == 5) res = -res;
		}
		std::swap(a, n);
		if ((a & 3) == 3 && (n & 3) == 3) res = -res;
		a %= n;
	}
	return n == 1 ? res : 0;
}

u64 powmod64(u64 a, u64 e, u64 p)
{
	u64 r = 1 % p;
	a %= p;
	for (; e; e >>= 1, a = mulmod64(a, a, p)) if (e & 1) r = mulmod64(r, a, p);
	return r;
}

u64 mod_u64(const bigint& x, u64 p)
{
	bigint r = x % bigint(p);
	if (r < 0) r += p;
	return r.convert_to<u64>();
}
