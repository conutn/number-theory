// polynomials over Z and F_p: gcd, powers mod f, roots and factorization degrees mod p, irreducibility test

#include "polymod.hpp"
#include "numtheory.hpp"

poly reduce(poly a, const bigint& p)
{
	for (bigint& c : a.data()) c = modp(c, p);
	a.normalize();
	return a;
}

// scale polynomial so leading coefficient is 1
static void make_monic(poly& a, const bigint& p)
{
	auto& A = a.data();
	if (A.empty() || A.back() == 1) return;

	const bigint inv = inv_mod(A.back(), p);
	for (bigint& c : A) c = c * inv % p;
}

// a%=b over F_p
// requires:
//  a reduced
//  b reduced & monic
static void rem_monic(poly& a, const poly& b, const bigint& p)
{
	auto& A = a.data();
	const auto& B = b.data();
	if (B.empty()) throw domain_error("rem_monic: division by the zero polynomial");

	const size_t nb = B.size();

	while (A.size() >= nb)
	{
		const size_t shift = A.size() - nb;
		const bigint m = A.back();

		for (size_t i = 0; i + 1 < nb; ++i)
			A[shift + i] = modp(A[shift + i] - m * B[i], p);

		A.pop_back();
		a.normalize();
	}
}

poly poly_gcd(poly a, poly b, const bigint& p)
{
	a = reduce(std::move(a), p);
	b = reduce(std::move(b), p);

	while (!b.data().empty())
	{
		make_monic(b, p);
		rem_monic(a, b, p);
		std::swap(a, b);
	}

	make_monic(a, p);
	return a;
}

poly mul_mod(const poly& a, const poly& b, const poly& f, const bigint& p)
{
	const auto& A = a.data();
	const auto& B = b.data();
	if (A.empty() || B.empty()) return poly();

	vector<bigint> c(A.size() + B.size() - 1);
	for (size_t i = 0; i < A.size(); ++i)
		for (size_t j = 0; j < B.size(); ++j)
			c[i + j] += A[i] * B[j];

	poly r(std::move(c));
	for (bigint& x : r.data()) x %= p;
	r.normalize();
	rem_monic(r, f, p);
	return r;
}

// x^e mod (f, p)
// requires:
//  f reduced & monic & deg>=1
static poly x_pow_mod(const bigint& e, const poly& f, const bigint& p)
{
	poly r{ 1 };
	if (e == 0) return r;

	for (int i = (int)msb(e); i >= 0; --i)
	{
		r = mul_mod(r, r, f, p);

		if (bit_test(e, i) && !r.data().empty())
		{
			r.data().insert(r.data().begin(), bigint(0));
			rem_monic(r, f, p);
		}
	}

	return r;
}

// a(x-s) over F_p
// requires:
//  a reduced
static poly poly_shift(const poly& a, const bigint& s, const bigint& p)
{
	poly g(a);
	auto& c = g.data();
	const size_t n = c.size();
	if (n < 2) return g;

	const bigint t = modp(-s, p);

	for (size_t i = 0; i + 1 < n; ++i)
		for (size_t j = n - 1; j > i; --j)
			c[j - 1] = modp(c[j - 1] + t * c[j], p);

	return g;
}

// append to out roots of a
// requires:
//  a monic & squarefree & splits completely into linear factors
static void split_roots(const poly& a, const bigint& p, vector<bigint>& out)
{
	if (a.size() <= 1) return;

	if (a.size() == 2)
	{
		out.push_back(modp(-a[0], p));
		return;
	}

	const bigint power = (p - 1) / 2;

	for (bigint offset = 0; offset < p; ++offset)
	{
		poly g = poly_shift(a, offset, p);
		vector<bigint> local;

		if (g[0] == 0)
		{
			local.push_back(0);
			g.data().erase(g.data().begin());
			split_roots(g, p, local);
		}
		else
		{
			const poly r = x_pow_mod(power, g, p);
			if (r.size() <= 1) continue;

			poly plus(r), minus(r);
			plus.data()[0]  = modp(r[0] + 1, p);
			minus.data()[0] = modp(r[0] - 1, p);

			split_roots(poly_gcd(g, plus, p), p, local);
			split_roots(poly_gcd(g, minus, p), p, local);
		}

		for (const bigint& x : local) out.push_back(modp(x - offset, p));
		return;
	}

	throw logic_error("split_roots: polynomial does not split over F_p");
}

vector<bigint> roots(const poly& f, const bigint& p)
{
	if (p == 2)
	{
		vector<bigint> r;

		if (f(0) % 2 == 0) r.push_back(0);
		if (f(1) % 2 == 0) r.push_back(1);

		return r;
	}

	poly a = reduce(f, p);

	if (a.data().empty()) throw domain_error("roots: f is identically 0 mod p");
	if (a.size() == 1) return {};

	make_monic(a, p);

	poly h = x_pow_mod(p, a, p);
	{
		auto& H = h.data();
		if (H.size() < 2) H.resize(2, 0);
		H[1] = modp(H[1] - 1, p);
		h.normalize();
	}

	vector<bigint> out;
	split_roots(poly_gcd(a, h, p), p, out);

	sort(out.begin(), out.end());
	return out;
}

poly pow_mod(const poly& b, const bigint& e, const poly& f, const bigint& p)
{
	poly r{ 1 };
	if (e == 0) return r;
 
	for (int i = (int)msb(e); i >= 0; --i)
	{
		r = mul_mod(r, r, f, p);
		if (bit_test(e, i)) r = mul_mod(r, b, f, p);
	}
 
	return r;
}

poly derivative(const poly& f, const bigint& p)
{
	vector<bigint> d;
	for (size_t i = 1; i < f.size(); ++i) d.push_back(f[i] * i % p);
	return poly(std::move(d));
}

vector<size_t> factor_degrees(const poly& f, const bigint& p)
{
	const size_t d = f.size() - 1;
	vector<size_t> nk(d / 2 + 1, 0); // nk[k] = # irreducible factors of degree k
	poly h{ 0, 1 }; // x^(p^k) mod f 
	size_t found = 0;
 
	for (size_t k = 1; 2 * k <= d; ++k)
	{
		h = pow_mod(h, p, f, p);
 
		poly t(h); // h-x
		auto& T = t.data();
		if (T.size() < 2) T.resize(2, 0);
		T[1] = modp(T[1] - 1, p);
		t.normalize();
 
		// deg gcd(f, x^(p^k)-x) = sum of j*nk[j] over j | k
		const size_t e = poly_gcd(f, t, p).size() - 1;
		size_t known = 0;
		for (size_t j = 1; j < k; ++j)
			if (k % j == 0) known += j * nk[j];
 
		nk[k] = (e - known) / k;
		found += k * nk[k];
	}
 
	vector<size_t> degs;
	for (size_t k = 1; k < nk.size(); ++k)
		degs.insert(degs.end(), nk[k], k);
	if (found < d) degs.push_back(d - found);
	return degs;
}

irred is_irreducible(const poly& f0, size_t num_primes)
{
	const size_t d = f0.size() - 1;
	if (f0.size() < 2) throw domain_error("is_irreducible: constant polynomial");
	if (d == 1) return irred::proven;
	if (d > 63) throw domain_error("is_irreducible: degree too large");
 
	uint64_t possible = ((uint64_t(1) << d) - 1) & ~uint64_t(1); // candidate factor degrees
 
	size_t used = 0; // number of primes used
	for (unsigned q = 2; used < num_primes && q < 2000; ++q)
	{
		bool prime = true;
		for (unsigned r = 2; r * r <= q; ++r) if (q % r == 0) { prime = false; break; }
		if (!prime) continue;
 
		const bigint p = q;
		poly f = reduce(f0, p);
 
		if (f.size() != f0.size()) continue; // p | lc(f)
		make_monic(f, p);
		if (poly_gcd(f, derivative(f, p), p).size() != 1) continue; // not squarefree mod p
 
		++used;
 
		uint64_t sums = 1;
		for (size_t g : factor_degrees(f, p)) sums |= sums << g;
 
		possible &= sums;
		if (possible == 0) return irred::proven;
	}
 
	return irred::undetermined;
}

bigint eval_mod(const poly& f, const bigint& x, const bigint& p)
{
	bigint r = 0;
	for (size_t i = f.size(); i-- > 0;) r = modp(r * x + f[i], p);
	return r;
}

poly monic_transform(const poly& f)
{
	const size_t d = f.size() - 1;
	const bigint lc = f[d];
	if (lc == 1) return f;
	vector<bigint> c(d + 1);
	c[d] = 1;
	bigint pw = 1; // lc^(d-1-j), built from j = d-1 downwards
	for (size_t j = d; j-- > 0;)
	{
		c[j] = f[j] * pw;
		pw *= lc;
	}
	return poly(std::move(c));
}
