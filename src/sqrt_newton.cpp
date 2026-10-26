// square root in Z[alpha] by p-adic Newton lifting (any degree; fallback for even degree)

#include "sqrt_newton.hpp"
#include "numtheory.hpp"
#include "polymod.hpp"

// a mod f over Z
// requires:
//  f monic
static void rem_monic_z(poly& a, const poly& f)
{
	auto& A = a.data();
	const auto& B = f.data();
	const size_t nb = B.size();

	while (A.size() >= nb)
	{
		const size_t shift = A.size() - nb;
		const bigint mm = A.back();
		for (size_t i = 0; i + 1 < nb; ++i) A[shift + i] -= mm * B[i];
		A.pop_back();
		a.normalize();
	}
}

// a*b mod f over Z
static poly mul_z(const poly& a, const poly& b, const poly& f)
{
	const auto& A = a.data();
	const auto& B = b.data();
	if (A.empty() || B.empty()) return poly();

	vector<bigint> c(A.size() + B.size() - 1);
	for (size_t i = 0; i < A.size(); ++i)
		for (size_t j = 0; j < B.size(); ++j)
			c[i + j] += A[i] * B[j];

	poly r(std::move(c));
	r.normalize();
	rem_monic_z(r, f);
	return r;
}

static bool is_one(const poly& a) { return a.size() == 1 && a[0] == 1; }

// sqrt in F_p[x]/f
// requires:
//  f irred mod p
//  g nonzero
// returns false if g is non-residue
static bool field_sqrt(const poly& g, const poly& f, const bigint& p, poly& root)
{
	const size_t d = f.size() - 1;
	const bigint qm1 = pow(p, (unsigned)d) - 1;

	if (!is_one(pow_mod(g, qm1 / 2, f, p))) return false;

	bigint t = qm1;
	unsigned s = 0;
	while ((t & 1) == 0) { t >>= 1; ++s; }

	poly z;
	for (bigint c = 0;; ++c) // z=x+c: find non-residue
	{
		z = poly{ c, bigint(1) };
		const poly e = pow_mod(z, qm1 / 2, f, p);
		if (e.size() == 1 && e[0] == p - 1) break;
	}

	poly c = pow_mod(z, t, f, p), r = pow_mod(g, (t + 1) / 2, f, p), tt = pow_mod(g, t, f, p);
	unsigned mm = s;

	while (!is_one(tt))
	{
		unsigned i = 0;
		poly x = tt;
		while (!is_one(x)) { x = mul_mod(x, x, f, p); ++i; if (i >= mm) return false; }

		poly b = c;
		for (unsigned j = 0; j + i + 1 < mm; ++j) b = mul_mod(b, b, f, p);

		r = mul_mod(r, b, f, p);
		c = mul_mod(b, b, f, p);
		tt = mul_mod(tt, c, f, p);
		mm = i;
	}

	root = r;
	return true;
}

// beta with beta^2=g in Z[x]/f
// returns false if g is not square
static bool poly_sqrt(const poly& g, const poly& f, poly& beta)
{
	const size_t d = f.size() - 1;
	bigint p = 1000;
	poly r;

	for (int tries = 0;; )
	{
		do ++p; while (!miller_rabin_test(p, 25));
		if (++tries > 500) throw runtime_error("poly_sqrt: no prime with f irreducible found");

		const poly fp = reduce(f, p);
		if (poly_gcd(fp, derivative(fp, p), p).size() != 1) continue;
		const vector<size_t> degs = factor_degrees(fp, p);
		if (degs.size() != 1 || degs[0] != d) continue;

		const poly gp = reduce(g, p);
		if (gp.data().empty()) continue;

		if (!field_sqrt(gp, fp, p, r)) return false; // g is not a square
		break;
	}

	// z = r^-1 = g^(-1/2) mod p
	poly z = pow_mod(r, pow(p, (unsigned)d) - 2, reduce(f, p), p);
	bigint M = p;

	for (int step = 0; step < 40; ++step)
	{
		// candidate beta = g z (symmetric residues mod M)
		beta = mul_mod(reduce(g, M), z, f, M);
		for (bigint& c : beta.data()) if (c > M / 2) c -= M;
		beta.normalize();
		if (mul_z(beta, beta, f).data() == g.data()) return true;

		// lift to M^2
		const bigint M2 = M * M;
		const poly gM = reduce(g, M2);
		poly t = mul_mod(gM, mul_mod(z, z, f, M2), f, M2); // g z^2
		vector<bigint> u(std::max<size_t>(t.size(), 1));
		for (size_t i = 0; i < t.size(); ++i) u[i] = modp(-t[i], M2);
		u[0] = modp(u[0] + 3, M2);
		poly up(std::move(u));
		up.normalize();

		z = mul_mod(z, up, f, M2);
		const bigint half = (M2 + 1) / 2;
		for (bigint& c : z.data()) c = c * half % M2;
		M = M2;
	}

	return false;
}

bool try_dependency_newton(const vector<nfs_relation>& rels, const vector<size_t>& dep, const poly& f0, const bigint& m, const bigint& n, bigint& factor)
{
	const bigint lc = f0.data().back();
	if (lc != 1 && (dep.size() & 1)) return false;

	const poly f = monic_transform(f0); // work in Z[alpha], alpha = lc*theta
	vector<bigint> dv;
	for (size_t i = 1; i < f.size(); ++i) dv.push_back(f[i] * i);
	const poly fd(dv);

	poly g = mul_z(fd, fd, f);
	bigint R = 1;

	for (size_t idx : dep)
	{
		const nfs_relation& r = rels[idx];
		g = mul_z(g, poly{ lc * bigint(r.a), bigint(r.b) }, f); // lc*a + b*alpha
		R *= bigint(r.a) + bigint(r.b) * m;
	}

	if (R <= 0) return false;
	const bigint y0 = boost::multiprecision::sqrt(R);
	if (y0 * y0 != R) return false; // rational side not a square

	poly beta;
	if (!poly_sqrt(g, f, beta)) return false;

	const bigint lcm = lc * m % n;
	bigint x = 0, fdm = 0;
	for (size_t i = beta.size(); i-- > 0;) x = modp(x * lcm + beta[i] % n, n);
	for (size_t i = fd.size(); i-- > 0;) fdm = modp(fdm * lcm + fd[i], n);
	bigint y = modp(y0 % n * fdm, n);
	if (lc != 1) y = y * boost::multiprecision::powm(lc, bigint(dep.size() / 2), n) % n;

	for (const bigint& c : { modp(x - y, n), modp(x + y, n) })
	{
		const bigint d = boost::multiprecision::gcd(c, n);
		if (d > 1 && d < n) { factor = d; return true; }
	}

	return false;
}
