// square root in Z[alpha] by p-adic Newton lifting

#include "sqrt_newton.hpp"
#include "numtheory.hpp"
#include "polymod.hpp"

// reduces a mod f over Z in place
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
		for (size_t i = 0; i + 1 < nb; ++i) A[shift + i] -= mm * B[i]; // subtract mm*x^shift*f, the top term cancels
		A.pop_back();
		a.normalize(); // lower terms may have become zero
	}
}

// returns a*b mod f over Z
static poly mul_z(const poly& a, const poly& b, const poly& f)
{
	const auto& A = a.data();
	const auto& B = b.data();
	if (A.empty() || B.empty()) return poly();

	vector<bigint> c(A.size() + B.size() - 1);
	for (size_t i = 0; i < A.size(); ++i)
		for (size_t j = 0; j < B.size(); ++j)
			c[i + j] += A[i] * B[j];

	poly r(move(c));
	r.normalize();
	rem_monic_z(r, f);
	return r;
}

static bool is_one(const poly& a)
{
	return a.size() == 1 && a[0] == 1;
}

// sets root = sqrt(g) in F_p[x]/f (Tonelli-Shanks)
// requires:
//  f irred mod p
//  g nonzero
// returns false if g is non-residue
static bool field_sqrt(const poly& g, const poly& f, const bigint& p, poly& root)
{
	const size_t d = f.size() - 1;
	const bigint qm1 = pow(p, (unsigned)d) - 1; // q-1, q = p^d

	if (!is_one(pow_mod(g, qm1 / 2, f, p))) return false; // euler criterion

	bigint t = qm1; // q-1 = t * 2^s, t odd
	unsigned s = 0;
	while ((t & 1) == 0)
	{
		t >>= 1;
		++s;
	}

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
		// least i with tt^(2^i) = 1
		unsigned i = 0;
		poly x = tt;
		while (!is_one(x))
		{
			x = mul_mod(x, x, f, p);
			++i;
			if (i >= mm) return false;
		}

		poly b = c;
		for (unsigned j = 0; j + i + 1 < mm; ++j) b = mul_mod(b, b, f, p); // c^(2^(mm-i-1))

		r = mul_mod(r, b, f, p);
		c = mul_mod(b, b, f, p);
		tt = mul_mod(tt, c, f, p);
		mm = i;
	}

	root = r;
	return true;
}

// returns f'^2 * prod(lc*a + b*alpha) over the dependency, mod f over Z
// balanced product (a stack that merges equal sized products like a binary counter): the operands of every multiplication have about the same size,
// so a fast big integer multiplication pays off, while a running product (one small factor at a time) costs time quadratic in the dependency size
static poly dependency_product(const vector<nfs_relation>& rels, const vector<size_t>& dep, const poly& f, const poly& fd, const bigint& lc)
{
	vector<pair<poly, unsigned>> stack; // (product, level), a product of level l holds 2^l relations
	for (size_t idx : dep)
	{
		const nfs_relation& r = rels[idx];
		poly cur{ lc * bigint(r.a), bigint(r.b) };
		unsigned level = 0;
		while (!stack.empty() && stack.back().second == level)
		{
			cur = mul_z(stack.back().first, cur, f);
			stack.pop_back();
			++level;
		}
		stack.emplace_back(move(cur), level);
	}

	poly g = mul_z(fd, fd, f);
	while (!stack.empty()) // levels decrease towards the top, so g stays about as large as the next factor
	{
		g = mul_z(g, stack.back().first, f);
		stack.pop_back();
	}
	return g;
}

// sets y = sqrt(prod |a+b*m|) mod n from the factorizations stored with the relations, so the huge integer product is never built
// returns false if the product is not a square (odd exponent or odd number of negative values)
static bool rational_sqrt_mod_n(const vector<nfs_relation>& rels, const vector<size_t>& dep, const factor_bases& fb, const bigint& n, bigint& y)
{
	vector<long> acc(fb.rational.size(), 0); // exponent of every factor base prime
	map<u64, long> lp; // exponent of every large prime
	size_t negatives = 0;
	for (size_t idx : dep)
	{
		const nfs_relation& r = rels[idx];
		if (!r.cols.empty() && r.cols[0] == 0) ++negatives; // column 0 = sign
		for (const fac_entry& fe : r.rfac) acc[fe.idx] += fe.e;
		if (r.lp_r) ++lp[r.lp_r];
	}
	if (negatives & 1) return false;

	y = 1;
	for (size_t i = 0; i < acc.size(); ++i)
	{
		if (acc[i] & 1) return false;
		if (acc[i]) y = y * powm(bigint(fb.rational[i].p), bigint(acc[i] / 2), n) % n;
	}
	for (const auto& kv : lp)
	{
		if (kv.second & 1) return false;
		y = y * powm(bigint(kv.first), bigint(kv.second / 2), n) % n;
	}
	return true;
}

// sets beta with beta^2 = g in Z[x]/f
// returns false if g is not square
static bool poly_sqrt(const poly& g, const poly& f, poly& beta)
{
	const size_t d = f.size() - 1;
	bigint p = 1000;
	poly r;

	// find a prime with f irreducible mod p and g a square mod p
	for (int tries = 0;;)
	{
		do ++p; while (!miller_rabin_test(p, 25));
		if (++tries > 500) throw runtime_error("poly_sqrt: no prime with f irreducible found");

		const poly fp = reduce(f, p);
		if (poly_gcd(fp, derivative(fp, p), p).size() != 1) continue; // f must be squarefree mod p
		const vector<size_t> degs = factor_degrees(fp, p);
		if (degs.size() != 1 || degs[0] != d) continue;

		const poly gp = reduce(g, p);
		if (gp.data().empty()) continue;

		if (!field_sqrt(gp, fp, p, r)) return false; // g is not a square
		break;
	}

	// z = r^-1 = g^(-1/2) mod p
	poly z = pow_mod(r, pow(p, (unsigned)d) - 2, reduce(f, p), p);
	bigint M = p; // z is correct mod M, the precision doubles every step

	size_t gbits = 0; // beta has about half as many bits as g: nothing to test before M is that large
	for (const bigint& c : g.data()) gbits = std::max<size_t>(gbits, msb(c < 0 ? bigint(-c) : c) + 1);

	for (int step = 0; step < 40; ++step)
	{
		if (2 * (size_t)msb(M) + 64 >= gbits) // M is large enough for beta
		{
			// candidate beta = g z (symmetric residues mod M)
			beta = mul_mod(reduce(g, M), z, f, M);
			for (bigint& c : beta.data())
				if (c > M / 2) c -= M;
			beta.normalize();
			if (mul_z(beta, beta, f).data() == g.data()) return true;
		}

		// lift to M^2: z <- z (3 - g z^2) / 2
		const bigint M2 = M * M;
		const poly gM = reduce(g, M2);
		poly t = mul_mod(gM, mul_mod(z, z, f, M2), f, M2); // g z^2
		vector<bigint> u(std::max<size_t>(t.size(), 1)); // u = 3 - t
		for (size_t i = 0; i < t.size(); ++i) u[i] = modp(-t[i], M2);
		u[0] = modp(u[0] + 3, M2);
		poly up(move(u));
		up.normalize();

		z = mul_mod(z, up, f, M2);
		const bigint half = (M2 + 1) / 2; // 1/2 mod M2 (M2 odd)
		for (bigint& c : z.data()) c = c * half % M2;
		M = M2;
	}

	return false;
}

bool try_dependency_newton(const vector<nfs_relation>& rels, const vector<size_t>& dep, const poly& f0, const factor_bases& fb, const bigint& m, const bigint& n, bigint& factor)
{
	const bigint lc = f0.data().back();
	if (lc != 1 && (dep.size() & 1)) return false;

	const poly f = monic_transform(f0); // work in Z[alpha], alpha = lc*theta
	vector<bigint> dv;
	for (size_t i = 1; i < f.size(); ++i) dv.push_back(f[i] * i);
	const poly fd(dv); // f'

	// rational side first, it is cheap: y0 = sqrt(prod(a + b*m)) mod n
	bigint y0;
	if (!rational_sqrt_mod_n(rels, dep, fb, n, y0)) return false; // rational side not a square

	// gamma = f'^2 prod(lc*a + b*alpha)
	const poly g = dependency_product(rels, dep, f, fd, lc);

	poly beta;
	if (!poly_sqrt(g, f, beta)) return false;

	// x = beta(lc*m), y = f'(lc*m) * lc^(k/2) * sqrt(R)   (mod n)
	const bigint lcm = lc * m % n;
	bigint x = 0, fdm = 0;
	for (size_t i = beta.size(); i-- > 0;) x = modp(x * lcm + beta[i] % n, n);
	for (size_t i = fd.size(); i-- > 0;) fdm = modp(fdm * lcm + fd[i], n);
	bigint y = modp(y0 * fdm, n);
	if (lc != 1) y = y * powm(lc, bigint(dep.size() / 2), n) % n;

	for (const bigint& c : { modp(x - y, n), modp(x + y, n) })
	{
		const bigint d = boost::multiprecision::gcd(c, n);
		if (d > 1 && d < n)
		{
			factor = d;
			return true;
		}
	}

	return false;
}
