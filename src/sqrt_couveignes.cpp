// Couveignes' square root (odd degree): word-size primes, sign fixed by the norm, CRT by product tree

#include "sqrt_couveignes.hpp"
#include "numtheory.hpp"
#include "polymod.hpp"
#include "crt_tree.hpp"

// arithmetic in F_p[x]/f, elements are d words, f given by its d low coefficients (f monic)
struct fp_ring
{
	size_t d;
	u64 p;
	vector<u64> fl; // f_0 .. f_{d-1} mod p

	vector<u64> mul(const vector<u64>& a, const vector<u64>& b) const
	{
		vector<u64> t(2 * d - 1, 0);
		for (size_t i = 0; i < d; ++i)
			for (size_t j = 0; j < d; ++j)
			{
				const u64 x = mulmod64(a[i], b[j], p);
				t[i + j] = t[i + j] + x >= p ? t[i + j] + x - p : t[i + j] + x; // p<2^62: no overflow
			}
		for (size_t i = 2 * d - 2; i >= d; --i)
		{
			const u64 c = t[i];
			if (!c) continue;
			for (size_t j = 0; j < d; ++j)
			{
				const u64 x = mulmod64(c, fl[j], p);
				t[i - d + j] = t[i - d + j] >= x ? t[i - d + j] - x : t[i - d + j] + p - x;
			}
		}
		t.resize(d);
		return t;
	}

	// g*(a+b*x)
	vector<u64> mul_linear(const vector<u64>& g, u64 a, u64 b) const
	{
		vector<u64> t(d + 1, 0);
		for (size_t i = 0; i < d; ++i)
		{
			t[i] = (t[i] + mulmod64(g[i], a, p)) % p;
			t[i + 1] = (t[i + 1] + mulmod64(g[i], b, p)) % p;
		}
		const u64 c = t[d];
		for (size_t j = 0; j < d; ++j)
		{
			const u64 x = mulmod64(c, fl[j], p);
			t[j] = t[j] >= x ? t[j] - x : t[j] + p - x;
		}
		t.resize(d);
		return t;
	}

	vector<u64> pow(vector<u64> b, bigint e) const
	{
		vector<u64> r(d, 0);
		r[0] = 1;
		for (int i = (int)msb(e); i >= 0; --i)
		{
			r = mul(r, r);
			if (bit_test(e, i)) r = mul(r, b);
		}
		return r;
	}

	// Norm_{F_p^d/F_p}(r) = det of the multiplication-by-r matrix
	u64 norm(const vector<u64>& r) const
	{
		vector<vector<u64>> M(d); // M[j] = r*x^j
		vector<u64> cur = r;
		for (size_t j = 0; j < d; ++j)
		{
			M[j] = cur;
			cur = mul_linear(cur, 0, 1);
		}

		u64 det = 1;
		for (size_t c = 0; c < d; ++c)
		{
			size_t piv = d;
			for (size_t i = c; i < d; ++i) if (M[i][c]) { piv = i; break; }
			if (piv == d) return 0;
			if (piv != c) { std::swap(M[piv], M[c]); det = (p - det) % p; }
			det = mulmod64(det, M[c][c], p);
			const u64 inv = powmod64(M[c][c], p - 2, p);
			for (size_t i = c + 1; i < d; ++i)
			{
				if (!M[i][c]) continue;
				const u64 fct = mulmod64(M[i][c], inv, p);
				for (size_t k = c; k < d; ++k)
				{
					const u64 x = mulmod64(fct, M[c][k], p);
					M[i][k] = M[i][k] >= x ? M[i][k] - x : M[i][k] + p - x;
				}
			}
		}
		return det;
	}
};

// |Res(f, f')| = |disc f| for monic f, by Bareiss on the Sylvester matrix
static bigint abs_disc(const poly& f)
{
	const size_t d = f.size() - 1, N = 2 * d - 1;
	vector<bigint> fd;
	for (size_t i = 1; i <= d; ++i) fd.push_back(f[i] * (unsigned)i);

	vector<vector<bigint>> M(N, vector<bigint>(N, 0));
	for (size_t i = 0; i + 1 < d; ++i)       // d-1 shifted copies of f (highest coefficient first)
		for (size_t k = 0; k <= d; ++k) M[i][i + k] = f[d - k];
	for (size_t i = 0; i < d; ++i)           // d shifted copies of f'
		for (size_t k = 0; k < d; ++k) M[d - 1 + i][i + k] = fd[d - 1 - k];

	bigint prev = 1;
	int sign = 1;
	for (size_t k = 0; k + 1 < N; ++k)
	{
		if (M[k][k] == 0)
		{
			size_t s = k + 1;
			while (s < N && M[s][k] == 0) ++s;
			if (s == N) return 0;
			std::swap(M[k], M[s]);
			sign = -sign;
		}
		for (size_t i = k + 1; i < N; ++i)
			for (size_t j = k + 1; j < N; ++j)
				M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) / prev;
		prev = M[k][k];
	}
	bigint det = M[N - 1][N - 1] * sign;
	return det < 0 ? -det : det;
}

// exponents of the primes in prod |a+b*m| over the dependency, from the factorizations stored with the relations
// (no trial division).  false if some total exponent or the number of negative values is odd
static bool rational_exponents(const vector<nfs_relation>& rels, const vector<size_t>& dep, const factor_bases& fb, std::map<u64, long>& ex)
{
	vector<long> acc(fb.rational.size(), 0);
	size_t negatives = 0;
	for (size_t idx : dep)
	{
		const nfs_relation& r = rels[idx];
		if (!r.cols.empty() && r.cols[0] == 0) ++negatives; // column 0 = sign
		for (const fac_entry& fe : r.rfac) acc[fe.idx] += fe.e;
	}
	if (negatives & 1) return false;
	for (size_t i = 0; i < acc.size(); ++i)
	{
		if (acc[i] & 1) return false;
		if (acc[i]) ex[fb.rational[i].p] += acc[i];
	}
	return true;
}

// same for prod |N(a+b*theta)|, summed per prime p (a prime with several roots in the base contributes v_p(N) once per
// ideal that actually divides a+b*theta, which for coprime a,b is at most one)
static bool algebraic_exponents(const vector<nfs_relation>& rels, const vector<size_t>& dep, const factor_bases& fb, std::map<u64, long>& ex)
{
	vector<long> acc(fb.algebraic.size(), 0);
	for (size_t idx : dep)
		for (const fac_entry& fe : rels[idx].afac) acc[fe.idx] += fe.e;

	for (size_t i = 0; i < acc.size(); ++i)
		if (acc[i]) ex[fb.algebraic[i].p] += acc[i];
	for (auto& kv : ex) if (kv.second & 1) return false;
	return true;
}

// approximate roots of monic f (Durand-Kerner)
static vector<std::complex<long double>> approx_roots(const poly& f)
{
	typedef std::complex<long double> C;
	const size_t d = f.size() - 1;
	vector<long double> c(d + 1);
	long double rad = 1;
	for (size_t k = 0; k <= d; ++k) { c[k] = f[k].convert_to<long double>(); if (k < d) rad = std::max(rad, 1 + std::fabs(c[k])); }

	vector<C> z(d);
	C w(0.4L, 0.9L), cur = 1;
	for (size_t k = 0; k < d; ++k) { z[k] = cur * rad * 0.5L; cur *= w; }

	for (int it = 0; it < 2000; ++it)
	{
		long double delta = 0;
		for (size_t i = 0; i < d; ++i)
		{
			C val = 1, den = 1;
			for (size_t k = d; k-- > 0;) val = val * z[i] + c[k];
			for (size_t j = 0; j < d; ++j) if (j != i) den *= z[i] - z[j];
			const C step = val / den;
			z[i] -= step;
			delta = std::max(delta, std::abs(step));
		}
		if (delta < 1e-16L * rad) break;
	}
	return z;
}

// upper bound on log2 of the coefficients of beta:
// beta = sum_j sigma_j(beta) f(x)/((x-theta_j) f'(theta_j)),  sigma_j(beta)/f'(theta_j) = sqrt(prod(a+b*theta_j))
static double beta_coeff_bits(const vector<nfs_relation>& rels, const vector<size_t>& dep, const poly& f, const bigint& lc)
{
	const size_t d = f.size() - 1;
	const auto roots_c = approx_roots(f);
	long double sumf = 0;
	for (size_t k = 0; k <= d; ++k) sumf += std::fabs(f[k].convert_to<long double>());

	const long double lcld = lc.convert_to<long double>();
	long double best = -1e300L;
	for (const auto& th : roots_c)
	{
		long double acc = 0;
		for (size_t idx : dep)
		{
			const nfs_relation& r = rels[idx];
			acc += std::log2(std::max(std::abs(std::complex<long double>((long double)r.a * lcld, 0) + (long double)r.b * th), (long double)1e-300L));
		}
		const long double t = 0.5L * acc + (long double)d * std::log2(std::max<long double>(1, std::abs(th)));
		best = std::max(best, t);
	}
	return (double)(best + std::log2((long double)d) + std::log2(sumf) + 8);
}

// next prime below `state` with p = 3 mod 4 and f irreducible mod p, p < 2^62
static u64 next_couveignes_prime(const poly& f, u64& state, const bigint& skip_dividing)
{
	const size_t d = f.size() - 1;
	for (;;)
	{
		state -= 4; // state = 3 (mod 4) stays 3 (mod 4)
		if (state < (u64(1) << 40)) throw runtime_error("couveignes: ran out of primes");
		const bigint P = state;
		if (!miller_rabin_test(P, 20)) continue;
		if (skip_dividing % P == 0) continue; // T = 0 (mod p) would make the sign ambiguous
		const poly fp = reduce(f, P);
		if (poly_gcd(fp, derivative(fp, P), P).size() != 1) continue;
		const vector<size_t> degs = factor_degrees(fp, P);
		if (degs.size() == 1 && degs[0] == d) return state;
	}
}

// sqrt(gamma) mod (f, p) with Norm = T mod p.  false if gamma is not a square mod p
static bool sqrt_mod_prime(const fp_ring& R, const poly& fd, const vector<nfs_relation>& rels, const vector<size_t>& dep, const bigint& lc, u64 T, vector<u64>& out)
{
	const size_t d = R.d;
	const u64 p = R.p;

	vector<u64> fdv(d, 0), g;
	for (size_t i = 0; i < fd.size() && i < d; ++i) fdv[i] = mod_u64(fd[i], p);
	g = R.mul(fdv, fdv);

	const u64 lcm = mod_u64(lc, p);
	for (size_t idx : dep)
	{
		const nfs_relation& r = rels[idx];
		const int64_t am = r.a % (int64_t)p;
		const u64 a = am < 0 ? (u64)(am + (int64_t)p) : (u64)am;
		g = R.mul_linear(g, mulmod64(lcm, a, p), (u64)r.b % p); // lc*a + b*alpha
	}

	bigint q = 1;
	for (size_t i = 0; i < d; ++i) q *= p;
	vector<u64> r = R.pow(g, (q + 1) / 4);
	if (R.mul(r, r) != g) return false;

	const u64 nr = R.norm(r);
	if (nr != T)
	{
		if (nr != (p - T) % p) return false; // neither sign matches: inconsistent
		for (u64& c : r) c = c ? p - c : 0;
	}
	out = r;
	return true;
}

// beta^2 = gamma (mod p)?  (used with primes not in the CRT)
static bool check_mod_prime(const fp_ring& R, const poly& fd, const vector<nfs_relation>& rels, const vector<size_t>& dep, const bigint& lc, const vector<bigint>& beta)
{
	const size_t d = R.d;
	const u64 p = R.p;
	vector<u64> fdv(d, 0), g, bp(d, 0);
	for (size_t i = 0; i < fd.size() && i < d; ++i) fdv[i] = mod_u64(fd[i], p);
	for (size_t i = 0; i < d && i < beta.size(); ++i) bp[i] = mod_u64(beta[i], p);
	g = R.mul(fdv, fdv);
	const u64 lcm = mod_u64(lc, p);
	for (size_t idx : dep)
	{
		const nfs_relation& r = rels[idx];
		const int64_t am = r.a % (int64_t)p;
		g = R.mul_linear(g, mulmod64(lcm, am < 0 ? (u64)(am + (int64_t)p) : (u64)am, p), (u64)r.b % p);
	}
	return R.mul(bp, bp) == g;
}

static fp_ring make_ring(const poly& f, u64 p)
{
	fp_ring R;
	R.d = f.size() - 1;
	R.p = p;
	for (size_t i = 0; i < R.d; ++i) R.fl.push_back(mod_u64(f[i], p));
	return R;
}

// beta (coefficients in the power basis) with beta^2 = gamma = f'(theta)^2 prod(a+b*theta), via Couveignes.
// disc_abs = |disc f|.  returns false if the dependency is not a square
// f here is the MONIC polynomial (fhat), lc the original leading coefficient: elements are lc*a + b*alpha
static bool couveignes_sqrt(const vector<nfs_relation>& rels, const vector<size_t>& dep, const poly& f, const factor_bases& fb, const bigint& disc_abs, const bigint& lc, vector<bigint>& beta)
{
	const size_t d = f.size() - 1;
	if (dep.size() & 1) return false; // parity column guarantees even size for lc != 1
	vector<bigint> dv;
	for (size_t i = 1; i < f.size(); ++i) dv.push_back(f[i] * (unsigned)i);
	const poly fd(dv);

	// T = |disc| * prod q^(e_q/2) = |Norm(beta)|, T^2 = Norm(gamma)
	std::map<u64, long> ex;
	if (!algebraic_exponents(rels, dep, fb, ex)) return false;

	for (size_t idx : dep)
		if (rels[idx].lp_a) ++ex[rels[idx].lp_a];
	for (const auto& kv : ex)
		if (kv.second & 1) return false;

	// Norm(lc*a + b*alpha) = lc^(d-1) * F(a,b)  =>  T = |disc fhat| * lc^((d-1)k/2) * prod q^(e_q/2)
	vector<bigint> tf{ disc_abs };
	if (lc != 1) tf.push_back(boost::multiprecision::pow(lc, (unsigned)((d - 1) * (dep.size() / 2))));
	for (auto& kv : ex) if (kv.second >= 2) tf.push_back(boost::multiprecision::pow(bigint(kv.first), (unsigned)(kv.second / 2)));
	const bigint T = balanced_product(std::move(tf));

	size_t need = (size_t)std::ceil((beta_coeff_bits(rels, dep, f, lc) + 2) / 61.0) + 1;

	u64 state = (u64(1) << 62) - 1; // = 3 (mod 4)
	vector<u64> primes;
	vector<vector<u64>> res;

	for (int round = 0; round < 8; ++round)
	{
		vector<u64> batch;
		while (primes.size() + batch.size() < need) batch.push_back(next_couveignes_prime(f, state, disc_abs)); // never divides disc
		const vector<u64> Tm = reduce_all(T, batch);

		for (size_t k = 0; k < batch.size(); ++k)
		{
			const fp_ring R = make_ring(f, batch[k]);
			vector<u64> r;
			if (Tm[k] == 0 || !sqrt_mod_prime(R, fd, rels, dep, lc, Tm[k], r)) return false; // not a square
			primes.push_back(batch[k]);
			res.push_back(std::move(r));
		}

		vector<bigint> X = crt_tree(primes, res, d);

		// verify with fresh primes
		bool ok = true;
		for (int k = 0; k < 2 && ok; ++k)
		{
			const u64 p = next_couveignes_prime(f, state, disc_abs);
			ok = check_mod_prime(make_ring(f, p), fd, rels, dep, lc, X);
		}

		if (ok) { beta = std::move(X); return true; }
		need += need / 4 + 2; // bound was too optimistic
	}

	return false;
}

bool try_dependency_couveignes(const vector<nfs_relation>& rels, const vector<size_t>& dep, const poly& f, const factor_bases& fb, const bigint& m, const bigint& n, bigint& factor)
{
	const bigint lc = f.data().back();
	if (lc != 1 && (dep.size() & 1)) return false;

	std::map<u64, long> rex;
	if (!rational_exponents(rels, dep, fb, rex)) return false;

	for (size_t idx : dep)
		if (rels[idx].lp_r) ++rex[rels[idx].lp_r];
	for (const auto& kv : rex)
		if (kv.second & 1) return false;

	const poly fh = monic_transform(f);
	vector<bigint> beta;
	if (!couveignes_sqrt(rels, dep, fh, fb, abs_disc(fh), lc, beta)) return false;

	// beta(lc*m)^2 = fhat'(lc*m)^2 prod(lc*(a+b*m)) = (fhat'(lc*m) * lc^(k/2) * sqrt(prod(a+b*m)))^2   (mod n)
	const bigint lcm = lc * m % n;
	bigint x = 0, fdm = 0, y = 1;
	for (size_t i = beta.size(); i-- > 0;) x = modp(x * lcm + beta[i], n);
	for (size_t i = fh.size(); i-- > 1;) fdm = modp(fdm * lcm + fh[i] * (unsigned)i, n);
	for (auto& kv : rex)
		y = y * boost::multiprecision::powm(bigint(kv.first), bigint(kv.second / 2), n) % n;
	y = y * fdm % n;
	if (lc != 1) y = y * boost::multiprecision::powm(lc, bigint(dep.size() / 2), n) % n;

	for (const bigint& c : { modp(x - y, n), modp(x + y, n) })
	{
		const bigint g = boost::multiprecision::gcd(c, n);
		if (g > 1 && g < n) { factor = g; return true; }
	}
	return false;
}
