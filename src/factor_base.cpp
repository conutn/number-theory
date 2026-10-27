// factor bases (rational, algebraic incl. projective roots, quadratic characters), bound heuristics, norms

#include "factor_base.hpp"
#include "numtheory.hpp"
#include "polymod.hpp"

// returns smoothness bound from the size of n
uint64_t choose_bound(const bigint& n)
{
	// ln(bound) = 0.87 * (ln n)^(1/3) (ln ln n)^(2/3) + 0.7, L[1/3] shape fitted to degree-3 runs on 12..32 digit semiprimes (A = max(10000, 4*bound))
	// smallest bound with enough relations: ~800 (12-16 digits), ~3200 (20), ~6400 (24), ~12800 (28), ~25600 (32)
	// the offset puts the result slightly above that threshold (fastest runs are at or one doubling above); main() doubles it if still too small
	// NOT measured beyond ~32 digits: treat larger values as a starting point
	const double ln = (double)(msb(n) + 1) * 0.693147; // ln n from bit length
	const double lnln = std::log(ln);
	const double b = std::exp(0.87 * std::cbrt(ln) * std::pow(lnln, 2.0 / 3.0) + 0.7);

	return (uint64_t)std::min(std::max(b, 200.0), 4.0e9);
}

// returns number of quadratic characters
size_t choose_num_qc(const bigint& n)
{
	return 20 + msb(n) / 8;
}

// build all three bases
factor_bases make_factor_bases(const poly& f, const bigint& m, uint64_t rat_bound, uint64_t alg_bound, size_t num_qc)
{
	factor_bases fb;

	for (uint64_t p : gen_prime(rat_bound))
		fb.rational.push_back({ p, (m % p).convert_to<uint64_t>() });

	const bigint lc = f.data().back();
	fb.parity = (lc != 1);

	for (uint64_t p : gen_prime(alg_bound))
	{
		const bigint P = p;

		if (reduce(f, P).data().empty())
			throw runtime_error("make_factor_bases: " + to_string(p) + " divides every coefficient of f, so it divides n");

		if (lc % P == 0) fb.algebraic.push_back({ p, FB_PROJECTIVE }); // root at infinity

		for (const bigint& r : roots(f, P))
			fb.algebraic.push_back({ p, r.convert_to<uint64_t>() });
	}

	bigint q = std::max<uint64_t>(std::max(rat_bound, alg_bound), 2); // characters use primes above both bounds
	const bigint give_up = q + 100000000;

	while (fb.quadratic.size() < num_qc)
	{
		if (++q > give_up) throw runtime_error("make_factor_bases: could not find enough quadratic-character primes");
		if (!miller_rabin_test(q, 25)) continue;

		const poly fq = reduce(f, q);
		if (fq.size() != f.size()) continue; // q | lc(f)
		const poly dq = derivative(fq, q);

		for (const bigint& s : roots(f, q))
		{
			if (fb.quadratic.size() >= num_qc) break;
			if (eval_mod(dq, s, q) == 0) continue; // repeated root
			fb.quadratic.push_back({ q.convert_to<uint64_t>(), s.convert_to<uint64_t>() });
		}
	}

	return fb;
}

// N(a+b*theta) = (-b)^d f(-a/b), d = deg f
bigint norm(poly f, bigint a, bigint b)
{
	vector<bigint> a_pow(f.size()), b_pow(f.size()); // a_pow[i] = a^i, b_pow[i] = b^(d-i)
	a_pow[0] = 1, b_pow[b_pow.size() - 1] = 1;

	for (size_t i = 1; i < a_pow.size(); ++i)
	{
		a_pow[i] = a_pow[i - 1] * a;
		b_pow[b_pow.size() - i - 1] = b_pow[b_pow.size() - i] * b;
	}

	bigint res = 0;

	for (size_t n = 0; n < f.size(); ++n)
	{
		bigint x = f[n] * a_pow[n] * b_pow[n];

		if ((f.size() + n) & 1) res += x; // sign (-1)^(d+n)
		else res -= x;
	}

	return res;
}

// number of matrix columns: sign + all bases (+ parity)
size_t num_columns(const factor_bases& fb)
{
	return 1 + fb.rational.size() + fb.algebraic.size() + fb.quadratic.size() + (fb.parity ? 1 : 0);
}
