// GNFS polynomial selection: plain base-m polynomials + time-budgeted lattice (LLL) search for better, non-monic ones

// search runs until time budget is used up, budget is a fraction of predicted best time

// stages:
//   1. fixed grid of lattices (m = m0*k/2, skews from a grid that grows with n), 8 best go to planner
//   2. random global search aimed at target lc L (10^3..10^7), skew follows from L
//     m is drawn just above smallest usable value (n/L)^(1/d)
//     followed by search around 'good' polynomials
//     translations f(x+t) & rotations f + w(x-m) & neighbouring skews
//   3. final planner evaluation of best candidates

#include <cstdio>
#include "poly_select.hpp"
#include "numtheory.hpp"
#include "polymod.hpp"
#include "factor_base.hpp"
#include "planner.hpp"

// rules of thumb for d given n
static unsigned choose_degree(const bigint& n)
{
	static const struct { unsigned max_digits; unsigned d; } table[] = {
		{  45, 3 },
		{  95, 4 },
		{ 220, 5 },
	};

	const unsigned digits = (unsigned)(msb(n) * 0.30103) + 1; // log10(2) * bits

	for (const auto& row : table)
		if (digits <= row.max_digits) return row.d;

	return 6;
}

// n written in base m, coefficients c[i] with n = sum c[i] m^i
// digits in (-m/2, m/2] if balanced is true
// digits in [0, m) otherwise
static poly base_m_poly(bigint n, const bigint& m, bool balanced)
{
	vector<bigint> c;
	const bigint half = m / 2;

	while (n != 0)
	{
		bigint r = n % m;
		if (balanced && r > half) r -= m;
		c.push_back(r);
		n = (n - r) / m;
	}

	return poly(move(c));
}

gnfs_poly make_gnfs_poly(const bigint& n, unsigned d, bool balanced, unsigned max_tries)
{
	if (d == 0) d = choose_degree(n);

	const bigint m0 = iroot(n, d);
	if (m0 <= d + 1) throw domain_error("make_gnfs_poly: n is too small for this degree");

	for (unsigned k = 0; k < max_tries; ++k)
	{
		const bigint m = m0 - k;
		if (m <= d + 1) break;

		poly f = base_m_poly(n, m, balanced);
		if (f.size() != d + 1) f = base_m_poly(n, m, false); // balanced digits can change the length
		if (f.size() != d + 1) continue;

		if (f(m) != n) throw logic_error("make_gnfs_poly: f(m) != n");
		if (is_irreducible(f) == irred::proven) return { f, m };
	}

	throw runtime_error("make_gnfs_poly: no irreducible polynomial found");
}

#ifdef GNFS_NO_GMP
typedef cpp_rational bigrat;
#else
typedef mpq_rational bigrat;
#endif

namespace polysel
{
	typedef chrono::steady_clock clk;

	// value of env variable name, or def if unset | not a number
	static double env_double(const char* name, double def)
	{
		const char* s = getenv(name);
		if (!s || !*s) return def;
		char* end = nullptr;
		const double v = strtod(s, &end);
		return end == s ? def : v;
	}

	// round q (half rounds up)
	static bigint round_rat(const bigrat& q)
	{
		const bigint num = numerator(q), den = denominator(q); // den > 0
		const bigint a = 2 * num + den, b = 2 * den;
		bigint r = a / b;
		if (a % b != 0 && a < 0) --r; // division truncates towards zero
		return r;
	}

	// return LLL reduced basis (delta = 0.99) of rows of B
	// exact rational Gram-Schmidt: ok for dimension <= ~8
	static vector<vector<bigint>> lll(vector<vector<bigint>> B)
	{
		const size_t n = B.size(), dim = B[0].size();
		vector<vector<bigrat>> Bs(n, vector<bigrat>(dim)), mu(n, vector<bigrat>(n, bigrat(0))); // Bs = ortho rows, mu = GS coefficients
		vector<bigrat> nb(n); // |Bs[i]|^2

		auto gso = [&]()
		{
			for (size_t i = 0; i < n; ++i)
			{
				vector<bigrat> v(dim);
				for (size_t c = 0; c < dim; ++c) v[c] = bigrat(B[i][c]);
				for (size_t j = 0; j < i; ++j)
				{
					bigrat dp = 0;
					for (size_t c = 0; c < dim; ++c) dp += bigrat(B[i][c]) * Bs[j][c];
					mu[i][j] = (nb[j] == 0) ? bigrat(0) : bigrat(dp / nb[j]);
					for (size_t c = 0; c < dim; ++c) v[c] -= mu[i][j] * Bs[j][c];
				}
				nb[i] = 0;
				for (size_t c = 0; c < dim; ++c) nb[i] += v[c] * v[c];
				Bs[i] = move(v);
			}
		};

		gso();
		const bigrat delta(99, 100);
		size_t k = 1;
		while (k < n)
		{
			// size reduction
			for (size_t j = k; j-- > 0;)
			{
				const bigint q = round_rat(mu[k][j]);
				if (q != 0)
				{
					for (size_t c = 0; c < dim; ++c) B[k][c] -= q * B[j][c];
					gso();
				}
			}

			// Lovasz condition
			if (nb[k] >= (delta - mu[k][k - 1] * mu[k][k - 1]) * nb[k - 1]) ++k;
			else
			{
				std::swap(B[k], B[k - 1]);
				gso();
				k = std::max<size_t>(k - 1, 1);
			}
		}
		return B;
	}

	// limits for the leading coefficient of a candidate, overridable by env variables
	// GNFS_LC_MAX: largest lc, GNFS_LC_PRIME_MAX: largest prime factor of lc
	// GNFS_LC_POWER_PRIME_MAX: primes up to this may divide lc more than once, larger primes only once
	// GNFS_LC_POWER_MAX: largest multiplicity such a prime may have (published polynomials use lc with factors like 2^5 * 3^2)
	// for n of ~100 digits only a handful of lattice vectors have lc <= 1e7: the best one found for RSA-100 is 332810 = 2*5*23*1447,
	// so a prime bound of 1000 left no candidate at all and 10000 was needed; check how factor_base treats primes dividing lc before raising the prime bound further
	static double lc_max() { static const double v = env_double("GNFS_LC_MAX", 1e9); return v; }
	static double lc_prime_max() { static const double v = env_double("GNFS_LC_PRIME_MAX", 10000); return v; }
	static unsigned lc_power_prime_max() { static const unsigned v = (unsigned)env_double("GNFS_LC_POWER_PRIME_MAX", 100); return v; }
	static unsigned lc_power_max() { static const unsigned v = (unsigned)env_double("GNFS_LC_POWER_MAX", 6); return v; }

	// primes up to sqrt(lc_max), enough to factor any allowed lc by trial division
	static const vector<unsigned>& lc_primes()
	{
		static const vector<unsigned> ps = planner::small_primes((unsigned)sqrt(lc_max()) + 2);
		return ps;
	}

	// true if c in (0, lc_max()], all prime factors <= lc_prime_max(), and every prime divides c once
	// except primes <= lc_power_prime_max(), which may divide it up to lc_power_max() times
	static bool good_leading_coefficient(bigint c)
	{
		if (c <= 0 || c > bigint((unsigned long long)lc_max())) return false;
		uint64_t n = c.convert_to<uint64_t>();
		for (unsigned p : lc_primes())
		{
			if ((uint64_t)p * p > n) break; // what is left is 1 or a prime
			if (n % p != 0) continue;
			if (p > lc_prime_max()) return false;
			unsigned k = 0; // multiplicity of p in lc
			while (n % p == 0)
			{
				n /= p;
				++k;
			}
			if (k > 1 && (p > lc_power_prime_max() || k > lc_power_max())) return false; // repeated large prime, or too high a power
		}
		return n == 1 || n <= lc_prime_max(); // a leftover n > 1 is a prime that occurs once
	}

	// skews of the stage 1 grid: 2 per decade from 1 up to where the leading coefficient of a reduced vector falls to ~1e3
	// reduced vectors have lc ~ n^(1/(d+1)) / skew^(d/2), the grid never ends below 1e4 (the old grid) or above 4e9
	// input: ln n and degree d; output: skews in increasing order
	static vector<unsigned> stage1_skews(double ln_n, unsigned d)
	{
		const double ln_top = (2.0 / d) * (ln_n / (d + 1) - log(1e3));
		const double top = std::min(4e9, std::max(1e4, exp(ln_top)));
		vector<unsigned> out;
		for (int i = 0; i <= 40; ++i)
		{
			const double s = pow(10.0, 0.5 * i);
			if (s > top * 1.0001) break;
			const unsigned u = (unsigned)llround(s);
			if (out.empty() || u != out.back()) out.push_back(u);
		}
		return out;
	}

	// predicted seconds as text: plain integer below 1e15, scientific notation above, "infeasible" for the 1e300 sentinel
	// input: seconds as a double; output: string
	static string fmt_seconds(double s)
	{
		if (s >= 1e299) return "infeasible";
		char buf[32];
		if (s < 1e15) snprintf(buf, sizeof buf, "%.0f", s);
		else snprintf(buf, sizeof buf, "%.3g", s); // casting a double above 9.2e18 to long long is undefined, it printed garbage before
		return buf;
	}

	struct candidate
	{
		poly f;
		bigint m;
		double proxy = 0;
		double seconds = 0; // predicted total time, 0 = not evaluated
		unsigned d = 0;
		bool tried = false; // already given to planner
		bool proven = false; // irred already proven
		unsigned skew = 1; // skew of lattice it came from
	};

	// divides the integer content (gcd of all coefficients) out of a coefficient vector, in place
	// input: coefficient vector c; output: the content that was removed (1 if c was already primitive)
	// a lattice vector can be a multiple g * f' of another lattice vector, g | every coefficient
	// f' has the same roots and number field as f, and f'(m) = 0 (mod n) still holds since gcd(g, n) = 1
	// (if not, the caller's modp(f(m), n) check rejects the candidate)
	static bigint remove_content(vector<bigint>& c)
	{
		bigint g = 0;
		for (const bigint& x : c)
		{
			g = boost::multiprecision::gcd(g, x); // gcd(0, x) = |x|, sign of coefficients does not matter
			if (g == 1) return g; // already primitive
		}
		if (g > 1)
			for (bigint& x : c) x /= g; // exact division, g divides every coefficient
		return g;
	}

	// polynomials from the lattice { c : sum c_i m^i = 0 (mod n) } w/ coordinate i weighted by skew^i
	// first use_max reduced vectors, combined with coefficients in [-range, range]
	// combinations are rejected:
	//   by leading coefficient
	//   when skewed norm exceeds slack * norm(shortest reduced vector)
	static void lattice_candidates(const bigint& n, const bigint& m, unsigned d, unsigned skew, size_t use_max, int range, double slack, vector<poly>& out)
	{
		const size_t dim = d + 1;
		vector<bigint> w(dim); // w[i] = skew^i
		w[0] = 1;
		for (size_t i = 1; i < dim; ++i) w[i] = w[i - 1] * skew;

		// basis: (n, 0, ..., 0) and (-m^i mod n, 0, .., 1, .., 0) for i >= 1
		vector<vector<bigint>> rows(dim, vector<bigint>(dim, 0));
		rows[0][0] = n;
		bigint mp = 1;
		for (size_t i = 1; i < dim; ++i)
		{
			mp = mp * m % n;
			rows[i][0] = modp(-mp, n);
			rows[i][i] = 1;
		}
		for (auto& r : rows)
			for (size_t i = 0; i < dim; ++i) r[i] *= w[i];

		const auto red = lll(rows);
		vector<vector<bigint>> v(dim, vector<bigint>(dim)); // reduced vectors with the weights divided out again
		for (size_t r = 0; r < dim; ++r)
			for (size_t i = 0; i < dim; ++i) v[r][i] = red[r][i] / w[i];

		vector<double> wd(dim);
		for (size_t i = 0; i < dim; ++i) wd[i] = (double)w[i];
		const auto skewed_norm = [&](const vector<bigint>& c)
		{
			double s = 0;
			for (size_t i = 0; i < dim; ++i)
			{
				const double x = (double)c[i] * wd[i];
				s += x * x;
			}
			return sqrt(s);
		};
		// the limit is relative to the shortest reduced vector that has a leading coefficient
		// for large n the shortest reduced vectors are multiples of (x - m): lc = 0, norm far below every useful polynomial, a limit relative to them rejects everything
		double ref = 0;
		for (size_t r = 0; r < dim; ++r)
			if (v[r][d] != 0 && (ref == 0 || skewed_norm(v[r]) < ref)) ref = skewed_norm(v[r]);
		const double limit = slack > 0 && ref > 0 ? slack * ref : 1e300;

		const size_t use = std::min<size_t>(use_max, dim);
		const size_t base = 2 * (size_t)range + 1;
		size_t combos = 1;
		for (size_t j = 0; j < use; ++j) combos *= base;

		vector<int> e(use); // coefficients of the combination
		for (size_t code = 0; code < combos; ++code)
		{
			// digits of code in given base shifted by range
			size_t x = code;
			int first = 0;
			for (size_t j = 0; j < use; ++j, x /= base)
			{
				e[j] = (int)(x % base) - range;
				if (!first) first = e[j];
			}
			if (first <= 0) continue; // all zero, or negative of combination tried anyway

			bigint lc = 0;
			for (size_t j = 0; j < use; ++j)
				if (e[j]) lc += e[j] * v[j][d];
			if (lc < 0) lc = -lc;
			if (lc == 0 || !good_leading_coefficient(lc)) continue;

			vector<bigint> c(dim, bigint(0));
			for (size_t j = 0; j < use; ++j)
				if (e[j])
					for (size_t i = 0; i < dim; ++i) c[i] += e[j] * v[j][i];
			if (c[d] < 0)
				for (bigint& t : c) t = -t;
			remove_content(c); // f must be primitive: a content of 3 makes factor_base treat 3 as a divisor of n (lc only shrinks, stays good)
			if (boost::multiprecision::gcd(c[d], c[d - 1]) != 1) continue; // p | lc and p | a_(d-1): repeated root at infinity mod p, the factor base assumes a simple one
			if (skewed_norm(c) > limit) continue; // norm of the primitive polynomial, this is the one that gets sieved
			out.push_back(poly(move(c)));
		}
	}

	// score of (f, m), weighted sum of expected relations in sieve lines
	// B, A = sieve region, max_line = largest line used: lines 2, 8, 32, ... (x4) up to max_line, line 2*4^i stands for 3*4^i lines
	// max_line should be of the order of the number of lines the sieve will need, 512 was enough for small n but covers < 1% of the lines at 100 digits
	static double proxy_score(const poly& f, const bigint& m, double B, double A, uint64_t max_line)
	{
		const planner::model M = planner::make_model(f, m, 400);
		double sc = 0;
		for (uint64_t line = 2, width = 3; line <= std::max<uint64_t>(512, max_line); line *= 4, width *= 4)
			sc += (double)width * planner::line_yield(M, B, A, line, 1500);
		return sc;
	}
}

gnfs_poly select_polynomial(const bigint& n, unsigned forced_d, bool verbose, bool always)
{
	using namespace polysel;

	const auto t_start = clk::now();
	const auto elapsed = [&]() { return chrono::duration<double>(clk::now() - t_start).count(); };

	const unsigned d0 = forced_d ? forced_d : choose_degree(n);
	vector<unsigned> degrees;
	if (forced_d) degrees = { forced_d };
	else
		for (unsigned d = d0 > 3 ? d0 - 1 : 3; d <= d0 + 1 && d <= 6; ++d) degrees.push_back(d);

	double Bref = 3.0 * (double)choose_bound(n), Aref = 4.0 * Bref; // proxy region, reset when B0 is re-anchored
	vector<candidate> pool, finalists; // finalists: everything the planner has evaluated

	// the plain base-m polynomials are always candidates
	vector<candidate> basem;
	if (!getenv("GNFS_NONMONIC"))
		for (unsigned d : degrees)
		{
			try
			{
				const gnfs_poly P = make_gnfs_poly(n, d);
				candidate c;
				c.f = P.f; c.m = P.m; c.d = d; c.proven = true;
				basem.push_back(c);
			}
			catch (const exception&) {} // no usable base-m polynomial for this degree
		}

	// the default polynomial's plan fixes the bound the candidates are compared at
	const double env_budget = env_double("GNFS_POLY_SECONDS", -1);
	const bool explicit_budget = env_budget >= 0;

	double B0 = 0, default_seconds = 1e300;
	uint64_t default_lines = 0;
	for (const candidate& c : basem)
		if (c.d == d0)
		{
			const planner::plan_t P = planner::plan_parameters(c.f, c.m, 200);
			if (P.feasible) { B0 = (double)P.bound; default_seconds = P.seconds; default_lines = P.lines; }
		}
	if (!always && !explicit_budget && default_seconds < 25)
	{
		if (verbose) cerr << "polynomial selection skipped: default polynomial predicted " << fmt_seconds(default_seconds) << " s\n";
		for (const candidate& c : basem)
			if (c.d == d0) return { c.f, c.m };
	}
	bool b0_fallback = (B0 == 0); // default polynomial not plannable: B0 comes from choose_bound, which was only calibrated up to ~32 digits
	if (b0_fallback) B0 = 3.0 * (double)choose_bound(n);
	if (verbose)
	{
		cerr << "polynomial selection: candidates are compared at bound B0 = " << (uint64_t)B0 << (b0_fallback ? " (FALLBACK 3 * choose_bound, default polynomial infeasible; re-anchored on a lattice candidate after the grid search)" : " (from default polynomial)");
		if (!b0_fallback) cerr << ", default predicted " << fmt_seconds(default_seconds) << " s, " << default_lines << " lines";
		cerr << '\n';
	}

	// the proxy looks at lines up to ~1/8 of what the default polynomial needs (a good polynomial needs about 10 times fewer lines), at least 512
	uint64_t proxy_max_line = std::min<uint64_t>(1u << 20, std::max<uint64_t>(512, default_lines / 8)); // 512 when the default plan is infeasible (default_lines = 0)

	double best_s = 1e300, planner_time = 0; // best predicted total time among everything evaluated so far
	size_t planner_calls = 0, non_monic_eval = 0;

	// give one candidate to planner, return false if not usable (reducible | duplicate)
	auto evaluate = [&](candidate& c)
	{
		c.tried = true;
		for (const candidate& x : finalists)
			if (x.f.data() == c.f.data()) return false;
		if (!c.proven && is_irreducible(c.f) != irred::proven) return false;
		const auto t0 = clk::now();
		const planner::plan_t P = planner::plan_near(c.f, c.m, B0);
		planner_time += chrono::duration<double>(clk::now() - t0).count();
		++planner_calls;
		c.seconds = P.feasible ? P.seconds : 1e300;
		c.proven = true;
		if (c.seconds < best_s)
		{
			best_s = c.seconds;
			if (verbose)
				cerr << "  " << (long long)elapsed() << " s: new best, degree " << c.d << " lc " << c.f.data().back() << " predicted " << fmt_seconds(c.seconds) << " s\n";
		}
		finalists.push_back(c);
		return true;
	};

	for (candidate& c : basem) evaluate(c);

	// time budget: 1% of the predicted time of the best polynomial found so far, between 20 s and 10 m
	// follows best_s, shrinks when faster polynomial appears
	const double budget_fraction = 0.01, budget_min = 20, budget_max = 1200;
	const auto current_budget = [&]()
	{
		if (explicit_budget) return env_budget;
		if (best_s < 1e299) return std::min(budget_max, std::max(budget_min, budget_fraction * best_s));
		return 1200.0; // nothing evaluated yet
	};

	const double checkpoint_min = 15;
	const size_t min_eval = 8; // minimum number of non-monic candidates evaluated

	double budget = current_budget();
	if (verbose)
	{
		cerr << "polynomial selection: initial time budget " << (long long)budget << " s";
		if (best_s < 1e299) cerr << " (best base-m polynomial predicted at " << fmt_seconds(best_s) << " s)";
		cerr << '\n';
	}

	// degrees that can be searched w/ m0 = n^(1/d)
	vector<unsigned> sdeg;
	vector<bigint> sm0;
	for (unsigned d : degrees)
	{
		const bigint m0 = iroot(n, d);
		if (m0 < 1000) continue;
		sdeg.push_back(d);
		sm0.push_back(m0);
	}

	const size_t pool_keep = 2048;
	set<vector<bigint>> seen; // coefficient vectors already generated
	mt19937_64 rng(0x9E3779B97F4A7C15ull);
	uniform_real_distribution<double> uni(0.0, 1.0);

	// drop evaluated candidates, keep the pool_keep best by proxy
	auto trim_pool = [&]()
	{
		pool.erase(remove_if(pool.begin(), pool.end(), [](const candidate& c) { return c.tried; }), pool.end()); // already in finalists
		sort(pool.begin(), pool.end(), [](const candidate& a, const candidate& b) { return a.proxy > b.proxy; });
		if (pool.size() > pool_keep) pool.resize(pool_keep);
	};

	// evaluate the best count candidates by proxy that have not been given to the planner yet, for at most the given seconds
	auto evaluate_best = [&](size_t count, double seconds = 1e300)
	{
		const double stop = elapsed() + seconds;
		trim_pool();
		for (candidate& c : pool)
		{
			if (count == 0 || elapsed() >= stop) break;
			if (c.tried) continue;
			if (evaluate(c))
			{
				++non_monic_eval;
				--count;
			}
		}
	};

	size_t lattices = 0, tested = 0;

	// one lattice reduction: usable candidates are put in out
	auto run_lattice = [&](unsigned d, const bigint& m, unsigned skew, size_t use_max, int range, double slack, vector<candidate>& out)
	{
		++lattices;
		vector<poly> found;
		lattice_candidates(n, m, d, skew, use_max, range, slack, found);
		for (const poly& f : found)
		{
			if (!good_leading_coefficient(f.data().back()) || f[0] == 0) continue;
			if (boost::multiprecision::gcd(f.data().back(), n) != 1) continue;
			if (modp(f(m), n) != 0) continue;
			if (!seen.insert(f.data()).second) continue;
			candidate c;
			c.f = f; c.m = m; c.d = d; c.skew = skew;
			c.proxy = proxy_score(f, m, Bref, Aref, proxy_max_line);
			out.push_back(move(c));
			++tested;
		}
	};

	// global step: random lattice aimed at useful region
	// for skew s, reduced vectors have lc ~n^(1/(d+1)) / s^(d/2)
	// so a target lc L in [10^3, 10^7] fixes the skew
	// f(m) = 0 mod n with f(m) != 0 needs |f(m)| >= n => m > (n / L)^(1/d),
	// smaller m gives smaller rational side, m is drawn from 1 to 6 * bound
	const double ln_n = msb(n) * 0.6931471806 + 0.35;
	auto global_step = [&]()
	{
		const size_t di = lattices % sdeg.size(); // all degrees get same share of time
		const unsigned d = sdeg[di];
		const double lnL = log(1e3) + uni(rng) * (log(1e7) - log(1e3));
		const double ln_skew = (2.0 / d) * (ln_n / (d + 1) - lnL) + (uni(rng) * 2 - 1) * 0.6931471806;
		const unsigned skew = (unsigned)std::min(4e9, std::max(1.0, exp(ln_skew)));
		const bigint mmin = iroot(n / bigint((unsigned long)exp(lnL)), d);
		const unsigned long scale = (unsigned long)(exp2(2.58 * uni(rng)) * 1048576.0); // 1..6 in fixed point, 2^20 = 1
		const bigint m = mmin * scale / 1048576 + bigint(rng() >> 24); // random offset avoids exact powers
		if (m < 1000) return;
		run_lattice(d, m, skew, 6, 1, 8.0, pool);
		if (pool.size() > 4 * pool_keep) trim_pool();
	};

	// local step around one 'good' polynomial, the best 3 neighbours go to planner
	// translations f(x + t), rotations f + w (x - m), same m at neighbouring skews
	auto local_round = [&]()
	{
		vector<size_t> idx; // evaluated and feasible
		for (size_t i = 0; i < finalists.size(); ++i)
			if (finalists[i].seconds < 1e299) idx.push_back(i);
		if (idx.empty()) return;
		sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return finalists[a].seconds < finalists[b].seconds; });
		const size_t top = std::min<size_t>(5, idx.size());
		const candidate el = finalists[idx[(size_t)(top * uni(rng) * uni(rng))]]; // product of two uniforms: biased towards best

		vector<candidate> nb;

		// translations, |t| up to ~skew
		const double lim = log10(std::max(10.0, (double)el.skew)) + 0.5;
		for (int i = 0; i < 40; ++i)
		{
			long t = lround(pow(10.0, uni(rng) * lim)); // log-uniform
			if (rng() & 1) t = -t;
			if (t == 0) continue;
			vector<bigint> c = el.f.data();
			const size_t deg = c.size() - 1;
			const bigint T(t);
			for (size_t a = 0; a < deg; ++a) // Taylor shift: f(x + t)
				for (size_t j = deg; j-- > a;) c[j] += T * c[j + 1];
			const bigint m2 = el.m - T;
			if (m2 < 1000) continue;
			poly f2(move(c));
			if (!seen.insert(f2.data()).second) continue;
			if (modp(f2(m2), n) != 0) continue;
			candidate cc;
			cc.f = f2; cc.m = m2; cc.d = el.d; cc.skew = el.skew; cc.proven = true;
			cc.proxy = proxy_score(cc.f, cc.m, Bref, Aref, proxy_max_line);
			nb.push_back(move(cc));
		}

		// rotations f + w (x - m), different root behaviour mod small primes
		for (int w = -3; w <= 3; ++w)
		{
			if (w == 0 || el.f.data().size() < 3) continue;
			vector<bigint> c = el.f.data();
			c[0] -= w * el.m;
			c[1] += w;
			if (c[0] == 0) continue;
			remove_content(c); // a rotation can introduce a common factor, lc is unchanged by it so gcd(content, n) = 1 still holds
			poly f2(move(c));
			if (!seen.insert(f2.data()).second) continue;
			candidate cc;
			cc.f = f2; cc.m = el.m; cc.d = el.d; cc.skew = el.skew; // irred is checked by evaluate()
			cc.proxy = proxy_score(cc.f, cc.m, Bref, Aref, proxy_max_line);
			nb.push_back(move(cc));
		}

		// neighbouring skews
		for (double factor : { 0.5, 0.7, 1.4, 2.0 })
		{
			const unsigned sk = std::max(1u, (unsigned)(el.skew * factor));
			if (sk != el.skew) run_lattice(el.d, el.m, sk, 4, 3, 8.0, nb);
		}

		sort(nb.begin(), nb.end(), [](const candidate& a, const candidate& b) { return a.proxy > b.proxy; });
		size_t k = 0;
		for (candidate& c : nb)
		{
			if (k >= 3) break;
			if (evaluate(c))
			{
				++k;
				++non_monic_eval;
			}
		}
	};

	// the default polynomial could not be planned, so B0 is only a guess (3 * choose_bound, calibrated up to ~32 digits)
	// at 120 digits that guess is far too large (all 685 evaluations of one run were infeasible), no candidate is feasible there and plan_near only walks upwards when nothing is feasible
	// plans the best pool candidates by proxy with a free bound search; the first feasible plan sets B0 for all later evaluations
	// input: number of pool candidates to try; output: none (B0 and b0_fallback are updated, B0 unchanged if none is feasible)
	set<vector<bigint>> anchored; // candidates already tried as anchor
	auto try_anchor = [&](size_t count)
	{
		if (!b0_fallback) return;
		trim_pool();
		for (const candidate& c : pool)
		{
			if (count == 0) break;
			if (!anchored.insert(c.f.data()).second) continue;
			--count;
			const planner::plan_t P = planner::plan_parameters(c.f, c.m, 200); // free bound search, same call as for the default polynomial
			if (!P.feasible) continue;
			B0 = (double)P.bound;
			b0_fallback = false;
			if (verbose) cerr << "  " << (long long)elapsed() << " s: bound re-anchored at B0 = " << P.bound << " (degree " << c.d << " lc " << c.f.data().back() << ", predicted " << fmt_seconds(P.seconds) << " s)\n";

			Bref = (double)P.bound; // the proxy was looking at the fallback bound and 512 lines, score the pool again in the region that is actually sieved
			Aref = (double)P.A;
			proxy_max_line = std::min<uint64_t>(1u << 20, std::max<uint64_t>(512, P.lines / 8));
			for (candidate& x : pool) x.proxy = proxy_score(x.f, x.m, Bref, Aref, proxy_max_line);

			finalists.erase(remove_if(finalists.begin(), finalists.end(), [](const candidate& x) { return x.seconds >= 1e299; }), finalists.end()); // were judged at the wrong bound
			for (candidate& x : basem) { x.tried = false; evaluate(x); }
			return;
		}
	};

	// stage 1: fixed grid (m = m0 * k / 2 + 12345, skews from stage1_skews, 3 combined vectors), planner evaluation of best 8
	if (!sdeg.empty())
	{
		for (size_t di = 0; di < sdeg.size(); ++di)
			for (unsigned k : { 1u, 2u, 4u, 8u, 16u })
				for (unsigned skew : stage1_skews(ln_n, sdeg[di]))
					run_lattice(sdeg[di], sm0[di] * k / 2 + 12345, skew, 3, 1, 0, pool);
		try_anchor(5); // only does something if B0 is the unreliable fallback
		evaluate_best(8);
		budget = current_budget();
		if (verbose) {
			cerr << "  " << (long long)elapsed() << " s: grid search done (" << lattices << " lattices), best predicted " << fmt_seconds(best_s) << " s, time budget now " << (long long)budget << " s\n";
		}
	}

	// stage 2: global random search
	// best candidates by proxy (plus a few random ones from better half of pool)
	// go to planner for 30% of the interval, then local search for 20%
	double next_checkpoint = elapsed() + checkpoint_min, next_report = elapsed() + 30;
	while (!sdeg.empty() && elapsed() < 0.75 * budget)
	{
		if (elapsed() >= next_checkpoint)
		{
			const double interval = std::max(checkpoint_min, budget_fraction * budget);
			try_anchor(3); // retry with the new best candidates if the grid gave no feasible plan
			evaluate_best(1u << 30, 0.3 * interval);
			trim_pool();
			for (int t = 0; t < 3 && !pool.empty(); ++t)
			{
				candidate& c = pool[rng() % std::max<size_t>(1, pool.size() / 2)];
				if (!c.tried && evaluate(c)) ++non_monic_eval;
			}
			const double t_end = elapsed() + 0.2 * interval;
			do local_round(); while (elapsed() < t_end && elapsed() < 0.75 * budget);

			const double old_budget = budget;
			budget = current_budget();
			next_checkpoint = elapsed() + interval;
			if (verbose && budget < old_budget)
				cerr << "  " << (long long)elapsed() << " s: best predicted " << fmt_seconds(best_s) << " s, time budget now " << (long long)budget << " s\n";
			continue;
		}

		global_step();

		if (verbose && elapsed() >= next_report)
		{
			next_report += 30;
			cerr << "  " << (long long)elapsed() << " s: " << lattices << " lattices, " << tested << " candidates scored, " << planner_calls << " planner evaluations\n";
		}
	}

	// stage 3: evaluate best remaining candidates by proxy until budget is used up
	{
		trim_pool();
		for (candidate& c : pool)
		{
			if (non_monic_eval >= min_eval && elapsed() >= budget) break;
			if (c.tried) continue;
			if (evaluate(c)) ++non_monic_eval;
			budget = current_budget();
		}
	}

	if (finalists.empty()) throw runtime_error("select_polynomial: no polynomial found");

	vector<size_t> order(finalists.size()); // finalists sorted by predicted time
	for (size_t i = 0; i < order.size(); ++i) order[i] = i;
	sort(order.begin(), order.end(), [&](size_t a, size_t b) { return finalists[a].seconds < finalists[b].seconds; });
	const size_t best = order[0];

	if (verbose)
	{
		cerr << "polynomial selection: " << lattices << " lattices, " << tested << " lattice candidates, " << planner_calls << " planner evaluations ("
		     << fixed << setprecision(2) << (planner_calls ? planner_time / planner_calls : 0.0) << " s each), " << (long long)elapsed() << " s total\n";
		for (size_t r = 0; r < order.size() && r < 10; ++r)
		{
			const candidate& c = finalists[order[r]];
			cerr << "  " << (r == 0 ? "* " : "  ") << "degree " << c.d << " lc " << c.f.data().back() << " alpha " << fixed << setprecision(2) << planner::alpha_algebraic(c.f, 1000)
			     << " predicted " << fmt_seconds(c.seconds) << (c.seconds < 1e299 ? " s" : "")
			     << (c.f.data().back() == 1 ? "  (base-m)" : "") << '\n';
		}
		cerr.unsetf(ios::floatfield);
	}

	if (finalists[best].seconds >= 1e299) return make_gnfs_poly(n, d0); // nothing feasible according to model: keep default
	return { finalists[best].f, finalists[best].m };
}
