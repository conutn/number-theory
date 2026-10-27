// GNFS polynomial selection: plain base-m polynomials + time-budgeted lattice (LLL) search for better, non-monic ones
//
// the search runs until a time budget is used up. the budget is a fraction of the predicted total time of the BEST
// polynomial found so far (the default one's prediction can be far off) and shrinks whenever a faster one is found
//
// stages:
//  1. fixed grid of lattices (m = m0*k/2, a few skews), its 8 best by proxy score go to the planner
//  2. random global search aimed at a target leading coefficient L (10^3..10^7): the skew follows from L, m is drawn just
//     above the smallest usable value (n/L)^(1/d). at regular checkpoints the best candidates by proxy go to the planner,
//     followed by a local search around the polynomials it likes best: translations f(x+t) (change sizes, not alpha),
//     rotations f + w(x-m) (change alpha) and neighbouring skews at the same m
//  3. final planner evaluation of the best remaining candidates
//
//   GNFS_POLY_SECONDS=<s>   fixed budget in seconds (0 = only compare the base-m polynomials), also disables the early skip
//   GNFS_NONMONIC=1         leave the base-m polynomials out (for testing the non-monic code paths)

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

	return poly(std::move(c));
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
typedef boost::multiprecision::cpp_rational bigrat;
#else
typedef boost::multiprecision::mpq_rational bigrat;
#endif

namespace polysel
{
	typedef std::chrono::steady_clock clk;

	// value of environment variable name, or def if unset or not a number
	static double env_double(const char* name, double def)
	{
		const char* s = getenv(name);
		if (!s || !*s) return def;
		char* end = nullptr;
		const double v = strtod(s, &end);
		return end == s ? def : v;
	}

	// nearest integer to q (ties up)
	static bigint round_rat(const bigrat& q)
	{
		const bigint num = boost::multiprecision::numerator(q), den = boost::multiprecision::denominator(q); // den > 0
		const bigint a = 2 * num + den, b = 2 * den;
		bigint r = a / b;
		if (a % b != 0 && a < 0) --r; // division truncates towards zero, we want floor
		return r;
	}

	// returns LLL reduced basis (delta = 0.99) of the rows of B
	// exact rational Gram-Schmidt, recomputed after every change: fine for dimension <= ~8
	static vector<vector<bigint>> lll(vector<vector<bigint>> B)
	{
		const size_t n = B.size(), dim = B[0].size();
		vector<vector<bigrat>> Bs(n, vector<bigrat>(dim)), mu(n, vector<bigrat>(n, bigrat(0))); // Bs = orthogonalized rows, mu = GS coefficients
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
				Bs[i] = std::move(v);
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

	static const vector<unsigned>& lc_primes()
	{
		static const vector<unsigned> ps = planner::small_primes(1000);
		return ps;
	}

	// true if c is in (0, 10^7], squarefree and all its prime factors are < 1000
	static bool good_leading_coefficient(bigint c)
	{
		if (c <= 0 || c > 10000000) return false;
		for (unsigned p : lc_primes())
			if (c % p == 0)
			{
				c /= p;
				if (c % p == 0) return false; // p^2 | c
			}
		return c == 1;
	}

	struct candidate
	{
		poly f;
		bigint m;
		double proxy = 0;
		double seconds = 0; // predicted total time (planner), 0 = not evaluated
		unsigned d = 0;
		bool tried = false; // already given to the planner
		bool proven = false; // irreducibility already proven (base-m polynomials, translations of evaluated ones)
		unsigned skew = 1; // skew of the lattice it came from
	};

	// appends to out polynomials from the lattice { c : sum c_i m^i = 0 (mod n) } with coordinate i weighted by skew^i:
	// the first use_max reduced vectors, combined with coefficients in [-range, range]
	// combinations are rejected on the leading coefficient first (cheap) and, if slack > 0,
	// when their skewed norm exceeds slack times that of the shortest reduced vector
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
			return std::sqrt(s);
		};
		const double limit = slack > 0 ? slack * skewed_norm(v[0]) : 1e300;

		const size_t use = std::min<size_t>(use_max, dim);
		const size_t base = 2 * (size_t)range + 1;
		size_t combos = 1;
		for (size_t j = 0; j < use; ++j) combos *= base;

		vector<int> e(use); // coefficients of the combination
		for (size_t code = 0; code < combos; ++code)
		{
			// digits of code in the given base, shifted by range
			size_t x = code;
			int first = 0;
			for (size_t j = 0; j < use; ++j, x /= base)
			{
				e[j] = (int)(x % base) - range;
				if (!first) first = e[j];
			}
			if (first <= 0) continue; // all zero, or the negative of a combination that is tried anyway

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
			if (skewed_norm(c) > limit) continue;
			out.push_back(poly(std::move(c)));
		}
	}

	// cheap ranking score of (f, m): weighted sum of expected relations in a few sieve lines
	// B, A = reference sieve region
	static double proxy_score(const poly& f, const bigint& m, double B, double A)
	{
		const planner::model M = planner::make_model(f, m, 400);
		static const uint64_t lines[] = { 2, 8, 32, 128, 512 };
		static const double width[] = { 3, 12, 48, 192, 768 };
		double sc = 0;
		for (int i = 0; i < 5; ++i) sc += width[i] * planner::line_yield(M, B, A, lines[i], 1500);
		return sc;
	}
}

gnfs_poly select_polynomial(const bigint& n, unsigned forced_d, bool verbose, bool always)
{
	using namespace polysel;

	const auto t_start = clk::now();
	const auto elapsed = [&]() { return std::chrono::duration<double>(clk::now() - t_start).count(); };

	const unsigned d0 = forced_d ? forced_d : choose_degree(n);
	vector<unsigned> degrees;
	if (forced_d) degrees = { forced_d };
	else
		for (unsigned d = d0 > 3 ? d0 - 1 : 3; d <= d0 + 1 && d <= 6; ++d) degrees.push_back(d);

	const double Bref = 3.0 * (double)choose_bound(n), Aref = 4.0 * Bref;
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
	for (const candidate& c : basem)
		if (c.d == d0)
		{
			const planner::plan_t P = planner::plan_parameters(c.f, c.m, 200);
			if (P.feasible) { B0 = (double)P.bound; default_seconds = P.seconds; }
		}
	if (!always && !explicit_budget && default_seconds < 25)
	{
		if (verbose) cerr << "polynomial selection skipped: the default polynomial is predicted at " << (long long)default_seconds << " s\n";
		for (const candidate& c : basem)
			if (c.d == d0) return { c.f, c.m };
	}
	if (B0 == 0) B0 = 3.0 * (double)choose_bound(n);

	double best_s = 1e300, planner_time = 0; // best predicted total time among everything evaluated so far
	size_t planner_calls = 0, non_monic_eval = 0;

	// gives one candidate to the planner, returns false if it was not usable (reducible / duplicate)
	auto evaluate = [&](candidate& c)
	{
		c.tried = true;
		for (const candidate& x : finalists)
			if (x.f.data() == c.f.data()) return false;
		if (!c.proven && is_irreducible(c.f) != irred::proven) return false;
		const auto t0 = clk::now();
		const planner::plan_t P = planner::plan_near(c.f, c.m, B0);
		planner_time += std::chrono::duration<double>(clk::now() - t0).count();
		++planner_calls;
		c.seconds = P.feasible ? P.seconds : 1e300;
		c.proven = true;
		if (c.seconds < best_s)
		{
			best_s = c.seconds;
			if (verbose)
				cerr << "  " << (long long)elapsed() << " s: new best, degree " << c.d << " lc " << c.f.data().back() << " predicted " << (long long)c.seconds << " s\n";
		}
		finalists.push_back(c);
		return true;
	};

	for (candidate& c : basem) evaluate(c);

	// time budget: 1% of the predicted time of the best polynomial found so far, between 20 s and 4 h
	// (a polynomial that is only 10% faster already repays several times that much search time)
	// it follows best_s, so it shrinks when a faster polynomial turns up and a badly overestimated default does not blow it up
	const double budget_fraction = 0.01, budget_min = 20, budget_max = 4 * 3600;
	const auto current_budget = [&]()
	{
		if (explicit_budget) return env_budget;
		if (best_s < 1e299) return std::min(budget_max, std::max(budget_min, budget_fraction * best_s));
		return 60.0; // nothing evaluated yet
	};

	const double checkpoint_min = 15;
	const size_t min_eval = 8; // minimum number of non-monic candidates evaluated

	double budget = current_budget();
	if (verbose)
	{
		cerr << "polynomial selection: initial time budget " << (long long)budget << " s";
		if (best_s < 1e299) cerr << " (best base-m polynomial predicted at " << (long long)best_s << " s)";
		cerr << '\n';
	}

	// degrees that can be searched (m large enough) and their m0 = n^(1/d)
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
	std::set<vector<bigint>> seen; // coefficient vectors already generated
	std::mt19937_64 rng(0x9E3779B97F4A7C15ull);
	std::uniform_real_distribution<double> uni(0.0, 1.0);

	// drop evaluated candidates, keep the pool_keep best by proxy
	auto trim_pool = [&]()
	{
		pool.erase(std::remove_if(pool.begin(), pool.end(), [](const candidate& c) { return c.tried; }), pool.end()); // already in finalists
		std::sort(pool.begin(), pool.end(), [](const candidate& a, const candidate& b) { return a.proxy > b.proxy; });
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

	// one lattice reduction: its usable candidates go into out with their proxy score
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
			c.proxy = proxy_score(f, m, Bref, Aref);
			out.push_back(std::move(c));
			++tested;
		}
	};

	// global step: a random lattice aimed at the useful region
	// for skew s the reduced vectors have leading coefficient about n^(1/(d+1)) / s^(d/2),
	// so a target leading coefficient L in [10^3, 10^7] fixes the skew (up to a random factor 2)
	// f(m) = 0 mod n with f(m) != 0 needs |f(m)| >= n, i.e. m >= (n / L)^(1/d),
	// smaller m gives a smaller rational side, so m is drawn from 1 to 6 times that bound
	const double ln_n = msb(n) * 0.6931471806 + 0.35;
	auto global_step = [&]()
	{
		const size_t di = lattices % sdeg.size(); // all degrees get the same share of the time
		const unsigned d = sdeg[di];
		const double lnL = std::log(1e3) + uni(rng) * (std::log(1e7) - std::log(1e3));
		const double ln_skew = (2.0 / d) * (ln_n / (d + 1) - lnL) + (uni(rng) * 2 - 1) * 0.6931471806;
		const unsigned skew = (unsigned)std::min(4e9, std::max(1.0, std::exp(ln_skew)));
		const bigint mmin = iroot(n / bigint((unsigned long)std::exp(lnL)), d);
		const unsigned long scale = (unsigned long)(std::exp2(2.58 * uni(rng)) * 1048576.0); // 1..6 in fixed point, 2^20 = 1
		const bigint m = mmin * scale / 1048576 + bigint(rng() >> 24); // random offset avoids exact powers
		if (m < 1000) return;
		run_lattice(d, m, skew, 4, 3, 8.0, pool);
		if (pool.size() > 4 * pool_keep) trim_pool();
	};

	// local step around one of the best polynomials found so far, the best 3 neighbours by proxy go straight to the planner
	// translations f(x + t) with root m - t (same f(m) mod n, same alpha, different sizes), rotations f + w (x - m),
	// and the same m at neighbouring skews
	auto local_round = [&]()
	{
		vector<size_t> idx; // evaluated and feasible
		for (size_t i = 0; i < finalists.size(); ++i)
			if (finalists[i].seconds < 1e299) idx.push_back(i);
		if (idx.empty()) return;
		std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return finalists[a].seconds < finalists[b].seconds; });
		const size_t top = std::min<size_t>(5, idx.size());
		const candidate el = finalists[idx[(size_t)(top * uni(rng) * uni(rng))]]; // product of two uniforms: biased towards the best

		vector<candidate> nb;

		// translations, |t| up to about the skew
		const double lim = std::log10(std::max(10.0, (double)el.skew)) + 0.5;
		for (int i = 0; i < 40; ++i)
		{
			long t = std::lround(std::pow(10.0, uni(rng) * lim)); // log-uniform
			if (rng() & 1) t = -t;
			if (t == 0) continue;
			vector<bigint> c = el.f.data();
			const size_t deg = c.size() - 1;
			const bigint T(t);
			for (size_t a = 0; a < deg; ++a) // Taylor shift: f(x + t)
				for (size_t j = deg; j-- > a;) c[j] += T * c[j + 1];
			const bigint m2 = el.m - T;
			if (m2 < 1000) continue;
			poly f2(std::move(c));
			if (!seen.insert(f2.data()).second) continue;
			if (modp(f2(m2), n) != 0) continue;
			candidate cc;
			cc.f = f2; cc.m = m2; cc.d = el.d; cc.skew = el.skew; cc.proven = true;
			cc.proxy = proxy_score(cc.f, cc.m, Bref, Aref);
			nb.push_back(std::move(cc));
		}

		// rotations f + w (x - m): same f(m), but different root behaviour mod small primes (alpha changes, unlike with translations)
		for (int w = -3; w <= 3; ++w)
		{
			if (w == 0 || el.f.data().size() < 3) continue;
			vector<bigint> c = el.f.data();
			c[0] -= w * el.m;
			c[1] += w;
			if (c[0] == 0) continue;
			poly f2(std::move(c));
			if (!seen.insert(f2.data()).second) continue;
			candidate cc;
			cc.f = f2; cc.m = el.m; cc.d = el.d; cc.skew = el.skew; // irreducibility is checked by evaluate()
			cc.proxy = proxy_score(cc.f, cc.m, Bref, Aref);
			nb.push_back(std::move(cc));
		}

		// neighbouring skews
		for (double factor : { 0.5, 0.7, 1.4, 2.0 })
		{
			const unsigned sk = std::max(1u, (unsigned)(el.skew * factor));
			if (sk != el.skew) run_lattice(el.d, el.m, sk, 4, 3, 8.0, nb);
		}

		std::sort(nb.begin(), nb.end(), [](const candidate& a, const candidate& b) { return a.proxy > b.proxy; });
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

	// stage 1: fixed grid (m = m0 * k / 2 + 12345, 5 skews, 3 combined vectors), planner evaluation of its best 8
	// this reproduces the plain quick search, so the result can never be worse than that, everything after it only adds candidates
	if (!sdeg.empty())
	{
		for (size_t di = 0; di < sdeg.size(); ++di)
			for (unsigned k : { 1u, 2u, 4u, 8u, 16u })
				for (unsigned skew : { 1u, 10u, 100u, 1000u, 10000u })
					run_lattice(sdeg[di], sm0[di] * k / 2 + 12345, skew, 3, 1, 0, pool);
		evaluate_best(8);
		budget = current_budget();
		if (verbose)
			cerr << "  " << (long long)elapsed() << " s: grid search done (" << lattices << " lattices), best predicted " << (long long)best_s << " s, time budget now " << (long long)budget << " s\n";
	}

	// stage 2: global random search
	// at every checkpoint the best candidates by proxy (plus a few random ones from the better half of the pool, since the
	// proxy is only a rough guide) go to the planner for 30% of the interval, then local search rounds for 20%
	double next_checkpoint = elapsed() + checkpoint_min, next_report = elapsed() + 30;
	while (!sdeg.empty() && elapsed() < 0.75 * budget)
	{
		if (elapsed() >= next_checkpoint)
		{
			const double interval = std::max(checkpoint_min, budget_fraction * budget);
			evaluate_best(1u << 30, 0.3 * interval); // the planner is cheap: give it a share of the time instead of a fixed count
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
				cerr << "  " << (long long)elapsed() << " s: best predicted " << (long long)best_s << " s, time budget now " << (long long)budget << " s\n";
			continue;
		}

		global_step();

		if (verbose && elapsed() >= next_report)
		{
			next_report += 30;
			cerr << "  " << (long long)elapsed() << " s: " << lattices << " lattices, " << tested << " candidates scored, " << planner_calls << " planner evaluations\n";
		}
	}

	// stage 3: evaluate the best remaining candidates by proxy until the (current) budget is used up
	// at least min_eval non-monic candidates are always evaluated
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
	std::sort(order.begin(), order.end(), [&](size_t a, size_t b) { return finalists[a].seconds < finalists[b].seconds; });
	const size_t best = order[0];

	if (verbose)
	{
		cerr << "polynomial selection: " << lattices << " lattices, " << tested << " lattice candidates, " << planner_calls << " planner evaluations ("
		     << std::fixed << std::setprecision(2) << (planner_calls ? planner_time / planner_calls : 0.0) << " s each), " << (long long)elapsed() << " s total\n";
		for (size_t r = 0; r < order.size() && r < 10; ++r)
		{
			const candidate& c = finalists[order[r]];
			cerr << "  " << (r == 0 ? "* " : "  ") << "degree " << c.d << " lc " << c.f.data().back() << " alpha " << std::fixed << std::setprecision(2) << planner::alpha_algebraic(c.f, 1000)
			     << " predicted " << (c.seconds < 1e299 ? std::to_string((long long)c.seconds) + " s" : std::string("infeasible"))
			     << (c.f.data().back() == 1 ? "  (base-m)" : "") << '\n';
		}
		cerr.unsetf(std::ios::floatfield);
	}

	if (finalists[best].seconds >= 1e299) return make_gnfs_poly(n, d0); // nothing feasible according to the model: keep the default
	return { finalists[best].f, finalists[best].m };
}
