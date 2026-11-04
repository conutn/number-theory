// parameter planning: model of the number of relations (full relations plus cycles of the large-prime graph),
// choice of bound / sieve width / lines, predicted running time

#include "planner.hpp"
#include "numtheory.hpp"

namespace planner
{
	// calibration constants (fitted to the counts printed by nine runs of the current multi-threaded sieve, bounds 33k to 13.5M, see the notes at the end of the file)
	static const double KAPPA_NOLP = 0.55; // fraction of the ideal full relations found when large primes are off (old measurement, not re-checked)
	static const double KAPPA_FULL = 0.85; // same with large primes on (0.67 to 1.03 per run, no trend with the bound; the old sieve found 0.72 to 0.75)
	static const double KAPPA_PART = 0.76; // fraction of the ideal partial relations (one large prime per side) found (0.58 to 0.83 per run, the old sieve found 0.65 to 0.69)
	static const double U_STEP = 0.01; // grid step of the rho table
	static const int U_N = 6000; // u up to 60
	static const int NB = 10; // bins of the large prime size (log scale) used by the graph model

	// Ei(x) for x > 0, by the power series (x <= ~8 here)
	static double expi_pos(double x)
	{
		double sum = 0, term = 1;
		for (int k = 1; k < 200; ++k)
		{
			term *= x / k; // x^k / k!
			const double t = term / k;
			sum += t;
			if (t < 1e-17 * sum) break;
		}
		return 0.5772156649015329 + log(x) + sum; // euler gamma + ln x + series
	}

	// ln rho(u) on a grid of step U_STEP, u in [0, U_N*U_STEP]
	// ODE (u rho'(u) = -rho(u-1)) for u <= 5, saddle point approximation beyond
	static const vector<double>& log_rho_table()
	{
		static vector<double> tab;
		if (!tab.empty()) return tab;

		const double H = 0.001; // ODE step
		const size_t k1 = 1000, N = 5000 + 1; // k1 = 1/H (index offset for u-1), N covers u in [0,5]
		vector<double> rho(N, 1.0); // rho = 1 on [0,1]
		for (size_t i = k1 + 1; i < N; ++i)
		{
			const double u = i * H;
			rho[i] = rho[i - 1] - H * 0.5 * (rho[i - k1] / u + rho[i - 1 - k1] / (u - H)); // trapezoid rule
		}

		tab.resize(U_N + 1);
		for (int i = 0; i <= U_N; ++i)
		{
			const double u = i * U_STEP;
			if (u <= 1.0) { tab[i] = 0; continue; }
			if (u <= 5.0) { tab[i] = log(std::max(rho[(size_t)lround(u / H)], 1e-300)); continue; }

			// e^xi = 1 + u xi, xi > 0
			double lo = 1e-3, hi = 80;
			for (int it = 0; it < 200; ++it) // bisection, e^xi - 1 - u xi < 0 left of the root
			{
				const double mid = 0.5 * (lo + hi);
				if (exp(mid) - 1 - u * mid < 0) lo = mid; else hi = mid;
			}
			const double xi = 0.5 * (lo + hi);
			const double xip = xi / (exp(xi) - u); // xi'(u)
			tab[i] = 0.5 * log(xip / (2 * 3.14159265358979323846)) + expi_pos(xi) - u * xi - log(xi);
		}
		return tab;
	}

	// ln rho(u), linear interpolation in the table
	static inline double log_rho(double u)
	{
		const vector<double>& t = log_rho_table();
		if (u <= 1.0) return 0;
		const double x = std::min(u, U_N * U_STEP - 0.02) / U_STEP; // clamp so i+1 stays in the table
		const size_t i = (size_t)x;
		const double fr = x - (double)i;
		return t[i] * (1 - fr) + t[i + 1] * fr;
	}

	// rho(u), 0 for u < 0 (a cofactor cannot be larger than the number it divides)
	static inline double rho(double u)
	{
		return u < 0 ? 0.0 : exp(log_rho(u));
	}

	// all primes below limit
	vector<unsigned> small_primes(unsigned limit)
	{
		vector<char> comp(limit, 0);
		vector<unsigned> ps;
		for (unsigned i = 2; i < limit; ++i)
		{
			if (comp[i]) continue;
			ps.push_back(i);
			for (unsigned j = i * i; j < limit; j += i) comp[j] = 1;
		}
		return ps;
	}

	// Murphy's alpha of g = x - m (one root mod every p)
	static double alpha_rational()
	{
		double a = 0;
		for (unsigned p : small_primes(3000)) a += (1.0 / (p - 1) - (double)p / ((double)p * p - 1)) * log((double)p);
		return a;
	}

	// Murphy's alpha of f: sum_p (1/(p-1) - q_p p/(p^2-1)) ln p over the primes below pmax
	// q_p = number of roots of f mod p (the root at infinity counts when p | lc(f))
	double alpha_algebraic(const poly& f, unsigned pmax)
	{
		double a = 0;
		const size_t d = f.size() - 1;
		const bigint lc = f.data().back();
		for (unsigned p : small_primes(pmax))
		{
			vector<uint64_t> c(d + 1);
			for (size_t k = 0; k <= d; ++k) c[k] = mod_u64(f[k], p);
			unsigned q = 0;
			for (uint64_t x = 0; x < p; ++x) // count roots mod p by brute force
			{
				uint64_t v = 0;
				for (size_t k = d + 1; k-- > 0;) v = (v * x + c[k]) % p; // horner
				if (v == 0) ++q;
			}
			if (lc % p == 0) ++q; // root at infinity
			a += (1.0 / (p - 1) - (double)q * p / ((double)p * p - 1)) * log((double)p);
		}
		return a;
	}

	// u values (log of the norm / log B) of the rational and algebraic side at npts evenly spaced a in [-A, A] on line b
	// uR, uA are resized to npts
	static void line_u(const model& M, double B, double A, uint64_t b, int npts, vector<double>& uR, vector<double>& uA)
	{
		const size_t d = M.fc.size() - 1;
		const double lb = log(B), bd = (double)b;
		vector<double> bp(d + 1);
		bp[0] = 1;
		for (size_t j = 1; j <= d; ++j) bp[j] = bp[j - 1] * -bd; // (-b)^j

		uR.resize(npts);
		uA.resize(npts);
		for (int i = 0; i < npts; ++i)
		{
			const double a = -A + 2 * A * i / (npts - 1); // evenly spaced sample in [-A, A]
			const double x = -a;
			double apow = 1.0;
			double h = 0.0; // algebraic norm of a+b*alpha, up to sign
			for (size_t k = 0; k <= d; ++k)
			{
				h += M.fc[k] * apow * bp[d - k];
				apow *= x;
			}
			uR[i] = (log(std::max(fabs(a + bd * M.m), 1.0)) + M.alpha_r) / lb;
			uA[i] = (log(std::max(fabs(h), 1.0)) + M.alpha_a) / lb;
		}
	}

	// expected number of fully smooth relations (no large primes) in line b (a in [-A, A]) for the factor base bound B, before any sieve loss
	// npts = number of sampled values of a
	double line_yield(const model& M, double B, double A, uint64_t b, int npts)
	{
		vector<double> uR, uA;
		line_u(M, B, A, b, npts, uR, uA);
		double sum = 0;
		for (int i = 0; i < npts; ++i) sum += exp(log_rho(uR[i]) + log_rho(uA[i])); // both norms smooth
		return (2 * A + 1) * 6.0 / (M_PI * M_PI) * sum / npts; // 6/pi^2 = density of coprime (a,b)
	}

	// expected relations per line, split by type (already multiplied by the sieve efficiencies)
	// r1 / a1: large prime on the rational / algebraic side only, rb / ab: large primes on both sides (rational resp. algebraic prime)
	// the arrays hold the expected count per bin of log(large prime) / log(B) in [1, s]
	struct line_rates
	{
		double full;
		double r1[NB], a1[NB], rb[NB], ab[NB];
	};

	// rates of line b: full relations and, if s > 1, partial relations with one large prime <= B^s per side
	// a number x is B-smooth with probability rho(u); with one extra prime p = B^t (1 < t <= s) it is rho(u - t) / t per unit of t (density of primes ~ 1/(p ln p))
	static void rates_of_line(const model& M, double B, double A, uint64_t b, int npts, double s, double kf, double kp, line_rates& out)
	{
		vector<double> uR, uA;
		line_u(M, B, A, b, npts, uR, uA);

		const bool lp = s > 1.0;
		const double dt = lp ? (s - 1.0) / NB : 0.0;
		double sf = 0, sr1[NB] = {}, sa1[NB] = {}, srb[NB] = {}, sab[NB] = {};
		for (int i = 0; i < npts; ++i)
		{
			const double fr = rho(uR[i]), fa = rho(uA[i]); // no large prime
			sf += fr * fa;
			if (!lp) continue;

			double wr[NB], wa[NB], lr = 0, la = 0; // large prime weight per bin on each side, and the totals
			for (int j = 0; j < NB; ++j)
			{
				const double t = 1.0 + (j + 0.5) * dt; // bin centre
				wr[j] = rho(uR[i] - t) / t * dt;
				wa[j] = rho(uA[i] - t) / t * dt;
				lr += wr[j];
				la += wa[j];
			}
			for (int j = 0; j < NB; ++j)
			{
				sr1[j] += wr[j] * fa;
				sa1[j] += fr * wa[j];
				srb[j] += wr[j] * la; // both sides have a large prime: the edge touches a rational and an algebraic vertex
				sab[j] += wa[j] * lr;
			}
		}

		const double scale = (2 * A + 1) * 6.0 / (M_PI * M_PI) / npts;
		out.full = kf * scale * sf;
		for (int j = 0; j < NB; ++j)
		{
			out.r1[j] = kp * scale * sr1[j];
			out.a1[j] = kp * scale * sa1[j];
			out.rb[j] = kp * scale * srb[j];
			out.ab[j] = kp * scale * sab[j];
		}
	}

	// pi(x) ~ x / (ln x - 1 - 1/ln x)
	static double pi_approx(double x)
	{
		const double L = log(x);
		return x / (L - 1.0 - 1.0 / L);
	}

	// wall clock seconds per sieve position (multi-threaded sieve, large primes on, including candidate factoring and filtering)
	// input: factor base bound B and sieve half-width A; output: seconds for one position of one line
	static double ns_per_position(double B, double A)
	{
		double c = 1.9e-9 * pow(B / 1e6, -0.05); // fitted to the nine runs of the threaded sieve: 2.2 ns at B=118k, 1.9 ns at 1M, 1.64 ns at 13.5M (nearly flat in B, thread overhead is amortised better by the big runs)
		if (A > 4.0 * B) c *= pow(A / (4.0 * B), 0.22); // wider A costs more per position (cache); measured with the old single thread sieve (+17% at A = 8B), not re-measured
		return c;
	}

	// factor of the mean-field cycle count (the model ignores the small components of the graph, they join the big one later and add cycles)
	// input: factor base bound B; output: multiplier for the cycles the model predicts
	static double cycle_gain(double B)
	{
		(void)B; // the nine runs (B = 33k to 13.5M) show no trend with B beyond the scatter of 1.0 to 1.6 per run
		return 1.28; // minimises the error of the cycle counts with KAPPA_FULL / KAPPA_PART above (flat between 1.25 and 1.30)
	}

	// seconds for the work after the sieve: block Lanczos and the square root steps
	// input: number of matrix columns; output: seconds (relation processing and matrix building are not counted here, they are inside the sieve time of ns_per_position)
	static double other_seconds(double ncols)
	{
		const double lanczos = 3.45e-9 * ncols * ncols + 1.25e-5 * ncols; // within 15% of all nine runs (7k to 1.76M columns): 9312 s at 1.76M, 703 s at 453k, 5.0 s at 37k
		const double sqrt_step = 2.35e-7 * pow(ncols, 1.358); // seconds per dependency, within 30% of the runs that took >= 0.2 s: 61 s at 1.76M, 9 s at 453k, 0.3 s at 37k
		return lanczos + 2.0 * sqrt_step; // 1 to 5 dependencies were tried (2.3 on average), 2 is the expectation since each one works with probability 1/2
	}

	// relations (full + cycles) that main.cpp asks for before it builds the matrix
	static double relations_needed(double ncols)
	{
		return TARGET_PER_COLUMN * ncols + TARGET_EXTRA;
	}

	// expected state of the sieve after a number of lines: full relations, partial relations, cycles of the large-prime graph
	// the graph has the vertex 0 ("no large prime"), a vertex per rational prime and one per algebraic ideal
	// a partial relation is an edge, a new edge closes a cycle if both its ends are already connected
	// mean-field model: S_R[j] / S_A[j] = fraction of the vertices in bin j that are connected to vertex 0, grows with the arriving edges
	struct sieve_state
	{
		double full = 0, partial = 0, cycles = 0;
		double SR[NB] = {}, SA[NB] = {};
		double Np[NB]; // number of primes in each bin of the large prime size
	};

	// adds `width` lines of rates r to the state
	static void add_lines(sieve_state& st, const line_rates& r, double width, double gain)
	{
		st.full += r.full * width;

		double eR[NB], eA[NB], bR[NB], bA[NB], er = 0, ea = 0, eb = 0;
		for (int j = 0; j < NB; ++j)
		{
			eR[j] = r.r1[j] * width;
			eA[j] = r.a1[j] * width;
			bR[j] = r.rb[j] * width;
			bA[j] = r.ab[j] * width;
			er += eR[j];
			ea += eA[j];
			eb += bR[j];
		}
		if (er + ea + eb <= 0) return;

		double sRo = 0, sAo = 0, sRb = 0, sAb = 0; // chance that the large prime of a new edge is already connected to vertex 0
		for (int j = 0; j < NB; ++j)
		{
			if (er > 0) sRo += eR[j] / er * st.SR[j];
			if (ea > 0) sAo += eA[j] / ea * st.SA[j];
			if (eb > 0)
			{
				sRb += bR[j] / eb * st.SR[j];
				sAb += bA[j] / eb * st.SA[j];
			}
		}

		const double edges = er + ea + eb;
		const double closing = er * sRo + ea * sAo + eb * sRb * sAb; // one-sided edges join vertex 0; two-sided edges need both ends connected
		st.cycles += std::min(edges, gain * closing);
		st.partial += edges;

		for (int j = 0; j < NB; ++j) // an edge with one connected end connects the other one
		{
			st.SR[j] = std::min(1.0, st.SR[j] + (eR[j] + bR[j] * sAb) * (1 - st.SR[j]) / st.Np[j]);
			st.SA[j] = std::min(1.0, st.SA[j] + (eA[j] + bA[j] * sRb) * (1 - st.SA[j]) / st.Np[j]);
		}
	}

	// plan for the given bound B and half-width A: number of lines needed and predicted time
	// lp = large prime setting (lp.bound(B) is the bound per side, 0 = off)
	// feasible = false if B_CAP_LINES lines are not enough
	plan_t evaluate(const model& M, double B, int64_t A, const lp_opt& lp, int npts)
	{
		plan_t P;
		P.bound = (uint64_t)B;
		P.A = A;
		P.ncols = 2 * pi_approx(B) + 45; // rational + algebraic factor base, plus extra columns
		P.needed = relations_needed(P.ncols);

		const double lb = log(B);
		const double Lp = std::min(lp.bound(B), B * B); // a cofactor <= B^2 is automatically prime
		const double s = Lp > B ? log(Lp) / lb : 1.0; // large primes are <= B^s
		const double kf = s > 1.0 ? KAPPA_FULL : KAPPA_NOLP, kp = s > 1.0 ? KAPPA_PART : 0.0;

		const double gain = cycle_gain(B);
		sieve_state st;
		for (int j = 0; j < NB; ++j) // primes in bin j: pi(B^t_hi) - pi(B^t_lo)
		{
			const double t0 = 1.0 + (s - 1.0) * j / NB, t1 = 1.0 + (s - 1.0) * (j + 1) / NB;
			st.Np[j] = std::max(1.0, pi_approx(exp(t1 * lb)) - pi_approx(exp(t0 * lb)));
		}

		// rates on a geometric grid of b, integrated (piecewise linear rates) until `needed` relations are reached
		const int G = 40, NSUB = 3;
		vector<uint64_t> grid{ 1 };
		for (int g = 1; g < G; ++g)
		{
			const uint64_t b = (uint64_t)llround(pow((double)B_CAP_LINES, (double)g / (G - 1)));
			if (b > grid.back()) grid.push_back(b); // rounding can repeat small b
		}

		line_rates cur, nxt;
		rates_of_line(M, B, (double)A, 1, npts, s, kf, kp, cur);
		double done = 0; // lines processed so far

		// advances by `width` lines at rates r; returns true (and the exact number of lines) if `needed` is reached inside
		auto advance = [&](const line_rates& r, double width, double& lines) -> bool
		{
			const double before = st.full + st.cycles;
			sieve_state keep = st;
			add_lines(st, r, width, gain);
			const double after = st.full + st.cycles;
			if (after < P.needed) { done += width; return false; }
			const double frac = (P.needed - before) / std::max(after - before, 1e-300); // linear inside the step
			lines = done + frac * width;
			(void)keep;
			done += width;
			return true;
		};

		double lines = 0;
		bool reached = advance(cur, 1.0, lines); // line b = 1
		for (size_t g = 1; g < grid.size() && !reached; ++g)
		{
			const double b0 = (double)grid[g - 1], b1 = (double)grid[g];
			rates_of_line(M, B, (double)A, grid[g], npts, s, kf, kp, nxt);
			for (int k = 0; k < NSUB && !reached; ++k)
			{
				const double w = (k + 0.5) / NSUB;
				line_rates mid;
				mid.full = (1 - w) * cur.full + w * nxt.full;
				for (int j = 0; j < NB; ++j)
				{
					mid.r1[j] = (1 - w) * cur.r1[j] + w * nxt.r1[j];
					mid.a1[j] = (1 - w) * cur.a1[j] + w * nxt.a1[j];
					mid.rb[j] = (1 - w) * cur.rb[j] + w * nxt.rb[j];
					mid.ab[j] = (1 - w) * cur.ab[j] + w * nxt.ab[j];
				}
				reached = advance(mid, (b1 - b0) / NSUB, lines);
			}
			cur = nxt;
		}

		if (reached)
		{
			P.feasible = true;
			P.lines = (uint64_t)ceil(lines);
			P.full = st.full;
			P.partial = st.partial;
			P.cycles = st.cycles;
			P.b_max = std::min<uint64_t>(B_CAP_LINES, 3 * P.lines + 200); // generous: running out of b restarts the whole attempt with a larger bound

			const double positions = (double)P.lines * (2.0 * (double)A + 1);
			P.seconds = positions * ns_per_position(B, (double)A) + other_seconds(P.ncols);
		}
		return P;
	}

	// precompute the model data of (f, m): coefficients as doubles and Murphy's alphas
	model make_model(const poly& f, const bigint& m, unsigned alpha_pmax)
	{
		model M;
		for (size_t k = 0; k < f.size(); ++k) M.fc.push_back(f[k].convert_to<double>());
		M.m = m.convert_to<double>();
		M.alpha_r = alpha_rational();
		M.alpha_a = alpha_algebraic(f, alpha_pmax);
		return M;
	}

	// search the bound and the sieve width for the smallest predicted total time
	// max_A limits the sieve memory (two float arrays of 2A+1 entries, about 12 bytes per position)
	plan_t plan_parameters(const poly& f, const bigint& m, double start_bound, int64_t max_A, const lp_opt& lp)
	{
		const model M = make_model(f, m);
		plan_t best;
		int since_best = 0; // bounds tried since the last improvement (counted once a feasible plan exists)
		for (double B = std::max(200.0, start_bound); B < 4e9 && since_best < 5; B *= 1.2)
		{
			bool improved = false;
			for (int k : { 2, 4, 8 }) // A = k*B
			{
				const int64_t A = std::max<int64_t>(10000, std::min<int64_t>(max_A, (int64_t)(k * B)));
				const plan_t P = evaluate(M, B, A, lp, 1000);
				if (P.feasible && (!best.feasible || P.seconds < best.seconds)) { best = P; improved = true; }
			}
			since_best = improved ? 0 : since_best + (best.feasible ? 1 : 0);
		}
		return best;
	}

	// cheap comparison of polynomials: best predicted plan over a range of bounds and the sieve widths A = 2B, 4B, 8B
	// input: polynomial (f, m), centre bound B0, sieve width limit max_A, large prime setting lp; output: the fastest plan found (feasible = false if none)
	// coarser sampling than plan_parameters (used by select_polynomial)
	// B0 is only the starting point: the best bound moves when the time model changes, so the search keeps walking away from B0 while the best plan sits on the edge of the grid
	// and then refines around the best bound, so that two polynomials are not ranked by the luck of where a coarse grid point falls
	plan_t plan_near(const poly& f, const bigint& m, double B0, int64_t max_A, const lp_opt& lp)
	{
		const model M = make_model(f, m);
		const double STEP = 1.5; // ratio between neighbouring bounds when the grid is extended (the initial grid has ratios 1.4 to 1.55)
		const int MAX_WALK = 8; // extensions of the grid per side and width
		plan_t best;
		for (int k : { 2, 4, 8 }) // A = k*B
		{
			plan_t bk; // best plan for this k
			double bk_f = 0; // its bound as a multiple of B0
			auto consider = [&](double fct) // evaluates B = fct*B0, keeps the plan if it is feasible and faster
			{
				const double B = std::max(200.0, B0 * fct);
				const plan_t P = evaluate(M, B, std::max<int64_t>(10000, std::min<int64_t>(max_A, (int64_t)(k * B))), lp, 400);
				if (P.feasible && (!bk.feasible || P.seconds < bk.seconds)) { bk = P; bk_f = fct; }
			};

			double lo = 0.2, hi = 1.4; // smallest and largest multiple tried so far
			for (double fct : { 0.2, 0.3, 0.45, 0.7, 1.0, 1.4 }) consider(fct);
			for (int it = 0; it < MAX_WALK; ++it)
			{
				if (!bk.feasible || bk_f == hi) { hi *= STEP; consider(hi); } // nothing feasible yet or best at the top edge: larger bounds
				else if (bk_f == lo) { lo /= STEP; consider(lo); } // best at the bottom edge: smaller bounds (an infeasible point ends the walk, bk_f stays at the old edge)
				else break;
			}
			if (bk.feasible) // refine: half a grid step on each side of the best bound
			{
				const double c = bk_f;
				consider(c / sqrt(STEP));
				consider(c * sqrt(STEP));
			}
			if (bk.feasible && (!best.feasible || bk.seconds < best.seconds)) best = bk;
		}
		return best;
	}

	// best sieve width for a given bound (used when the bound is raised after a failed attempt)
	plan_t plan_for_bound(const poly& f, const bigint& m, double B, int64_t max_A, const lp_opt& lp)
	{
		const model M = make_model(f, m);
		plan_t best;
		for (int k : { 2, 4, 8 })
		{
			const plan_t P = evaluate(M, B, std::max<int64_t>(10000, std::min<int64_t>(max_A, (int64_t)(k * B))), lp, 1000);
			if (P.feasible && (!best.feasible || P.seconds < best.seconds)) best = P;
		}
		return best;
	}
}

// calibration notes (nine runs of the multi-threaded sieve, lc 2 to 16027, degree 3 to 5, A = 1.5B to 4B, large prime bound 30 B)
// - method: the model was replayed on the exact polynomials of the runs and evaluated at the number of lines the sieve actually used (not at the planned lines, see below)
//   with KAPPA = 1 the model counts give per-run fits of kf = 0.67 .. 1.03 (log mean 0.88 weighted by size) and kp = 0.58 .. 0.83 (0.76)
// - the old constants (kf 0.74, kp 0.67, cycle_gain 1.16 .. 1.38 growing with B) were fitted to the old sieve; with them the rms log error of the counts is 0.154 (full, partial) and 0.318 (cycles),
//   with kf 0.85, kp 0.76, cycle_gain 1.28 it is 0.119 and 0.207
// - the cycle gain fitted per run is 0.8 .. 1.6 with no trend in B, so it is a constant now
// - the lines the sieve used are multiples of 8 and exceed the planned lines for small runs (312 used for 91 planned at B = 40k, 328 for 138 at 33k, 432 for 329 at 118k, 1120 for 927 at 3.1M),
//   an overshoot of roughly 100 to 250 lines from the batch granularity of the threaded sieve; it is not modelled, it matters below ~1000 lines only (a second or less)
// - time per position, old single thread sieve (total sieve time of the attempt / (lines * (2A+1))): 30.6 ns at B = 204k growing like B^0.375, which also fits the 256 bit run
// - time per position, current multi-threaded sieve (same definition, lines = lines actually sieved): 1.64 ns at B = 13.5M (A = 1.5B), 1.73 ns at 6.5M (A = 2B), 1.92 ns at 3.1M (A = 4B),
//   1.66 ns at 1.05M (A = 4B), 2.15 ns at 610k, 1.94 ns at 204k, 2.20 ns at 118k; below 40k the sieve takes < 1 s and the timer resolution (0.1 s) dominates, so those points were not used
//   the fit 1.9 ns * (B/1e6)^-0.05 is within 10% of every point from 118k up; the A dependence of the threaded sieve could not be resolved from these runs
// - the earlier numbers (CYCLE_GAIN 4.5, kappa 0.9) came from runs whose line counts were rounded up to multiples of 100 by the old sieve; do not reuse them
// - when new measurements arrive, compare "predicted lines" with the "calibration" line printed by main.cpp; kf, kp are the first constants to adjust
// - block Lanczos (nine runs, first number is the time of the whole Lanczos step): 9311.6 s at 1.76M columns, 3187.9 s at 892k, 703.2 s at 453k, 106.8 s at 165k, 35.7 s at 100k, 5.0 s at 37k, 2.0 s at 22k, 0.3 s at 7k, 0.4 s at 8k
//   free power law gives n^1.96, so a*n^2 + b*n was used (a = 3.45e-9, b = 1.25e-5), ratio model/measured 0.87 to 1.15 for every run
// - square root, seconds per dependency: 61.2 s at 1.76M columns, 36.7 s at 892k, 9.0 s at 453k, 3.5 s at 165k, 1.6 s at 100k, 0.3 s at 37k, 0.2 s at 22k, < 0.05 s below; fit 2.35e-7 * n^1.358 (ratio 0.77 to 1.26)
//   dependencies tried: 3, 2, 3, 2, 5, 1, 1, 2, 2 (mean 2.3), the model uses 2
// - the old term 1e-6 * n^1.5 for relation processing and the matrix was dropped: the stamp "relation finding finished (matrix built)" already contains that work, so it is part of the fitted sieve time
