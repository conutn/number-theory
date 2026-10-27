// parameter planning: model of the number of relations, choice of bound / sieve width / lines, predicted running time

#include "planner.hpp"
#include "numtheory.hpp"

namespace planner
{
	static const double KAPPA = 0.55; // fraction of the ideal relation count the sieve actually finds
	static const double NEED_PER_COLUMN = 0.98; // relations needed per matrix column
	static const double U_STEP = 0.01; // grid step of the rho table
	static const int U_N = 6000; // u up to 60

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
		return 0.5772156649015329 + std::log(x) + sum; // euler gamma + ln x + series
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
			if (u <= 5.0) { tab[i] = std::log(std::max(rho[(size_t)std::lround(u / H)], 1e-300)); continue; }

			// e^xi = 1 + u xi, xi > 0
			double lo = 1e-3, hi = 80;
			for (int it = 0; it < 200; ++it) // bisection, e^xi - 1 - u xi < 0 left of the root
			{
				const double mid = 0.5 * (lo + hi);
				if (std::exp(mid) - 1 - u * mid < 0) lo = mid; else hi = mid;
			}
			const double xi = 0.5 * (lo + hi);
			const double xip = xi / (std::exp(xi) - u); // xi'(u)
			tab[i] = 0.5 * std::log(xip / (2 * 3.14159265358979323846)) + expi_pos(xi) - u * xi - std::log(xi);
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
		for (unsigned p : small_primes(3000)) a += (1.0 / (p - 1) - (double)p / ((double)p * p - 1)) * std::log((double)p);
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
			a += (1.0 / (p - 1) - (double)q * p / ((double)p * p - 1)) * std::log((double)p);
		}
		return a;
	}

	// expected number of relations in line b (a in [-A, A]) for the factor base bound B
	// npts = number of sampled values of a
	double line_yield(const model& M, double B, double A, uint64_t b, int npts)
	{
		const size_t d = M.fc.size() - 1;
		const double lb = std::log(B), bd = (double)b;
		vector<double> bp(d + 1);
		bp[0] = 1;
		for (size_t j = 1; j <= d; ++j) bp[j] = bp[j - 1] * -bd; // (-b)^j

		double sum = 0;
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
			const double uR = (std::log(std::max(std::fabs(a + bd * M.m), 1.0)) + M.alpha_r) / lb;
			const double uA = (std::log(std::max(std::fabs(h), 1.0)) + M.alpha_a) / lb;

			sum += std::exp(log_rho(uR) + log_rho(uA)); // both norms smooth
		}
		return (2 * A + 1) * 6.0 / (M_PI * M_PI) * sum / npts; // 6/pi^2 = density of coprime (a,b)
	}

	// pi(x) ~ x / (ln x - 1 - 1/ln x)
	static double pi_approx(double x)
	{
		const double L = std::log(x);
		return x / (L - 1.0 - 1.0 / L);
	}

	// seconds per sieve position
	// measured (A = 4B): 14.7 ns at B=28047, 27 ns at 602534, 36 ns at 1205068, ~41 ns at 2.2M
	// a wider A costs a bit more per position (cache): +17% at A = 8B
	static double ns_per_position(double B, double A)
	{
		double c = (B <= 602534 ? 14.7e-9 * std::pow(B / 28047.0, 0.2) : 27e-9 * std::pow(B / 602534.0, 0.42));
		if (A > 4.0 * B) c *= std::pow(A / (4.0 * B), 0.22);
		return c;
	}

	// seconds for relation processing, matrix, Lanczos, square root
	static double other_seconds(double ncols)
	{
		return 1.0e-6 * std::pow(ncols, 1.5) + 1.5e-9 * ncols * ncols;
	}

	// plan for the given bound B and half-width A: number of lines needed and predicted time
	// feasible = false if B_CAP_LINES lines are not enough
	plan_t evaluate(const model& M, double B, int64_t A, int npts = 6000)
	{
		plan_t P;
		P.bound = (uint64_t)B;
		P.A = A;
		P.ncols = 2 * pi_approx(B) + 45; // rational + algebraic factor base, plus extra columns
		P.needed = NEED_PER_COLUMN * P.ncols + 100;

		// yields on a geometric grid of b, integrated (trapezoid) until `needed` relations are reached
		const int G = 48;
		double prev_b = 1, prev_y = KAPPA * line_yield(M, B, (double)A, 1, npts);
		double cum = 0.5 * prev_y;
		if (cum >= P.needed) { P.feasible = true; P.lines = 1; }
		for (int g = 1; g < G && !P.feasible; ++g)
		{
			const double b = std::round(std::pow((double)B_CAP_LINES, (double)g / (G - 1)));
			if (b <= prev_b) continue; // rounding can repeat small b
			const double y = KAPPA * line_yield(M, B, (double)A, (uint64_t)b, npts);
			const double seg = 0.5 * (y + prev_y) * (b - prev_b);
			if (cum + seg >= P.needed)
			{
				// solve for the point inside the segment (linear yield)
				const double need = P.needed - cum;
				const double slope = (y - prev_y) / (b - prev_b);
				double t;
				if (std::fabs(slope) < 1e-18) t = need / std::max(prev_y, 1e-18); // flat segment
				else t = (-prev_y + std::sqrt(std::max(prev_y * prev_y + 2 * slope * need, 0.0))) / slope; // prev_y*t + slope*t^2/2 = need
				P.lines = (uint64_t)std::ceil(prev_b + t);
				P.feasible = true;
			}

			cum += seg;
			prev_b = b;
			prev_y = y;
		}

		if (P.feasible)
		{
			P.b_max = std::min<uint64_t>(B_CAP_LINES, 2 * P.lines + 100);
			P.seconds = (double)P.lines * (2.0 * (double)A + 1) * ns_per_position(B, (double)A) + other_seconds(P.ncols);
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
	plan_t plan_parameters(const poly& f, const bigint& m, double start_bound, int64_t max_A)
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
				const plan_t P = evaluate(M, B, A);
				if (P.feasible && (!best.feasible || P.seconds < best.seconds)) { best = P; improved = true; }
			}
			since_best = improved ? 0 : since_best + (best.feasible ? 1 : 0);
		}
		return best;
	}

	// cheap comparison of polynomials: best predicted time among a few bounds around B0
	// coarser sampling than plan_parameters (used by select_polynomial)
	plan_t plan_near(const poly& f, const bigint& m, double B0, int64_t max_A)
	{
		const model M = make_model(f, m);
		plan_t best;
		for (double fct : { 0.2, 0.3, 0.45, 0.7, 1.0, 1.4 }) // B = fct*B0
			for (int k : { 2, 4, 8 })
			{
				const double B = std::max(200.0, B0 * fct);
				const plan_t P = evaluate(M, B, std::max<int64_t>(10000, std::min<int64_t>(max_A, (int64_t)(k * B))), 2500);
				if (P.feasible && (!best.feasible || P.seconds < best.seconds)) best = P;
			}
		return best;
	}

	// best sieve width for a given bound (used when the bound is raised after a failed attempt)
	plan_t plan_for_bound(const poly& f, const bigint& m, double B, int64_t max_A)
	{
		const model M = make_model(f, m);
		plan_t best;
		for (int k : { 2, 4, 8 })
		{
			const plan_t P = evaluate(M, B, std::max<int64_t>(10000, std::min<int64_t>(max_A, (int64_t)(k * B))));
			if (P.feasible && (!best.feasible || P.seconds < best.seconds)) best = P;
		}
		return best;
	}
}
