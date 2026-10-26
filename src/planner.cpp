// parameter planning: model of the number of relations, choice of bound / sieve width / lines, predicted running time

#include "planner.hpp"
#include "numtheory.hpp"

namespace planner
{
	static const double KAPPA = 0.55;

	static const double NEED_PER_COLUMN = 0.98;

	static const double U_STEP = 0.01;

	static const int U_N = 6000; // u up to 60

	static double expi_pos(double x) // Ei(x), x > 0 (series; x <= ~8 here)
	{
		double sum = 0, term = 1;
		for (int k = 1; k < 200; ++k)
		{
			term *= x / k;
			const double t = term / k;
			sum += t;
			if (t < 1e-17 * sum) break;
		}
		return 0.5772156649015329 + std::log(x) + sum;
	}

	// ln rho(u) on a grid: ODE (u rho'(u) = -rho(u-1)) for u <= 5, saddle point approximation beyond
	static const vector<double>& log_rho_table()
	{
		static vector<double> tab;
		if (!tab.empty()) return tab;

		const double H = 0.001;
		const size_t k1 = 1000, N = 5000 + 1;
		vector<double> rho(N, 1.0);
		for (size_t i = k1 + 1; i < N; ++i)
		{
			const double u = i * H;
			rho[i] = rho[i - 1] - H * 0.5 * (rho[i - k1] / u + rho[i - 1 - k1] / (u - H));
		}

		tab.resize(U_N + 1);
		for (int i = 0; i <= U_N; ++i)
		{
			const double u = i * U_STEP;
			if (u <= 1.0) { tab[i] = 0; continue; }
			if (u <= 5.0) { tab[i] = std::log(std::max(rho[(size_t)std::lround(u / H)], 1e-300)); continue; }

			// e^xi = 1 + u xi, xi > 0
			double lo = 1e-3, hi = 80;
			for (int it = 0; it < 200; ++it)
			{
				const double mid = 0.5 * (lo + hi);
				if (std::exp(mid) - 1 - u * mid < 0) lo = mid; else hi = mid;
			}
			const double xi = 0.5 * (lo + hi);
			const double xip = xi / (std::exp(xi) - u);
			tab[i] = 0.5 * std::log(xip / (2 * 3.14159265358979323846)) + expi_pos(xi) - u * xi - std::log(xi);
		}
		return tab;
	}

	static inline double log_rho(double u)
	{
		const vector<double>& t = log_rho_table();
		if (u <= 1.0) return 0;
		const double x = std::min(u, U_N * U_STEP - 0.02) / U_STEP;
		const size_t i = (size_t)x;
		const double fr = x - (double)i;
		return t[i] * (1 - fr) + t[i + 1] * fr;
	}

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

	static double alpha_rational()
	{
		double a = 0;
		for (unsigned p : small_primes(3000)) a += (1.0 / (p - 1) - (double)p / ((double)p * p - 1)) * std::log((double)p);
		return a;
	}

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
			for (uint64_t x = 0; x < p; ++x)
			{
				uint64_t v = 0;
				for (size_t k = d + 1; k-- > 0;) v = (v * x + c[k]) % p;
				if (v == 0) ++q;
			}
			if (lc % p == 0) ++q; // root at infinity
			a += (1.0 / (p - 1) - (double)q * p / ((double)p * p - 1)) * std::log((double)p);
		}
		return a;
	}

	static double coprime_density(uint64_t b)
	{
		double r = 1;
		for (uint64_t p = 2; p * p <= b; ++p)
			if (b % p == 0) { r *= 1 - 1.0 / (double)p; while (b % p == 0) b /= p; }
		if (b > 1) r *= 1 - 1.0 / (double)b;
		return r;
	}

	double line_yield(const model& M, double B, double A, uint64_t b, int npts)
	{
		const size_t d = M.fc.size() - 1;
		const double lb = std::log(B), bd = (double)b;
		vector<double> bp(d + 1);
		bp[0] = 1;
		for (size_t j = 1; j <= d; ++j) bp[j] = bp[j - 1] * -bd;

		double sum = 0;
		for (int i = 0; i < npts; ++i)
		{
			const double a = -A + 2 * A * i / (npts - 1);
			const double x = -a;
			double apow = 1.0;
			double h = 0.0;
			for (size_t k = 0; k <= d; ++k)
			{
				h += M.fc[k] * apow * bp[d - k];
				apow *= x;
			}
			const double uR = (std::log(std::max(std::fabs(a + bd * M.m), 1.0)) + M.alpha_r) / lb;
			const double uA = (std::log(std::max(std::fabs(h), 1.0)) + M.alpha_a) / lb;
			
			sum += std::exp(log_rho(uR) + log_rho(uA));
		}
		return (2 * A + 1) * 6.0 / (M_PI * M_PI) * sum / npts;
	}

	static double pi_approx(double x)
	{
		const double L = std::log(x);
		return x / (L - 1.0 - 1.0 / L);
	}

	// seconds per sieve position.  measured (A = 4B): 14.7 ns at B=28047, 27 ns at 602534, 36 ns at 1205068, ~41 ns at 2.2M;
	// a wider A costs a bit more per position (cache): +17% at A = 8B
	static double ns_per_position(double B, double A)
	{
		double c = (B <= 602534 ? 14.7e-9 * std::pow(B / 28047.0, 0.2) : 27e-9 * std::pow(B / 602534.0, 0.42));
		if (A > 4.0 * B) c *= std::pow(A / (4.0 * B), 0.22);
		return c;
	}

	static double other_seconds(double ncols) // relation processing, matrix, Lanczos, square root
	{
		return 1.0e-6 * std::pow(ncols, 1.5) + 1.5e-9 * ncols * ncols;
	}

	plan_t evaluate(const model& M, double B, int64_t A, int npts = 6000)
	{
		plan_t P;
		P.bound = (uint64_t)B;
		P.A = A;
		P.ncols = 2 * pi_approx(B) + 45;
		P.needed = NEED_PER_COLUMN * P.ncols + 100;

		// yields on a geometric grid of b, integrated (trapezoid) until `needed` relations are reached
		const int G = 48;
		double prev_b = 1, prev_y = KAPPA * line_yield(M, B, (double)A, 1, npts);
		double cum = 0.5 * prev_y;
		if (cum >= P.needed) { P.feasible = true; P.lines = 1; }
		for (int g = 1; g < G && !P.feasible; ++g)
		{
			const double b = std::round(std::pow((double)B_CAP_LINES, (double)g / (G - 1)));
			if (b <= prev_b) continue;
			const double y = KAPPA * line_yield(M, B, (double)A, (uint64_t)b, npts);
			const double seg = 0.5 * (y + prev_y) * (b - prev_b);
			if (cum + seg >= P.needed)
			{
				// solve for the point inside the segment (linear yield)
				const double need = P.needed - cum;
				const double slope = (y - prev_y) / (b - prev_b);
				double t;
				if (std::fabs(slope) < 1e-18) t = need / std::max(prev_y, 1e-18);
				else t = (-prev_y + std::sqrt(std::max(prev_y * prev_y + 2 * slope * need, 0.0))) / slope;
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

	model make_model(const poly& f, const bigint& m, unsigned alpha_pmax)
	{
		model M;
		for (size_t k = 0; k < f.size(); ++k) M.fc.push_back(f[k].convert_to<double>());
		M.m = m.convert_to<double>();
		M.alpha_r = alpha_rational();
		M.alpha_a = alpha_algebraic(f, alpha_pmax);
		return M;
	}

	plan_t plan_parameters(const poly& f, const bigint& m, double start_bound, int64_t max_A)
	{
		const model M = make_model(f, m);
		plan_t best;
		int since_best = 0;
		for (double B = std::max(200.0, start_bound); B < 4e9 && since_best < 5; B *= 1.2)
		{
			bool improved = false;
			for (int k : { 2, 4, 8 })
			{
				const int64_t A = std::max<int64_t>(10000, std::min<int64_t>(max_A, (int64_t)(k * B)));
				const plan_t P = evaluate(M, B, A);
				if (P.feasible && (!best.feasible || P.seconds < best.seconds)) { best = P; improved = true; }
			}
			since_best = improved ? 0 : since_best + (best.feasible ? 1 : 0);
		}
		return best;
	}

	plan_t plan_near(const poly& f, const bigint& m, double B0, int64_t max_A)
	{
		const model M = make_model(f, m);
		plan_t best;
		for (double fct : { 0.2, 0.3, 0.45, 0.7, 1.0, 1.4 })
			for (int k : { 2, 4, 8 })
			{
				const double B = std::max(200.0, B0 * fct);
				const plan_t P = evaluate(M, B, std::max<int64_t>(10000, std::min<int64_t>(max_A, (int64_t)(k * B))), 2500);
				if (P.feasible && (!best.feasible || P.seconds < best.seconds)) best = P;
			}
		return best;
	}

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
