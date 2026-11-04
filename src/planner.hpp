#pragma once
// parameter planning: model of the number of relations, choice of bound / sieve width / lines, predicted running time
//
// model, relations in line b (a in [-A, A], gcd(a,b) = 1):
//   u_R = (ln|a+b*m| + alpha_R) / ln B,   u_A = (ln|N(a,b)| + alpha_A) / ln B,   rho = Dickman function
//   alpha = Murphy's alpha: how many fewer small prime factors the norms have than random integers (computed from
//   the root counts of f mod small primes)
//   a side is B-smooth with probability rho(u); with one extra prime p = B^t (1 < t <= s, B^s = large prime bound) the
//   density is rho(u - t) / t per unit of t, so every line yields four kinds of relations:
//     full (no large prime)               kappa_full * sum rho(u_R) rho(u_A)
//     one large prime on one side         kappa_part * sum rho(u - t)/t ...
//     large primes on both sides          (partial relations are edges of the large-prime graph)
//   kappa_full ~ 0.74, kappa_part ~ 0.67 with large primes on (0.55 without): the sieve loses candidates
//   (algebraic prime powers are not sieved, slack thresholds ...)
// cycles: partial relations only count once they close a cycle in the graph (vertex 0 = no large prime, one vertex per
//   large prime / ideal). The number of cycles is estimated by a mean-field model of the part of the graph that is connected
//   to vertex 0, times a small correction (1.1 - 1.6, see cycle_gain in planner.cpp)
// needed: TARGET_PER_COLUMN * columns + TARGET_EXTRA of (full relations + cycles), columns ~ 2 pi(B)
//   (the same target main.cpp gives the sieve)
// cost: sieve time (lines * positions * ns per position, large prime processing included) + a term growing with the matrix size
//   + the square root steps
// the best (B, A) is the one with the smallest predicted total time.
// calibration: a 192 bit run (degree 4, lc 570, bounds 204k - 758k, A = 2B..4B, large prime bound 30 B): full and partial counts
//   within 5% at nine checkpoints, predicted 349 lines against 346 used at B = 757941. The constants are in planner.cpp,
//   where the measurements they come from are listed; main.cpp prints a "calibration:" line (predicted lines against lines
//   used) at the end of the relation search so they can be refitted. Little data exists for other sizes.

#include "common.hpp"

namespace planner
{
	// relations (full + cycles) wanted before the matrix is built, shared by main.cpp (sieve target) and the planner
	static const double TARGET_PER_COLUMN = 0.93;
	static const double TARGET_EXTRA = 64;

	// large prime setting: bound per side = abs if abs != 0, else mult * B; both 0 = no large primes
	struct lp_opt
	{
		double mult = 30.0;
		uint64_t abs = 0;
		double bound(double B) const { return abs ? (double)abs : mult * B; }
	};

	// all primes below limit
	vector<unsigned> small_primes(unsigned limit);

	// Murphy's alpha of f: sum_p (1/(p-1) - q_p p/(p^2-1)) ln p over the primes below pmax
	// q_p = number of roots of f mod p (the root at infinity counts when p | lc(f))
	double alpha_algebraic(const poly& f, unsigned pmax = 3000);

	// what the yield model needs to know about a polynomial pair (f, g = x - m)
	struct model
	{
		vector<double> fc; // coefficients of f as doubles
		double m;
		double alpha_r, alpha_a;
	};

	// expected number of fully smooth relations (no large primes, no sieve loss) in line b (a in [-A, A]) for the factor base bound B
	// npts = number of sampled values of a
	double line_yield(const model& M, double B, double A, uint64_t b, int npts = 6000);

	// a plan: bound, sieve half-width A, and the predicted number of lines, relations and seconds
	struct plan_t
	{
		bool feasible = false;
		uint64_t bound = 0;
		int64_t A = 0;
		uint64_t lines = 0; // predicted number of lines (b range) to collect enough relations
		uint64_t b_max = 0; // limit to give the sieve (generous: unused lines cost nothing, running out restarts the attempt)
		double ncols = 0, needed = 0, seconds = 0;
		double full = 0, partial = 0, cycles = 0; // predicted full relations, partial relations and cycles when `lines` are done
	};

	inline constexpr uint64_t B_CAP_LINES = 2000000;

	// precompute the model data of (f, m): coefficients as doubles and Murphy's alphas
	model make_model(const poly& f, const bigint& m, unsigned alpha_pmax = 3000);

	// plan for the bound B and half-width A with the given large prime setting (lp.bound(B) per side)
	// npts = number of sampled values of a per line
	plan_t evaluate(const model& M, double B, int64_t A, const lp_opt& lp = lp_opt(), int npts = 1000);

	// search the bound and the sieve width for the smallest predicted total time
	// max_A limits the sieve memory (two float arrays of 2A+1 entries, about 12 bytes per position)
	plan_t plan_parameters(const poly& f, const bigint& m, double start_bound, int64_t max_A = 20000000, const lp_opt& lp = lp_opt());

	// cheap comparison of polynomials: best predicted time among a few bounds around B0
	// coarser sampling than plan_parameters (used by select_polynomial)
	plan_t plan_near(const poly& f, const bigint& m, double B0, int64_t max_A = 20000000, const lp_opt& lp = lp_opt());

	// best sieve width for a given bound (used when the bound is raised after a failed attempt)
	plan_t plan_for_bound(const poly& f, const bigint& m, double B, int64_t max_A = 20000000, const lp_opt& lp = lp_opt());
}
