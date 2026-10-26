#pragma once
// parameter planning: model of the number of relations, choice of bound / sieve width / lines, predicted running time
// parameter planning: choose the factor base bound B, the sieve half-width A and the number of lines (b range)
//
// model (kappa fitted on 5 real runs, 32-77 digits, degree 3 and 4; then checked window by window against the real sieve on a
// 56-digit number at bounds 602534, 1205068 and 2200000 and A = 2B..8B: 24 windows, median model/observed 1.00,
// windows with >= 400 relations within 0.87-1.07.  An end-to-end 56-digit run took 1193 s against a predicted 1253 s):
//   relations in line b  ~  kappa * sum_{|a|<=A, gcd(a,b)=1} rho(u_R) rho(u_A)
//   u_R = (ln|a+b*m| + alpha_R) / ln B,   u_A = (ln|N(a,b)| + alpha_A) / ln B,   rho = Dickman function
//   alpha = Murphy's alpha: how many fewer small prime factors the norms have than random integers (computed from
//   the root counts of f mod small primes).  kappa ~ 0.62 (the sieve loses candidates: algebraic prime powers are
//   not sieved, slack thresholds ...)
// needed: ~0.95 * columns relations (measured 0.86-0.98), columns ~ 2 pi(B)
// cost: sieve time (lines * positions * ns per position) + a term growing with the matrix size
// the best (B, A) is the one with the smallest predicted total time.
// accuracy: predicted relation counts were within ~15% on the test runs; predicted times within ~2x (machine dependent)

#include "common.hpp"

namespace planner
{
	// all primes below limit
	vector<unsigned> small_primes(unsigned limit);

	// Murphy's alpha of f (f monic: no projective roots): sum_p (1/(p-1) - q_p p/(p^2-1)) ln p, q_p = number of roots mod p
	// Murphy's alpha of f, summed over the primes below pmax (the root at infinity counts when p | lc(f))
	double alpha_algebraic(const poly& f, unsigned pmax = 3000);

	// what the yield model needs to know about a polynomial pair (f, g = x - m)
	struct model
	{
		vector<double> fc; // coefficients of f as doubles
		double m;
		double alpha_r, alpha_a;
	};

	// expected number of relations in line b for bound B and half-width A
	// expected number of relations in line b (a in [-A, A]) for the factor base bound B; npts = sampled values of a
	double line_yield(const model& M, double B, double A, uint64_t b, int npts = 6000);

	// a plan: bound, sieve half-width A, and the predicted number of lines, relations and seconds
	struct plan_t
	{
		bool feasible = false;
		uint64_t bound = 0;
		int64_t A = 0;
		uint64_t lines = 0;     // predicted number of lines (b range) to collect enough relations
		uint64_t b_max = 0;     // limit to give the sieve (generous: unused lines cost nothing)
		double ncols = 0, needed = 0, seconds = 0;
	};

	inline constexpr uint64_t B_CAP_LINES = 200000;

	// precompute the model data of (f, m): coefficients as doubles and Murphy's alphas
	model make_model(const poly& f, const bigint& m, unsigned alpha_pmax = 3000);

	// search the bound and width: the cheapest predicted total time.  A <= max_A (memory: two float arrays of 2A+1 entries)
	// search the bound and the sieve width for the smallest predicted total time.
	// max_A limits the sieve memory (two float arrays of 2A+1 entries, about 12 bytes per position)
	plan_t plan_parameters(const poly& f, const bigint& m, double start_bound, int64_t max_A = 20000000);

	// cheap comparison of polynomials: best predicted time among a few bounds around B0 (coarser sampling than plan_parameters)
	// cheap comparison of polynomials: best predicted time among a few bounds around B0 (used by select_polynomial)
	plan_t plan_near(const poly& f, const bigint& m, double B0, int64_t max_A = 20000000);

	// plan for a given bound (used when the bound has to be raised after a failed attempt)
	// best sieve width for a given bound (used when the bound is raised after a failed attempt)
	plan_t plan_for_bound(const poly& f, const bigint& m, double B, int64_t max_A = 20000000);
}
