// command line driver

#include "common.hpp"
#include "colreduce.hpp"
#include "numtheory.hpp"
#include "polymod.hpp"
#include "factor_base.hpp"
#include "poly_select.hpp"
#include "sieve.hpp"
#include "matrix.hpp"
#include "lanczos.hpp"
#include "sqrt_newton.hpp"
#include "crt_tree.hpp"
#include "sqrt_couveignes.hpp"
#include "sqrt.hpp"
#include "planner.hpp"
#include "dump.hpp"
#include "state.hpp"

// prints a progress line on stderr: [wall clock | total elapsed | since previous stamp] msg
static void stamp(const std::string& msg)
{
	using clock = std::chrono::steady_clock;
	static const clock::time_point t_start = clock::now();
	static clock::time_point t_last = t_start;
	const clock::time_point now = clock::now();
	const double total = std::chrono::duration<double>(now - t_start).count();
	const double lap = std::chrono::duration<double>(now - t_last).count();
	t_last = now;

	const std::time_t wall = std::time(nullptr);
	char buf[16];
	std::strftime(buf, sizeof buf, "%H:%M:%S", std::localtime(&wall));

	std::ostringstream os;
	os << '[' << buf << " | " << std::fixed << std::setprecision(1) << total << "s total | +" << lap << "s] " << msg << '\n';
	std::cerr << os.str();
}

// usage: ./sieve [factor base bound] [sieve half-width A] [max b] < n
//        ./sieve plan < n      only print the parameter plan and the predicted running time
// reads n from stdin, prints n = p * q and returns 0 on success, 1 if no factor was found
// environment: GNFS_POLY=base | select    polynomial construction (default: select from 108 bits up)
//              GNFS_DEGREE=d              force the polynomial degree
//              GNFS_NONMONIC=1            allow non-monic polynomials
//              GNFS_LP=L                  large prime bound per side (default 30*bound, 0 = no large primes)
//              GNFS_DUMP_PREFIX=p         dump the state (p.rels) and the matrix (p.mat) for finish / postprocess
// the planner still assumes no large primes, so its line count is pessimistic when they are on
int main(int argc, char** argv)
{
	bigint n;
	cin >> n;
	stamp("start");

	const char* env_poly = getenv("GNFS_POLY");
	const char* env_deg = getenv("GNFS_DEGREE");
	const char* env_lp = getenv("GNFS_LP");
	const unsigned forced_deg = env_deg ? (unsigned)atoi(env_deg) : 0;
	const bool do_select = env_poly ? std::string(env_poly) == "select" : msb(n) + 1 >= 108;
	const gnfs_poly P = do_select ? select_polynomial(n, forced_deg, true, env_poly != nullptr) : make_gnfs_poly(n, forced_deg);

	const bool plan_only = argc > 1 && std::string(argv[1]) == "plan";
	if (plan_only) argc = 1;

	const bool bound_given = argc > 1;
	const bool A_given = argc > 2;
	const bool bmax_given = argc > 3;
	const long EXCESS = 64; // rows beyond the column count needed before linear algebra
	const bool use_planner = !A_given && (plan_only || msb(n) + 1 >= 108);

	stamp("polynomial selection finished");
	cerr << "f = " << P.f.data() << " m = " << P.m << '\n';

	uint64_t bound = bound_given ? strtoull(argv[1], nullptr, 10) : choose_bound(n);
	int64_t plan_A = 0;
	uint64_t plan_b_max = 0;

	if (use_planner)
	{
		const planner::plan_t plan = planner::plan_parameters(P.f, P.m, 200);
		if (plan.feasible)
		{
			bound = plan.bound;
			plan_A = plan.A;
			plan_b_max = plan.b_max;
			cerr << "plan: bound " << plan.bound << ", A = " << plan.A << ", about " << plan.lines << " lines for ~" << (uint64_t)plan.needed << " relations, "
			     << (uint64_t)plan.ncols << " columns; predicted ~" << (uint64_t)plan.seconds << " s on a machine like the reference (sieve memory ~"
			     << (uint64_t)((2.0 * plan.A + 1) * 12 / 1e6) << " MB)\n";
			if (plan.seconds > 3600) cerr << "warning: that is " << plan.seconds / 3600 << " hours; this implementation (base-m polynomial, one thread) is far from the state of the art at this size\n";
		}
		else
			cerr << "plan: no bound reaches enough relations within " << planner::B_CAP_LINES << " lines for this polynomial; trying the formula bound " << bound << " anyway\n";
	}

	if (plan_only) { stamp("planning finished"); return 0; }

	factor_bases fb;
	size_t ncols = 0;
	vector<nfs_relation> rels;      // full relations (no large prime)
	vector<nfs_relation> partials;  // relations with a large prime on at least one side
	vector<nfs_relation> mat;       // matrix rows: rels followed by one combined row per cycle
	filtered_matrix F;
	bool have_matrix = false;

	// matrix rows = full relations + one combined row per independent cycle of the large-prime graph
	auto build_rows = [&]() {
		mat = rels;
		combine_partials(partials, mat);
	};

	for (int attempt = 0; attempt < (bound_given ? 1 : 6) && !have_matrix; ++attempt)
	{
		if (attempt > 0)
		{
			bound = (uint64_t)(bound * 1.3) + 1; // retry with a 30% larger bound
			if (use_planner)
			{
				const planner::plan_t plan = planner::plan_for_bound(P.f, P.m, (double)bound);
				plan_A = plan.feasible ? plan.A : 0;
				plan_b_max = plan.feasible ? plan.b_max : 0;
			}
		}

		fb = make_factor_bases(P.f, P.m, bound, bound, choose_num_qc(n));
		ncols = num_columns(fb);
		stamp("building factor bases (attempt " + std::to_string(attempt + 1) + ")");
		cerr << "bound " << bound << ": base sizes " << fb.rational.size() << " / " << fb.algebraic.size() << " / " << fb.quadratic.size() << '\n';

		sieve_params sp;
		if (A_given) sp.A = strtoll(argv[2], nullptr, 10); else if (plan_A) sp.A = plan_A;
		if (bmax_given) sp.b_max = strtoull(argv[3], nullptr, 10); else if (plan_b_max) sp.b_max = plan_b_max;
		sp.lp_bound = env_lp ? strtoull(env_lp, nullptr, 10) : 30 * bound;
		const uint64_t b_limit = sp.b_max ? sp.b_max : DEFAULT_B_MAX;

		rels.clear();
		partials.clear(); // relations of an earlier attempt belong to another factor base
		uint64_t next_b = 1;
		sp.target = (size_t)(0.93 * ncols) + 64; // full relations + cycles together must reach this

		for (;;)
		{
			stamp("relation finding: sieving batch");
			find_relations(P.f, P.m, fb, sp, rels, partials, next_b);
			build_rows();
			F = filter_matrix(mat, ncols, 2 * EXCESS);
			const long excess = (long)F.B.nrows - (long)F.B.ncols;
			cerr << rels.size() << " full + " << partials.size() << " partial relations -> " << mat.size() << " rows; after filtering "
			     << F.B.nrows << " rows x " << F.B.ncols << " columns (excess " << excess << ")\n";
			if (excess >= EXCESS) { have_matrix = true; stamp("relation finding finished (matrix built)"); break; }
			if (next_b > b_limit) break; // b range exhausted
			sp.target = mat.size() + (size_t)(1.5 * (EXCESS - excess)) + ncols / 50 + 50; // each new full relation or cycle adds ~1 to the excess
		}

		if (!have_matrix) cerr << "not enough relations" << (bound_given ? "\n" : ", retrying with a larger bound\n");
	}
	if (!have_matrix) return 1;

	remove_dependent_columns(F.B);
	if (const char* p = getenv("GNFS_DUMP_PREFIX"))
	{
		ostringstream sn, sm, sf;
		sn << n;
		sm << P.m;
		for (const bigint& c : P.f.data()) sf << c << ' ';
		save_state(string(p) + ".rels", sn.str(), sm.str(), sf.str(), bound, rels, partials, mat.size());
		save_filtered(string(p) + ".mat", F);
	}

	// flat pool of the real (a,b) relations: full relations first, then the partials
	vector<nfs_relation> pool;
	const size_t nfull = rels.size();
	pool.reserve(nfull + partials.size());
	for (nfs_relation& r : rels) pool.push_back(std::move(r));
	for (nfs_relation& r : partials) pool.push_back(std::move(r));
	vector<nfs_relation>().swap(rels);
	vector<nfs_relation>().swap(partials);

	for (size_t i = 0; i < mat.size(); ++i)
	{
		nfs_relation& m = mat[i];
		if (i < nfull)
		{
			m = nfs_relation();
			continue;
		}
		vector<uint32_t>().swap(m.cols);
		vector<fac_entry>().swap(m.rfac);
		vector<fac_entry>().swap(m.afac);
	}
	// maps a dependency over rows of the filtered matrix to one over the flat pool
	// a combined row stands for the product of its members, and a relation used an even number of times is a square: it drops out
	auto expand = [&](const vector<size_t>& dep) {
		vector<char> flip(pool.size(), 0);
		for (size_t i : dep)
		{
			const size_t r = F.row_map[i]; // index into mat
			if (r < nfull) flip[r] ^= 1;
			else for (uint32_t j : mat[r].members) flip[nfull + j] ^= 1;
		}
		vector<size_t> out;
		for (size_t i = 0; i < flip.size(); ++i) if (flip[i]) out.push_back(i);
		return out;
	};

	// tries every distinct dependency (square root step), name is only used in the progress output
	// returns true and prints n = p * q as soon as one gives a factor
	auto try_all = [&](const vector<vector<size_t>>& deps_in, const char* name) {
		vector<vector<size_t>> deps;
		for (const auto& dep : deps_in) { vector<size_t> e = expand(dep); if (!e.empty()) deps.push_back(std::move(e)); }
		std::sort(deps.begin(), deps.end());
		deps.erase(std::unique(deps.begin(), deps.end()), deps.end());
		cerr << name << ": " << deps.size() << " distinct dependencies\n";
		for (const auto& dep : deps)
		{
			stamp(std::string(name) + ": trying a dependency (square root step)");
			bigint factor;
			if (try_dependency(pool, dep, P.f, fb, P.m, n, factor))
			{
				cout << n << " = " << factor << " * " << n / factor << '\n';
				return true;
			}
		}
		return false;
	};

	stamp("linear algebra: block Lanczos starting");
	const auto lanczos_deps = lanczos_dependencies(F.B);
	stamp("block Lanczos finished (" + std::to_string(lanczos_deps.size()) + " raw dependencies)");
	if (try_all(lanczos_deps, "block Lanczos")) { stamp("factoring finished"); return 0; }

	if (F.B.ncols <= 40000) // dense elimination is only feasible for small matrices
	{
		cerr << "falling back to dense elimination\n";
		const auto dense_deps = dense_dependencies(F.B);
		stamp("dense elimination finished");
		if (try_all(dense_deps, "dense elimination")) { stamp("factoring finished"); return 0; }
	}

	stamp("factoring failed");
	cerr << "no dependency gave a factor; collect more relations\n";
	return 1;
}
