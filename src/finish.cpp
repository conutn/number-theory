// replays the end of the pipeline: saved state + dependencies -> square root -> factor

#include "common.hpp"
#include "numtheory.hpp"
#include "polymod.hpp"
#include "factor_base.hpp"
#include "sieve.hpp"
#include "sqrt.hpp"
#include "state.hpp"

// usage: ./finish state.rels deps.txt
// state.rels is written by main.cpp when GNFS_DUMP_PREFIX is set (save_state)
// deps.txt is written by `postprocess -o deps.txt` (without --raw): one dependency per line, indices into `mat`
// prints n = p * q and returns 0 on success, 1 if no dependency gave a factor, 2 on bad input
int main(int argc, char** argv)
{
	if (argc < 3) { cerr << "usage: finish state.rels deps.txt\n"; return 2; }

	gnfs_state S;
	if (!load_state(argv[1], S)) { cerr << "cannot read state " << argv[1] << '\n'; return 2; }
	cerr << S.rels.size() << " full + " << S.partials.size() << " partial relations, bound " << S.bound << '\n';

	const bigint n(S.n), m(S.m);
	vector<bigint> coef;
	{
		istringstream is(S.f_text);
		bigint c;
		while (is >> c) coef.push_back(c);
	}
	if (coef.size() < 2) { cerr << "polynomial in the state file is empty or broken\n"; return 2; }
	const poly f(coef.begin(), coef.end());

	// m must still be a root of f mod n, otherwise the polynomial was read incorrectly
	{
		bigint acc = 0;
		for (size_t i = coef.size(); i-- > 0;) acc = (acc * m + coef[i]) % n; // horner
		if (acc != 0) { cerr << "f(m) != 0 mod n: the polynomial was not restored correctly\n"; return 2; }
	}

	const factor_bases fb = make_factor_bases(f, m, S.bound, S.bound, choose_num_qc(n));
	cerr << "factor bases: " << fb.rational.size() << " / " << fb.algebraic.size() << " / " << fb.quadratic.size() << '\n';

	// same construction as in main.cpp
	vector<nfs_relation> mat = S.rels;
	combine_partials(S.partials, mat);
	if (mat.size() != S.mat_rows)
	{
		cerr << "rebuilt matrix has " << mat.size() << " rows, the original had " << S.mat_rows << ": the dependencies would not fit\n";
		return 2;
	}

	vector<nfs_relation> pool = S.rels; // full relations, then partials
	const size_t nfull = S.rels.size();
	pool.insert(pool.end(), S.partials.begin(), S.partials.end());

	// maps a dependency over rows of mat to one over the flat pool
	// a relation used an even number of times drops out
	auto expand = [&](const vector<size_t>& dep) {
		vector<char> flip(pool.size(), 0);
		for (size_t r : dep)
		{
			if (r < nfull) flip[r] ^= 1;
			else for (uint32_t j : mat[r].members) flip[nfull + j] ^= 1; // combined row = product of its partials
		}
		vector<size_t> out;
		for (size_t i = 0; i < flip.size(); ++i) if (flip[i]) out.push_back(i);
		return out;
	};

	std::ifstream in(argv[2]);
	if (!in) { cerr << "cannot read " << argv[2] << '\n'; return 2; }
	vector<vector<size_t>> deps;
	for (string line; std::getline(in, line);)
	{
		istringstream is(line);
		vector<size_t> d;
		size_t x;
		while (is >> x)
		{
			if (x >= mat.size()) { cerr << "dependency index " << x << " is out of range (" << mat.size() << " rows)\n"; return 2; }
			d.push_back(x);
		}
		d = expand(d);
		if (!d.empty()) deps.push_back(std::move(d)); // empty = everything cancelled
	}
	std::sort(deps.begin(), deps.end());
	deps.erase(std::unique(deps.begin(), deps.end()), deps.end());
	cerr << deps.size() << " distinct dependencies\n";

	for (size_t k = 0; k < deps.size(); ++k)
	{
		cerr << "trying dependency " << k + 1 << " of " << deps.size() << " (" << deps[k].size() << " relations)\n";
		bigint factor;
		if (try_dependency(pool, deps[k], f, fb, m, n, factor))
		{
			cout << n << " = " << factor << " * " << n / factor << '\n';
			return 0;
		}
	}
	cerr << "no dependency gave a factor\n";
	return 1;
}
