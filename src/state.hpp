#pragma once
// binary dump / load of everything the square-root step needs besides the matrix (which dump.hpp handles):
// the polynomial and n as text, the factor base bound, the full and the partial relations, and the number of rows of `mat`.
// With a state file and the matrix file, the end of the pipeline (dependencies -> square root) can be rerun
// without sieving.

#include "common.hpp"
#include "sieve.hpp"
#include <string>

struct gnfs_state
{
	std::string n, m, f_text; // decimal n, decimal m, f as printed by `cerr << P.f.data()`
	uint64_t bound = 0;
	vector<nfs_relation> rels, partials;
	uint64_t mat_rows = 0;    // rows of `mat` (rels + one row per cycle) when the dump was made, to check the rebuild
};

void save_state(const std::string& path, const std::string& n, const std::string& m, const std::string& f_text, uint64_t bound,
                const vector<nfs_relation>& rels, const vector<nfs_relation>& partials, uint64_t mat_rows);
bool load_state(const std::string& path, gnfs_state& S);
