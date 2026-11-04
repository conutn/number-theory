#pragma once
// sparse matrix over F_2: construction, singleton/surplus-row filtering, products with B and B^T

#include "common.hpp"
#include "sieve.hpp"

typedef vector<uint64_t> bitvec; // word i of a block vector holds entry i of 64 F_2 vectors

// rows = relations, columns = factor base elements
struct sparse_matrix
{
	size_t nrows = 0, ncols = 0;
	vector<vector<uint32_t>> rows; // columns w/ 1 in each row
};

// one row per relation
sparse_matrix build_matrix(const vector<nfs_relation>& rels, size_t ncols);

struct filtered_matrix
{
	sparse_matrix B;
	vector<size_t> row_map; // row of B -> index into rels
};

// build matrix w/o rows that cannot be in any dependency
// singleton removal: column occuring in exactly one row removes row
// surplus rows are dropped (heaviest first then singleton removal) until rows - columns <= keep_excess
filtered_matrix filter_matrix(const vector<nfs_relation>& rels, size_t ncols, long keep_excess = -1);

// symmetric operator A=BB^T (nrows x nrows)
// requires v.size() == B.nrows
vector<uint64_t> mul_A(const sparse_matrix& B, const vector<uint64_t>& v);

// B^Ty for a single vector y (bit i of y = entry i)
bitvec mul_Bt_bits(const sparse_matrix& B, const bitvec& y);
