#pragma once
// sparse matrix over F_2: construction, singleton/surplus-row filtering, products with B and B^T

#include "common.hpp"
#include "sieve.hpp"

typedef vector<uint64_t> bitvec; // word i of a block vector holds entry i of 64 F_2 vectors

// rows = relations, columns = factor base elements
struct sparse_matrix
{
	size_t nrows = 0, ncols = 0;
	vector<vector<uint32_t>> rows; // columns with a 1 in each row
};

// one row per relation
sparse_matrix build_matrix(const vector<nfs_relation>& rels, size_t ncols);

struct filtered_matrix
{
	sparse_matrix B;
	vector<size_t> row_map; // row of B -> index into rels
};

// build the matrix without rows that cannot be in any dependency
// singleton removal: a column that occurs in exactly one row forces that row out of every dependency,
// so drop the row, repeat (removing a row can create new singletons). Empty columns are dropped and the rest renumbered
// without this the matrix is very rank deficient (hundreds of unused columns), V^T A V in Lanczos has tiny rank and the iteration crawls
// if keep_excess >= 0, surplus rows are also dropped (heaviest first, then singleton removal again) until rows - columns <= keep_excess
// block Lanczos is designed for only a few dozen more rows than columns; with thousands of surplus rows A=BB^T has a huge kernel and it did not converge
filtered_matrix filter_matrix(const vector<nfs_relation>& rels, size_t ncols, long keep_excess = -1);

// symmetric operator A=BB^T (nrows x nrows)
// requires v.size() == B.nrows
vector<uint64_t> mul_A(const sparse_matrix& B, const vector<uint64_t>& v);

// B^Ty for a single vector y (bit i of y = entry i)
bitvec mul_Bt_bits(const sparse_matrix& B, const bitvec& y);
