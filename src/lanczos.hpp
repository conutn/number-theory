#pragma once
// linear algebra over F_2: block Lanczos on B*B^T and the dense fallback; both return dependencies among the rows

#include "common.hpp"
#include "matrix.hpp"

// dependencies among rows of B (x with x^TB=0) by block Lanczos on BB^T
// ker(BB^T) is a superset of ker(B^T):
//  1. combine Lanczos outputs [X+Y | V_m] into vectors y with BB^Ty=0
//  2. combine those y so that B^Ty=0
// returns lists of row indices, all linearly independent
// B should be filtered first (see filter_matrix). Tries a few random starts if Lanczos does not converge
vector<vector<size_t>> lanczos_dependencies(const sparse_matrix& B, uint64_t seed = 12345, int max_attempts = 4);

// Gaussian elimination fallback, same output
vector<vector<size_t>> dense_dependencies(const sparse_matrix& B);

// step 1: independent vectors y with BB^Ty=0 from one Lanczos output (xy = X+Y, vm = last V)
vector<bitvec> lanczos_kernel_vectors(const sparse_matrix& B, const vector<uint64_t>& xy, const vector<uint64_t>& vm);

// step 2: dependencies (lists of row indices of B) from a pool of vectors in ker(BB^T), possibly from several runs
vector<vector<size_t>> dependencies_from_pool(const sparse_matrix& B, const vector<bitvec>& pool);

// steps 1 and 2 for one Lanczos output on its own
vector<vector<size_t>> postprocess_lanczos(const sparse_matrix& B, const vector<uint64_t>& xy, const vector<uint64_t>& vm);
