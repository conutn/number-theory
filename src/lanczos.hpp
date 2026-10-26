#pragma once
// linear algebra over F_2: block Lanczos on B*B^T and the dense fallback; both return dependencies among the rows

#include "common.hpp"
#include "matrix.hpp"

// dependencies among rows of B (x with x^TB=0) by block Lanczos on BB^T
// ker(BB^T) is a superset of  ker(B^T):
//  1. combine Lanczos outputs [X+Y | V_m] into vectors y with BB^Ty=0
//  2. combine those y so that B^Ty=0
// returns lists of row indices, all linearly independent.
// B should be filtered first (see filter_matrix).  Tries a few random starts if Lanczos does not converge.
vector<vector<size_t>> lanczos_dependencies(const sparse_matrix& B, uint64_t seed = 12345, int max_attempts = 4);

// Gaussian elimination fallback, same output
vector<vector<size_t>> dense_dependencies(const sparse_matrix& B);
