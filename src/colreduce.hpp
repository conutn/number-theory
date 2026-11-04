#pragma once
// removes cols that are linear combinations of other cols

#include "matrix.hpp"

// removes cols that are sums of other cols in place
//   one column per cycle of the graph formed by weight-2 columns
//   every copy of a column w/ same support as earlier one
// rows unchanged (only column idxs are renumbered)
// returns number of columns removed, B.ncols is updated
size_t remove_dependent_columns(sparse_matrix& B);
