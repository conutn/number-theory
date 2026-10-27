#pragma once
// removal of columns that are linear combinations of other columns

#include "matrix.hpp"

// removes columns that are sums of other columns, in place:
//  - one column per independent cycle of the graph formed by the weight-2 columns
//  - every extra copy of a column with the same support as an earlier one
// rows are unchanged (only their column indices are renumbered), so row dependencies of the result
// are exactly those of the input and row_map / expand() need no change
// returns number of columns removed, B.ncols is updated
size_t remove_dependent_columns(sparse_matrix& B);
