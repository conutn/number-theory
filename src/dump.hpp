#pragma once
// load/dump filtered matrix and block Lanczos output
// checkpointing so post-processing can be run without resieving or redoing Lanczos

#include "common.hpp"
#include "matrix.hpp"

// writes B (rows, ncols) and row_map of F to path, errors are reported on cerr
void save_filtered(const std::string& path, const filtered_matrix& F);
// reads F back from path, returns false if the file is missing, truncated or corrupt
bool load_filtered(const std::string& path, filtered_matrix& F);

// the two n x 64 blocks returned by block_lanczos
// save writes xy and vm to path, errors are reported on cerr
void save_vectors(const std::string& path, const vector<u64>& xy, const vector<u64>& vm);
// returns false if the file is missing, truncated or has the wrong type
bool load_vectors(const std::string& path, vector<u64>& xy, vector<u64>& vm);
