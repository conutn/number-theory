#pragma once
// load/dump filtered matrix + block Lanczos output
// checkpointing so post-processing can be run without resieving or redoing Lanczos

#include "common.hpp"
#include "matrix.hpp"

// writes B (rows, ncols) and row_map of F
void save_filtered(const string& path, const filtered_matrix& F);
// reads F, returns false if file missing | truncated | corrupt
bool load_filtered(const string& path, filtered_matrix& F);

// the two n x 64 blocks block_lanczos returns
// write xy + vm to path
void save_vectors(const string& path, const vector<u64>& xy, const vector<u64>& vm);
// reads xy + vm, returns false if file missing | truncated | corrupt
bool load_vectors(const string& path, vector<u64>& xy, vector<u64>& vm);
