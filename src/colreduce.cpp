// removal of dependent columns

#include "colreduce.hpp"

// removes dependent columns of B (see colreduce.hpp), renumbers the remaining ones
// returns number of columns removed
size_t remove_dependent_columns(sparse_matrix& B)
{
	const size_t nr = B.nrows, nc = B.ncols;
	vector<uint32_t> w(nc, 0), first(nc, 0), second(nc, 0); // column weight, and its first two rows
	// two independent 64-bit hashes of the support of each column: a false match is ~2^-128 likely
	vector<uint64_t> h1(nc, 0x9e3779b97f4a7c15ULL), h2(nc, 0xc2b2ae3d27d4eb4fULL);

	for (size_t r = 0; r < nr; ++r)
		for (uint32_t c : B.rows[r])
		{
			if (w[c] == 0) first[c] = (uint32_t)r; else if (w[c] == 1) second[c] = (uint32_t)r;
			++w[c];
			h1[c] = (h1[c] ^ (r + 1)) * 0x100000001b3ULL;
			h2[c] = (h2[c] + (r + 1)) * 0xff51afd7ed558ccdULL;
			h2[c] ^= h2[c] >> 29;
		}

	vector<char> drop(nc, 0);

	// weight-2 column = edge between two rows, the edge closing a cycle is the sum of the others
	vector<uint32_t> parent(nr);
	std::iota(parent.begin(), parent.end(), 0);
	auto find = [&](uint32_t x) { while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; }; // path halving
	size_t cycles = 0;
	for (size_t c = 0; c < nc; ++c)
		if (w[c] == 2)
		{
			const uint32_t a = find(first[c]), b = find(second[c]);
			if (a == b) { drop[c] = 1; ++cycles; } else parent[a] = b;
		}

	// identical columns of weight >= 3 have equal hashes, sorting puts them next to each other (smallest index first)
	vector<std::pair<std::pair<uint64_t, uint64_t>, uint32_t>> keys;
	for (size_t c = 0; c < nc; ++c)
		if (w[c] >= 3) keys.push_back({ { h1[c], h2[c] }, (uint32_t)c });
	std::sort(keys.begin(), keys.end());
	size_t dups = 0;
	for (size_t i = 1; i < keys.size(); ++i)
		if (keys[i].first == keys[i - 1].first) { drop[keys[i].second] = 1; ++dups; }

	const size_t removed = cycles + dups;
	if (removed == 0) return 0;

	// renumber the surviving columns and compact each row in place
	vector<uint32_t> newid(nc, 0);
	uint32_t k = 0;
	for (size_t c = 0; c < nc; ++c) if (!drop[c]) newid[c] = k++;
	for (auto& row : B.rows)
	{
		size_t o = 0;
		for (uint32_t c : row) if (!drop[c]) row[o++] = newid[c];
		row.resize(o);
	}
	B.ncols = k;
	cerr << "removed " << removed << " dependent columns (" << cycles << " from cycles of weight-2 columns, " << dups << " duplicates): "
	     << nc << " -> " << B.ncols << " columns\n";
	return removed;
}
