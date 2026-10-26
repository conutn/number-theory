// sparse matrix over F_2: construction, singleton/surplus-row filtering, products with B and B^T

#include "matrix.hpp"

sparse_matrix build_matrix(const vector<nfs_relation>& rels, size_t ncols)
{
	sparse_matrix B;
	B.nrows = rels.size();
	B.ncols = ncols;
	for (const nfs_relation& r : rels) B.rows.push_back(r.cols);
	return B;
}

filtered_matrix filter_matrix(const vector<nfs_relation>& rels, size_t ncols, long keep_excess)
{
	vector<vector<uint32_t>> col_rows(ncols);
	vector<uint32_t> w(ncols, 0);
	for (size_t i = 0; i < rels.size(); ++i)
		for (uint32_t c : rels[i].cols) { col_rows[c].push_back((uint32_t)i); ++w[c]; }

	vector<char> alive(rels.size(), 1);
	vector<uint32_t> todo;
	for (size_t c = 0; c < ncols; ++c) if (w[c] == 1) todo.push_back((uint32_t)c);

	auto propagate = [&]() {
		while (!todo.empty())
		{
			const uint32_t c = todo.back();
			todo.pop_back();
			if (w[c] != 1) continue;
			for (uint32_t i : col_rows[c])
			{
				if (!alive[i]) continue;
				alive[i] = 0;
				for (uint32_t c2 : rels[i].cols)
					if (--w[c2] == 1) todo.push_back(c2);
				break;
			}
		}
	};
	propagate();

	while (keep_excess >= 0)
	{
		vector<uint32_t> rows;
		for (size_t i = 0; i < rels.size(); ++i) if (alive[i]) rows.push_back((uint32_t)i);
		size_t cols = 0;
		for (size_t c = 0; c < ncols; ++c) if (w[c] > 0) ++cols;
		const long excess = (long)rows.size() - (long)cols;
		if (excess <= keep_excess) break;

		const size_t kill = (size_t)std::max<long>(1, (excess - keep_excess) / 2);
		std::sort(rows.begin(), rows.end(), [&](uint32_t x, uint32_t y) {
			return rels[x].cols.size() != rels[y].cols.size() ? rels[x].cols.size() > rels[y].cols.size() : x < y; });
		for (size_t k = 0; k < kill; ++k)
		{
			alive[rows[k]] = 0;
			for (uint32_t c2 : rels[rows[k]].cols) if (--w[c2] == 1) todo.push_back(c2);
		}
		propagate();
	}

	vector<int64_t> newcol(ncols, -1);
	filtered_matrix F;
	for (size_t c = 0; c < ncols; ++c) if (w[c] > 0) newcol[c] = (int64_t)F.B.ncols++;
	for (size_t i = 0; i < rels.size(); ++i)
	{
		if (!alive[i]) continue;
		vector<uint32_t> row;
		for (uint32_t c : rels[i].cols) row.push_back((uint32_t)newcol[c]);
		F.B.rows.push_back(std::move(row));
		F.row_map.push_back(i);
	}
	F.B.nrows = F.B.rows.size();
	return F;
}

// out=B^Tv (ncols words)
// requires:
//  v.size() same as B.nrows
static vector<uint64_t> mul_Bt(const sparse_matrix& B, const vector<uint64_t>& v)
{
	vector<uint64_t> out(B.ncols, 0);
	for (size_t i = 0; i < B.nrows; ++i)
		for (uint32_t c : B.rows[i]) out[c] ^= v[i];
	return out;
}

// out=Bu (nrows words)
// requires:
//  u.size() same as B.ncols
static vector<uint64_t> mul_B(const sparse_matrix& B, const vector<uint64_t>& u)
{
	vector<uint64_t> out(B.nrows, 0);
	for (size_t i = 0; i < B.nrows; ++i)
		for (uint32_t c : B.rows[i]) out[i] ^= u[c];
	return out;
}

vector<uint64_t> mul_A(const sparse_matrix& B, const vector<uint64_t>& v)
{
	return mul_B(B, mul_Bt(B, v));
}

bitvec mul_Bt_bits(const sparse_matrix& B, const bitvec& y)
{
	bitvec out((B.ncols + 63) / 64, 0);
	for (size_t i = 0; i < B.nrows; ++i)
		if (y[i / 64] >> (i % 64) & 1)
			for (uint32_t c : B.rows[i]) out[c / 64] ^= uint64_t(1) << (c % 64);
	return out;
}
