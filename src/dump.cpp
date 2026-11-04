// binary ckpt files for filtered matrix + Lanczos vectors

#include "dump.hpp"

// file type tags
static const u64 MAGIC_MAT = 0x58544d41534e4647ULL;
static const u64 MAGIC_VEC = 0x43455653474e4647ULL;

// binary write/read one value
template <class T> static void put(ostream& o, const T& x) { o.write(reinterpret_cast<const char*>(&x), sizeof x); }
template <class T> static bool get(istream& i, T& x) { i.read(reinterpret_cast<char*>(&x), sizeof x); return (bool)i; }

// layout of file:
//   magic # <> nrows <> ncols <> for each row write (len, col idxs) <> row_map
void save_filtered(const string& path, const filtered_matrix& F)
{
	ofstream o(path, ios::binary);
	if (!o) { cerr << "cannot write " << path << '\n'; return; }
	put<u64>(o, MAGIC_MAT);
	put<u64>(o, F.B.nrows);
	put<u64>(o, F.B.ncols);
	for (size_t i = 0; i < F.B.nrows; ++i)
	{
		const auto& r = F.B.rows[i];
		put<u32>(o, (u32)r.size());
		if (!r.empty()) o.write(reinterpret_cast<const char*>(r.data()), r.size() * sizeof(u32));
	}
	for (size_t i = 0; i < F.B.nrows; ++i) put<u64>(o, F.row_map[i]);
	o.flush();
	if (!o) cerr << "error while writing " << path << '\n';
	else cerr << "saved matrix " << F.B.nrows << " x " << F.B.ncols << " to " << path << '\n';
}

// read layout written by save_filtered, return true on success
bool load_filtered(const string& path, filtered_matrix& F)
{
	ifstream in(path, ios::binary);
	u64 magic = 0, nr = 0, nc = 0;
	if (!in || !get(in, magic) || magic != MAGIC_MAT || !get(in, nr) || !get(in, nc)) return false;

	F = filtered_matrix();
	F.B.nrows = nr;
	F.B.ncols = nc;
	F.B.rows.resize(nr);
	for (u64 i = 0; i < nr; ++i)
	{
		u32 len;
		if (!get(in, len)) return false;
		F.B.rows[i].resize(len);
		if (len) in.read(reinterpret_cast<char*>(F.B.rows[i].data()), (streamsize)len * sizeof(u32));
		if (!in) return false;
		for (u32 c : F.B.rows[i]) if (c >= nc) return false; // corrupt file
	}
	F.row_map.resize(nr);
	for (u64 i = 0; i < nr; ++i)
	{
		u64 v;
		if (!get(in, v)) return false;
		F.row_map[i] = (size_t)v;
	}
	return true;
}

// layout of file:
//   magic # <> |xy| <> |vm| <> entries of xy + vm
void save_vectors(const string& path, const vector<u64>& xy, const vector<u64>& vm)
{
	ofstream o(path, ios::binary);
	if (!o) { cerr << "cannot write " << path << '\n'; return; }
	put<u64>(o, MAGIC_VEC);
	put<u64>(o, xy.size());
	put<u64>(o, vm.size());
	if (!xy.empty()) o.write(reinterpret_cast<const char*>(xy.data()), xy.size() * sizeof(u64));
	if (!vm.empty()) o.write(reinterpret_cast<const char*>(vm.data()), vm.size() * sizeof(u64));
	o.flush();
	if (!o) cerr << "error while writing " << path << '\n';
	else cerr << "saved Lanczos vectors to " << path << '\n';
}

// read layout written by save_vectors, return true on success
bool load_vectors(const string& path, vector<u64>& xy, vector<u64>& vm)
{
	ifstream in(path, ios::binary);
	u64 magic = 0, n1 = 0, n2 = 0;
	if (!in || !get(in, magic) || magic != MAGIC_VEC || !get(in, n1) || !get(in, n2)) return false;
	xy.resize(n1);
	vm.resize(n2);
	if (n1) in.read(reinterpret_cast<char*>(xy.data()), (streamsize)n1 * sizeof(u64));
	if (n2) in.read(reinterpret_cast<char*>(vm.data()), (streamsize)n2 * sizeof(u64));
	return (bool)in;
}
