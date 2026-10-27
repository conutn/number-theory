// standalone linear algebra post-processing

#include "common.hpp"
#include "matrix.hpp"
#include "lanczos.hpp"
#include "dump.hpp"
#include "colreduce.hpp"

// prints column statistics of B to stderr: weight histogram, cycles among weight-2 columns, duplicate columns
// ends with a lower bound on dim ker B, the number of dependent columns (Lanczos on B B^T needs it to stay well below 64 per run)
static void print_stats(const sparse_matrix& B)
{
	const size_t nr = B.nrows, nc = B.ncols;
	vector<uint32_t> w(nc, 0), first(nc, 0), second(nc, 0); // weight and first two rows of each column
	vector<uint64_t> h1(nc, 0x9e3779b97f4a7c15ULL), h2(nc, 0xc2b2ae3d27d4eb4fULL); // two independent hashes of each column's support
	size_t nnz = 0;
	for (size_t r = 0; r < nr; ++r)
	{
		nnz += B.rows[r].size();
		for (uint32_t c : B.rows[r])
		{
			if (w[c] == 0) first[c] = (uint32_t)r; else if (w[c] == 1) second[c] = (uint32_t)r;
			++w[c];
			h1[c] = (h1[c] ^ (r + 1)) * 0x100000001b3ULL;
			h2[c] = (h2[c] + (r + 1)) * 0xff51afd7ed558ccdULL;
			h2[c] ^= h2[c] >> 29;
		}
	}

	size_t hist[6] = {};
	for (size_t c = 0; c < nc; ++c) ++hist[std::min<uint32_t>(w[c], 5)];
	cerr << "stats: " << nr << " rows, " << nc << " columns, " << nnz << " nonzeros (" << (double)nnz / nr << " per row)\n";
	cerr << "stats: column weights: 0:" << hist[0] << " 1:" << hist[1] << " 2:" << hist[2] << " 3:" << hist[3] << " 4:" << hist[4] << " >=5:" << hist[5] << '\n';

	// a weight-2 column is an edge between two rows, every cycle of that graph is a set of columns summing to 0
	vector<uint32_t> parent(nr);
	std::iota(parent.begin(), parent.end(), 0);
	auto find = [&](uint32_t x) { while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; }; // union-find with path halving
	size_t edges = 0, cycles = 0;
	for (size_t c = 0; c < nc; ++c)
		if (w[c] == 2)
		{
			++edges;
			const uint32_t a = find(first[c]), b = find(second[c]);
			if (a == b) ++cycles; else parent[a] = b; // edge inside a component closes a cycle
		}

	// two columns with the same support (weight >= 3) give e_j + e_j' in ker B
	vector<std::pair<std::pair<uint64_t, uint64_t>, uint32_t>> keys;
	for (size_t c = 0; c < nc; ++c)
		if (w[c] >= 3) keys.push_back({ { h1[c], h2[c] }, (uint32_t)c });
	std::sort(keys.begin(), keys.end());
	size_t dups = 0;
	for (size_t i = 1; i < keys.size(); ++i)
		if (keys[i].first == keys[i - 1].first) ++dups; // equal hashes taken as equal supports

	cerr << "stats: weight-2 columns: " << edges << ", independent cycles among them: " << cycles << " (each is a column dependency)\n";
	cerr << "stats: duplicate columns of weight >= 3: " << dups << '\n';
	cerr << "stats: dim ker B >= " << cycles + dups << "  (Lanczos on B B^T needs this to stay well below 64 per run)\n";
}

// usage: ./postprocess [options] run.mat [run.vec0 run.vec1 ...]
// run.mat is written by save_filtered
// run.vecN are written by lanczos_dependencies when GNFS_DUMP_PREFIX=run is set
// several vec files are pooled: every Lanczos run gives another random slice of ker A
// the vectors only fit the matrix they were computed with: do not combine --reduce with old vec files
// options (--stats and --reduce run in the order given, before anything else):
//   --stats             column statistics
//   --reduce            remove columns that are combinations of other columns (colreduce.cpp)
//   --save-reduced f    write the reduced matrix to f, so that vec files from --solve can be reused later
//   --solve N           run block Lanczos, up to N attempts, pooling the kernel vectors
//   -o file             write one dependency per line, as row indices into `mat` in main.cpp (before filtering)
//   --raw               with -o: write indices into the filtered matrix instead
// example: ./postprocess --stats --reduce --stats --solve 2 -o deps.txt run2.mat
// returns 0 on success, 1 if no dependencies were found, 2 on bad input
int main(int argc, char** argv)
{
	bool raw = false, stats = false, reduce = false;
	int solve = 0;
	std::string out_path, reduced_path;
	vector<std::string> args;
	vector<int> order; // 1 = --stats, 2 = --reduce, in the order given, so that "--stats --reduce --stats" works
	for (int i = 1; i < argc; ++i)
	{
		const std::string a = argv[i];
		if (a == "--raw") raw = true;
		else if (a == "--stats") { stats = true; order.push_back(1); }
		else if (a == "--reduce") { reduce = true; order.push_back(2); }
		else if (a == "--save-reduced" && i + 1 < argc) reduced_path = argv[++i];
		else if (a == "--solve" && i + 1 < argc) solve = atoi(argv[++i]);
		else if (a == "-o" && i + 1 < argc) out_path = argv[++i];
		else args.push_back(a);
	}
	if (args.empty())
	{
		cerr << "usage: postprocess [--stats] [--reduce] [--save-reduced f] [--solve N] [--raw] [-o deps.txt] matrix.bin [vectors.bin ...]\n";
		return 2;
	}

	const auto t0 = std::chrono::steady_clock::now();

	filtered_matrix F;
	if (!load_filtered(args[0], F)) { cerr << "cannot read matrix " << args[0] << '\n'; return 2; }
	cerr << "matrix: " << F.B.nrows << " rows x " << F.B.ncols << " columns\n";

	for (int step : order)
	{
		if (step == 1) print_stats(F.B);
		else remove_dependent_columns(F.B);
	}
	(void)stats; (void)reduce; // only the order matters
	if (!reduced_path.empty()) save_filtered(reduced_path, F);

	vector<vector<size_t>> deps;

	if (solve > 0)
	{
		deps = lanczos_dependencies(F.B, 1, solve); // pools the kernel vectors of all attempts
	}
	else
	{
		vector<bitvec> pool;
		for (size_t k = 1; k < args.size(); ++k)
		{
			vector<uint64_t> xy, vm;
			if (!load_vectors(args[k], xy, vm)) { cerr << "cannot read vectors " << args[k] << '\n'; return 2; }
			cerr << "--- " << args[k] << '\n';
			for (bitvec& y : lanczos_kernel_vectors(F.B, xy, vm)) pool.push_back(std::move(y));
			deps = dependencies_from_pool(F.B, pool);
			if (!deps.empty()) break; // enough already; drop the break to see the numbers for every file
		}
		if (args.size() < 2) return 0; // only --stats / --reduce
	}

	const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
	cerr << deps.size() << " dependencies in " << secs << " s\n";
	if (deps.empty()) return 1;

	if (!out_path.empty())
	{
		std::ofstream out(out_path);
		for (const auto& d : deps)
		{
			for (size_t k = 0; k < d.size(); ++k) out << (k ? " " : "") << (raw ? d[k] : F.row_map[d[k]]);
			out << '\n';
		}
		cerr << "written to " << out_path << '\n';
	}
	return 0;
}
