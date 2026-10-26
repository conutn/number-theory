// linear algebra over F_2: block Lanczos on B*B^T and the dense fallback; both return dependencies among the rows

#include "lanczos.hpp"
#include "numtheory.hpp"

// 64x64 matrices over F_2, r[i] = row i
struct mat64 { uint64_t r[64] = {}; };

static mat64 mat_identity() { mat64 I; for (int i = 0; i < 64; ++i) I.r[i] = uint64_t(1) << i; return I; }

static mat64 mat_xor(const mat64& A, const mat64& B) { mat64 C; for (int i = 0; i < 64; ++i) C.r[i] = A.r[i] ^ B.r[i]; return C; }

static mat64 mat_mask_cols(const mat64& A, uint64_t mask) { mat64 C; for (int i = 0; i < 64; ++i) C.r[i] = A.r[i] & mask; return C; }

static bool mat_is_zero(const mat64& A) { for (int i = 0; i < 64; ++i) if (A.r[i]) return false; return true; }

static mat64 mat_mul(const mat64& A, const mat64& B)
{
	mat64 C;
	for (int i = 0; i < 64; ++i)
	{
		uint64_t x = A.r[i], acc = 0;
		while (x) { acc ^= B.r[ctz64(x)]; x &= x - 1; }
		C.r[i] = acc;
	}
	return C;
}

// V^TU (V, U: n x 64)
static mat64 inner(const vector<uint64_t>& V, const vector<uint64_t>& U)
{
	mat64 R;
	for (size_t i = 0; i < V.size(); ++i)
	{
		uint64_t x = V[i];
		while (x) { R.r[ctz64(x)] ^= U[i]; x &= x - 1; }
	}
	return R;
}

// V*M (n x 64 times 64 x 64)
static vector<uint64_t> vec_mat(const vector<uint64_t>& V, const mat64& M)
{
	vector<uint64_t> out(V.size(), 0);
	for (size_t i = 0; i < V.size(); ++i)
	{
		uint64_t x = V[i], acc = 0;
		while (x) { acc ^= M.r[ctz64(x)]; x &= x - 1; }
		out[i] = acc;
	}
	return out;
}

// find S (bit mask) so T[S,S] is invertible and maximal (rank T = |S|)
// every column in `prefer` (the ones left out by the previous iteration) must be selected for the recurrence to stay valid.
// Gram entries between two such columns vanish by construction, so a preferred column k usually can only
// enter S as part of a 2x2 block {k, j} with some partner j, which may be ANY unselected column (even one with
// a nonzero diagonal).  therefore all preferred columns are handled first, before any other pivot can use
// up their partners; the old code took the first diagonal pivot among all columns and often got stuck.
// returns false if some preferred column could not be selected (breakdown)
// Winv = (T[S,S])^-1 as a 64x64 matrix (zero outside S)
static bool select_cols(const mat64& T, uint64_t prefer, uint64_t& S, mat64& Winv)
{
	mat64 t = T; // Schur complement of T w.r.t. the pivots chosen so far
	S = 0;

	// eliminate a single pivot p (t[p][p] = 1)
	auto pivot1 = [&](int p) {
		const uint64_t rp = t.r[p];
		for (int r = 0; r < 64; ++r)
			if (r != p && (t.r[r] >> p & 1)) t.r[r] ^= rp;
		t.r[p] = 0;
		S |= uint64_t(1) << p;
	};

	// eliminate the 2x2 block {p,q}, t[p][p] = 0, t[p][q] = 1, c = t[q][q] arbitrary; block inverse = [[c,1],[1,0]]
	auto pivot2 = [&](int p, int q) {
		const uint64_t rp = t.r[p], rq = t.r[q];
		const bool c = t.r[q] >> q & 1;
		for (int r = 0; r < 64; ++r)
		{
			if (r == p || r == q) continue;
			const bool x = t.r[r] >> p & 1, y = t.r[r] >> q & 1;
			if ((x && c) != y) t.r[r] ^= rp; // coefficient on row p: x*c + y
			if (x) t.r[r] ^= rq;             // coefficient on row q: x
		}
		t.r[p] = t.r[q] = 0;
		S |= (uint64_t(1) << p) | (uint64_t(1) << q);
	};

	// phase 1: the preferred columns
	bool ok = true;
	for (;;)
	{
		const uint64_t todo = prefer & ~S;
		if (!todo) break;

		int k = ctz64(todo);
		if (t.r[k] >> k & 1) { pivot1(k); continue; }

		// partner: prefer one that must be selected anyway
		const uint64_t row = t.r[k] & ~S & ~(uint64_t(1) << k);
		if (!row) { ok = false; prefer &= ~(uint64_t(1) << k); continue; } // k is in the radical: cannot be selected
		const uint64_t both = row & todo;
		pivot2(k, ctz64(both ? both : row));
	}

	// phase 2: extend to a maximal S (diagonal pivots first, then pairs of zero-diagonal columns)
	for (;;)
	{
		const uint64_t cand = ~S;
		int p = -1, q = -1;
		for (int i = 0; i < 64 && p < 0; ++i)
			if ((cand >> i & 1) && (t.r[i] >> i & 1)) p = i;
		if (p >= 0) { pivot1(p); continue; }

		for (int i = 0; i < 64 && p < 0; ++i)
			if (cand >> i & 1)
				for (int j = i + 1; j < 64; ++j)
					if ((cand >> j & 1) && (t.r[i] >> j & 1)) { p = i; q = j; break; }
		if (p < 0) break;
		pivot2(p, q);
	}

	// invert T[S,S] using Gauss-Jordan on [T[S,S]+I_(not S) | I]
	mat64 a, inv = mat_identity();
	for (int i = 0; i < 64; ++i) a.r[i] = (S >> i & 1) ? (T.r[i] & S) : (uint64_t(1) << i);

	for (int c = 0; c < 64; ++c)
	{
		int piv = -1;
		for (int r = c; r < 64; ++r) if (a.r[r] >> c & 1) { piv = r; break; }
		if (piv < 0) throw logic_error("select_cols: singular submatrix");
		std::swap(a.r[c], a.r[piv]); std::swap(inv.r[c], inv.r[piv]);
		for (int r = 0; r < 64; ++r)
			if (r != c && (a.r[r] >> c & 1)) { a.r[r] ^= a.r[c]; inv.r[r] ^= inv.r[c]; }
	}

	for (int i = 0; i < 64; ++i) Winv.r[i] = (S >> i & 1) ? (inv.r[i] & S) : 0;
	return ok;
}

// block Lanczos for symmetric matrix A=BB^T
// returns false if it did not converge (V_i^T A V_i never became 0 within a generous iteration cap).
// For a well conditioned matrix it takes ~n/62 iterations; far more means the matrix is badly
// conditioned (filter it first) or the random start was unlucky (try another seed)
static int g_lanczos_status = 0; // 1 = converged, 2 = breakdown, 3 = iteration cap (for diagnostics)

static size_t g_lanczos_iters = 0;

static bool block_lanczos(const sparse_matrix& B, vector<uint64_t>& xy, vector<uint64_t>& vm, uint64_t seed)
{
	const size_t n = B.nrows;
	std::mt19937_64 rng(seed);

	vector<uint64_t> Y(n);
	for (auto& y : Y) y = rng();

	const vector<uint64_t> V0 = mul_A(B, Y);
	vector<uint64_t> X(n, 0), v = V0, v1(n, 0), v2(n, 0);
	mat64 Winv1, Winv2, VAV1, C1;
	const mat64 I = mat_identity();
	uint64_t excluded = 0;
	bool converged = false;
	const size_t cap = std::min<size_t>(n + 64, 3 * (n / 60) + 200);

	g_lanczos_status = 3;
	for (size_t it = 0; it < cap; ++it)
	{
		g_lanczos_iters = it;
		const vector<uint64_t> Av = mul_A(B, v);
		const mat64 VAV = inner(v, Av);
		if (mat_is_zero(VAV)) { converged = true; g_lanczos_status = 1; break; }

		uint64_t S;
		mat64 Winv;
		if (!select_cols(VAV, excluded, S, Winv)) { g_lanczos_status = 2; break; }

		// the columns left out last time must all be selected now; if V^T A V is too degenerate for that
		// the recurrence is no longer valid (unlucky random start): give up so the caller can reseed
		if (excluded & ~S) { g_lanczos_status = 2; break; }

		// X+=V_i Winv_i (V_i^T V_0)
		const vector<uint64_t> dx = vec_mat(v, mat_mul(Winv, inner(v, V0)));
		for (size_t i = 0; i < n; ++i) X[i] ^= dx[i];

		const mat64 VAAV = inner(Av, Av);
		const mat64 C = mat_xor(mat_mask_cols(VAAV, S), VAV); // V^T A^2 V S S^T + V^T A V
		const mat64 D = mat_xor(I, mat_mul(Winv, C));
		const mat64 E = mat_mul(Winv1, mat_mask_cols(VAV, S));
		const mat64 F = mat_mul(Winv2, mat_mul(mat_xor(I, mat_mul(VAV1, Winv1)), mat_mask_cols(C1, S)));

		// V_{i+1} = A V_i S S^T + V_i D + V_{i-1} E + V_{i-2} F
		vector<uint64_t> vn = vec_mat(v, D);
		const vector<uint64_t> t1 = vec_mat(v1, E), t2 = vec_mat(v2, F);
		for (size_t i = 0; i < n; ++i) vn[i] ^= (Av[i] & S) ^ t1[i] ^ t2[i];

		v2 = std::move(v1); v1 = std::move(v); v = std::move(vn);
		Winv2 = Winv1; Winv1 = Winv; VAV1 = VAV; C1 = C;
		excluded = ~S;
	}

	xy.resize(n);
	for (size_t i = 0; i < n; ++i) xy[i] = X[i] ^ Y[i];
	vm = std::move(v);
	// a breakdown (status 2) is not fatal: whatever is returned is verified exactly (A z = 0, B^T x = 0)
	// by lanczos_dependencies, the breakdown only limits how many dependencies come out
	return converged || g_lanczos_status == 2;
}

// column j of n x 64 block as bit vector of length n
static bitvec column_of(const vector<uint64_t>& V, int j)
{
	bitvec c((V.size() + 63) / 64, 0);
	for (size_t i = 0; i < V.size(); ++i)
		if (V[i] >> j & 1) c[i / 64] |= uint64_t(1) << (i % 64);
	return c;
}

static bool bitvec_zero(const bitvec& a) { for (uint64_t w : a) if (w) return false; return true; }

// all independent combinations c with XOR_{j in c} v[j] = 0
// returns bit vectors over j
static vector<bitvec> dependencies(const vector<bitvec>& v)
{
	const size_t k = v.size(), kw = (k + 63) / 64;
	struct row { bitvec x, tag; size_t pivot; };
	vector<row> basis;
	vector<bitvec> out;

	for (size_t j = 0; j < k; ++j)
	{
		bitvec x = v[j], tag(kw, 0);
		tag[j / 64] |= uint64_t(1) << (j % 64);

		for (const row& b : basis)
			if (x[b.pivot / 64] >> (b.pivot % 64) & 1)
			{
				for (size_t w = 0; w < x.size(); ++w) x[w] ^= b.x[w];
				for (size_t w = 0; w < kw; ++w) tag[w] ^= b.tag[w];
			}

		size_t piv = 0;
		bool found = false;
		for (size_t w = 0; w < x.size() && !found; ++w)
			if (x[w]) { piv = w * 64 + ctz64(x[w]); found = true; }

		if (found) basis.push_back({ std::move(x), std::move(tag), piv });
		else out.push_back(std::move(tag));
	}

	return out;
}

static vector<size_t> bits_to_indices(const bitvec& a, size_t n)
{
	vector<size_t> idx;
	for (size_t i = 0; i < n; ++i) if (a[i / 64] >> (i % 64) & 1) idx.push_back(i);
	return idx;
}

// XOR of the vectors v[j] for j in combo
static bitvec combine(const vector<bitvec>& v, const bitvec& combo, size_t len)
{
	bitvec r(len, 0);
	for (size_t j = 0; j < v.size(); ++j)
		if (combo[j / 64] >> (j % 64) & 1)
			for (size_t w = 0; w < len; ++w) r[w] ^= v[j][w];
	return r;
}

// indices of a maximal linearly independent subset of v (greedy, in order)
static vector<size_t> independent_indices(const vector<bitvec>& v)
{
	struct row { bitvec x; size_t pivot; };
	vector<row> basis;
	vector<size_t> keep;
	for (size_t j = 0; j < v.size(); ++j)
	{
		bitvec x = v[j];
		for (const row& b : basis)
			if (x[b.pivot / 64] >> (b.pivot % 64) & 1)
				for (size_t w = 0; w < x.size(); ++w) x[w] ^= b.x[w];

		for (size_t w = 0; w < x.size(); ++w)
			if (x[w])
			{
				const size_t piv = w * 64 + (size_t)ctz64(x[w]);
				basis.push_back({ std::move(x), piv });
				keep.push_back(j);
				break;
			}
	}
	return keep;
}

vector<vector<size_t>> lanczos_dependencies(const sparse_matrix& B, uint64_t seed, int max_attempts)
{
	const size_t n = B.nrows, nw = (n + 63) / 64;

	for (int attempt = 0; attempt < max_attempts; ++attempt)
	{
		vector<uint64_t> xy, vm;
		if (!block_lanczos(B, xy, vm, seed + 7919 * attempt))
		{
			cerr << "block Lanczos did not converge (seed " << attempt << ")\n";
			continue;
		}

		vector<bitvec> Z, AZ;
		for (const auto* blk : { &xy, &vm })
		{
			const vector<uint64_t> Ablk = mul_A(B, *blk);
			for (int j = 0; j < 64; ++j)
			{
				bitvec z = column_of(*blk, j);
				if (bitvec_zero(z)) continue;
				Z.push_back(std::move(z));
				AZ.push_back(column_of(Ablk, j));
			}
		}

		// step 1
		vector<bitvec> ys;
		for (const bitvec& c : dependencies(AZ))
		{
			bitvec y = combine(Z, c, nw);
			if (!bitvec_zero(y)) ys.push_back(std::move(y));
		}

		// keep ys independent: Z is usually linearly dependent, which used to give the same y many times
		{
			vector<bitvec> indep;
			for (size_t j : independent_indices(ys)) indep.push_back(ys[j]);
			ys = std::move(indep);
		}

		// step 2
		vector<bitvec> ws;
		for (const bitvec& y : ys) ws.push_back(mul_Bt_bits(B, y));

		vector<vector<size_t>> result;
		for (const bitvec& c : dependencies(ws))
		{
			const bitvec x = combine(ys, c, nw);
			if (bitvec_zero(x)) continue;
			if (!bitvec_zero(mul_Bt_bits(B, x))) continue; // sanity check
			result.push_back(bits_to_indices(x, n));
		}

		if (!result.empty()) return result;
	}

	return {};
}

vector<vector<size_t>> dense_dependencies(const sparse_matrix& B)
{
	vector<bitvec> rows;
	for (size_t i = 0; i < B.nrows; ++i)
	{
		bitvec r((B.ncols + 63) / 64, 0);
		for (uint32_t c : B.rows[i]) r[c / 64] |= uint64_t(1) << (c % 64);
		rows.push_back(std::move(r));
	}

	vector<vector<size_t>> result;
	for (const bitvec& c : dependencies(rows)) result.push_back(bits_to_indices(c, B.nrows));
	return result;
}
