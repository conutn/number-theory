// relation finding: line sieve on both sides, candidate factoring, quadratic chars,
//
// segmented line sieve over a in [-A,A] for b = 1,2,3...

#include "sieve.hpp"
#include "numtheory.hpp"

// return true if cofactor v is 1 | large prime <= lp_bound, set lp to it (0 if v = 1)
// if v <= lp_bound <= B^2 it must be prime | 1
static bool take_cofactor(const bigint& v, uint64_t lp_bound, uint64_t& lp)
{
	lp = 0;
	if (v == 1) return true;
	if (lp_bound && v <= bigint(lp_bound))
	{
		lp = v.convert_to<uint64_t>();
		return true;
	}
	return false;
}


// check if a+b*m is smooth w/ at most one large prime
// primes dividing value are already known
// hits = ascending indices from base
static bool factor_rational_hits(int64_t a, int64_t b, const bigint& m, const factor_bases& fb, const vector<uint32_t>& hits, uint64_t lp_bound, uint64_t& lp, vector<uint32_t>& cols, vector<fac_entry>& fac)
{
	bigint v = bigint(a) + bigint(b) * m;
	if (v < 0)
	{
		cols.push_back(0);
		v = -v;
	}

	for (uint32_t i : hits)
	{
		const bigint P = fb.rational[i].p;
		unsigned e = 0;
		while (v % P == 0)
		{
			v /= P;
			++e;
		}
		if (e) fac.push_back({ i, e });
		if (e & 1) cols.push_back(uint32_t(1 + i));
	}

	return take_cofactor(v, lp_bound, lp);
}

static bool factor_algebraic_hits(const poly& f, int64_t a, int64_t b, const factor_bases& fb, size_t offset, const vector<uint32_t>& hits, uint64_t lp_bound, uint64_t& lp, vector<uint32_t>& cols, vector<fac_entry>& fac)
{
	bigint N = abs(norm(f, a, b));
	if (N == 0) return false;

	for (uint32_t i : hits)
	{
		const bigint P = fb.algebraic[i].p;
		unsigned e = 0;
		while (N % P == 0)
		{
			N /= P;
			++e;
		}
		if (e) fac.push_back({ i, e });
		if (e & 1) cols.push_back(uint32_t(offset + i));
	}

	return take_cofactor(N, lp_bound, lp);
}

// quadratic chars
// update cols if ((a+b*s)/q)=-1
// return false if a+b*s = 0 mod q
static bool quadratic_chars(int64_t a, int64_t b, const factor_bases& fb, size_t offset, vector<uint32_t>& cols)
{
	for (size_t k = 0; k < fb.quadratic.size(); ++k)
	{
		const uint64_t q = fb.quadratic[k].p;
		const uint64_t v = lin_mod(a, b, q, fb.quadratic[k].r);
		if (v == 0) return false;
		if (jacobi(v, q) == -1) cols.push_back(uint32_t(offset + k));
	}
	return true;
}

static uint64_t mulmod_u64(uint64_t a, uint64_t b, uint64_t m)
{
	return (uint64_t)((unsigned __int128)a * b % m);
}

// return a^-1 mod p (a != 0 mod p) by extended Euclid
static uint64_t inv_mod_u64(uint64_t a, uint64_t p)
{
	__int128 r0 = p, r1 = a % p, t0 = 0, t1 = 1;
	while (r1 != 0)
	{
		const __int128 q = r0 / r1;
		const __int128 r2 = r0 - q * r1;
		r0 = r1;
		r1 = r2;
		const __int128 t2 = t0 - q * t1;
		t0 = t1;
		t1 = t2;
	}
	if (t0 < 0) t0 += p;
	return (uint64_t)t0;
}

// set r to root of deg 1 algebraic ideal above p dividing a+b*theta (r = -a/b mod p)
// return false if p | b
static bool alg_root(int64_t a, int64_t b, uint64_t p, uint64_t& r)
{
	const uint64_t bm = uint64_t(b) % p;
	if (bm == 0) return false;
	int64_t am = a % (int64_t)p;
	if (am < 0) am += (int64_t)p;
	r = (p - mulmod_u64((uint64_t)am, inv_mod_u64(bm, p), p)) % p;
	return true;
}

// large-prime graph
// vertex 0 = no large prime, other vertices: rational p or algebraic ideal (p,r)
// partial relation is edge between two vertices
struct lp_graph
{
	map<tuple<int, uint64_t, uint64_t>, uint32_t> id; // (type, p, r) -> vertex
	vector<array<uint32_t, 2>> edges;

	explicit lp_graph(const vector<nfs_relation>& rels)
	{
		auto vid = [&](int side, uint64_t p, uint64_t r) -> uint32_t
		{
			if (p == 0) return 0;
			auto it = id.find({ side, p, r });
			if (it != id.end()) return it->second;
			const uint32_t v = uint32_t(id.size() + 1);
			id.emplace(make_tuple(side, p, r), v);
			return v;
		};
		for (const nfs_relation& x : rels)
			edges.push_back({ vid(0, x.lp_r, 0), vid(1, x.lp_a, x.lp_a_root) });
	}

	size_t nvert() const
	{
		return id.size() + 1;
	}
};

// union-find w/ path halving
static uint32_t uf_find(vector<uint32_t>& p, uint32_t x)
{
	while (p[x] != x)
	{
		p[x] = p[p[x]];
		x = p[x];
	}
	return x;
}

// # of independent cycles = # of full relations obtainable by combining partials
size_t count_lp_cycles(const vector<nfs_relation>& partials)
{
	if (partials.empty()) return 0;
	lp_graph g(partials);
	vector<uint32_t> uf(g.nvert());
	for (uint32_t i = 0; i < uf.size(); ++i) uf[i] = i;
	size_t cycles = 0;
	for (const auto& e : g.edges)
	{
		const uint32_t u = uf_find(uf, e[0]), v = uf_find(uf, e[1]);
		if (u == v) ++cycles; // edge closes a cycle
		else uf[u] = v;
	}
	return cycles;
}

// incremental count_lp_cycles: partial relations added online, cycle count updated
// vertices numbered as in lp_graph (0 = no large prime)
struct cycle_counter
{
	struct key
	{
		uint64_t p, r; // large prime and (algebraic side) root, r = 0 on rational side
		int side;
		bool operator==(const key& o) const { return p == o.p && r == o.r && side == o.side; }
	};
	struct key_hash
	{
		size_t operator()(const key& k) const
		{
			uint64_t h = k.p * 0x9E3779B97F4A7C15ULL ^ (k.r + 0x7F4A7C15ULL) * 0xC2B2AE3D27D4EB4FULL ^ (uint64_t)k.side;
			return (size_t)(h ^ (h >> 29));
		}
	};

	unordered_map<key, uint32_t, key_hash> id; // (side, p, r) -> vertex
	vector<uint32_t> uf{ 0 }; // union-find parents, vertex 0 exists from the start
	size_t cycles = 0; // independent cycles so far

	uint32_t vertex(int side, uint64_t p, uint64_t r)
	{
		if (p == 0) return 0;
		auto it = id.find({ p, r, side });
		if (it != id.end()) return it->second;
		const uint32_t v = (uint32_t)uf.size();
		uf.push_back(v);
		id.emplace(key{ p, r, side }, v);
		return v;
	}

	// adds edge of partial relation, return true if closed a cycle
	bool add(const nfs_relation& x)
	{
		const uint32_t u = uf_find(uf, vertex(0, x.lp_r, 0)), v = uf_find(uf, vertex(1, x.lp_a, x.lp_a_root));
		if (u == v)
		{
			++cycles;
			return true;
		}
		uf[u] = v;
		return false;
	}
};

// add one combined relation per cycle of lp graph to out
void combine_partials(const vector<nfs_relation>& partials, vector<nfs_relation>& out)
{
	if (partials.empty()) return;
	lp_graph g(partials);
	const size_t V = g.nvert(), E = g.edges.size();

	vector<vector<uint32_t>> adj(V); // edge indices
	for (size_t e = 0; e < E; ++e)
	{
		adj[g.edges[e][0]].push_back((uint32_t)e);
		adj[g.edges[e][1]].push_back((uint32_t)e);
	}

	const uint32_t NONE = UINT32_MAX;
	vector<uint32_t> pv(V, NONE), pe(V, NONE), depth(V, 0); // parent vertex, parent edge, and depth in the spanning forest
	vector<char> seen(V, 0), tree(E, 0); // tree[e] = edge is in forest

	for (uint32_t root = 0; root < V; ++root) // BFS spanning forest
	{
		if (seen[root]) continue;
		seen[root] = 1;
		vector<uint32_t> q{ root };
		for (size_t h = 0; h < q.size(); ++h)
		{
			const uint32_t u = q[h];
			for (uint32_t e : adj[u])
			{
				const uint32_t v = g.edges[e][0] == u ? g.edges[e][1] : g.edges[e][0];
				if (seen[v]) continue;
				seen[v] = 1;
				pv[v] = u;
				pe[v] = e;
				depth[v] = depth[u] + 1;
				tree[e] = 1;
				q.push_back(v);
			}
		}
	}

	// every non-tree edge closes a cycle: edge + tree paths from ends up to ancestor
	for (size_t e0 = 0; e0 < E; ++e0)
	{
		if (tree[e0]) continue;
		vector<uint32_t> used{ (uint32_t)e0 };
		uint32_t u = g.edges[e0][0], v = g.edges[e0][1];
		while (u != v)
		{
			if (depth[u] >= depth[v])
			{
				used.push_back(pe[u]);
				u = pv[u];
			}
			else
			{
				used.push_back(pe[v]);
				v = pv[v];
			}
		}

		nfs_relation c;
		c.a = partials[used[0]].a;
		c.b = partials[used[0]].b;
		vector<uint32_t> cols;
		for (uint32_t e : used)
			cols.insert(cols.end(), partials[e].cols.begin(), partials[e].cols.end());
		sort(cols.begin(), cols.end());
		for (size_t i = 0; i < cols.size();) // keep cols w/ odd multiplicity
		{
			size_t j = i;
			while (j < cols.size() && cols[j] == cols[i]) ++j;
			if ((j - i) & 1) c.cols.push_back(cols[i]);
			i = j;
		}
		c.members = move(used);
		out.push_back(move(c));
	}
}

// approximate log(x) for x >= 1
// x = 2^e * (1+t), t in [0,1), polynomial given by equioscillation theorem
static inline double fast_log(double x)
{
	uint64_t bits;
	memcpy(&bits, &x, sizeof bits);
	const int e = int((bits >> 52) & 0x7ff) - 1023; // exponent
	bits = (bits & 0x000fffffffffffffULL) | 0x3ff0000000000000ULL; // mantissa with exponent 0: 1+t
	double m;
	memcpy(&m, &bits, sizeof m);
	const double t = m - 1.0;
	return ((double)e + 0.0006+t*(1.4189+t*(-0.5771+t*0.1582))) * 0.6931471805599453;
}

// segmented, multithreaded line sieve
//
// every line b, interval a in [-A,A] (index i = a + A) is split into segments of size SEG
//   sieve values are int16 in units of 1/SC nat (val = log(not yet sieved norm))
//   primes w/ p^k < SEG ("small") walk segments, position survives across segments
//   primes w/ p^k >= SEG ("large") hit segments <=1 times, hits bucketed per segment
//   algebraic log-norm evaluated every NORM_STRIDE positions, interpolated in between
//     except when not smooth (sign change | big jump) every position is evaluated
//   lines assigned in chunks of consecutive b to threads

namespace
{

constexpr double SC = 16.0; // sieve values are in units of 1/SC nat
#ifndef SIEVE_SEG_LOG
#define SIEVE_SEG_LOG 16
#endif
#ifndef SIEVE_SEGS_PER_SUPER
#define SIEVE_SEGS_PER_SUPER 64
#endif
#ifndef SIEVE_TRIAL_LIMIT
#define SIEVE_TRIAL_LIMIT 1024
#endif
constexpr int SEG_LOG = SIEVE_SEG_LOG; // segment = 2^SEG_LOG positions, offsets must fit in uint16_t
constexpr uint32_t SEG = 1u << SEG_LOG;
constexpr uint32_t SEG_MASK = SEG - 1;
constexpr uint32_t SEGS_PER_SUPER = SIEVE_SEGS_PER_SUPER; // segments per bucket round
constexpr uint32_t TRIAL_LIMIT = SIEVE_TRIAL_LIMIT; // p^k below this: find divisors of candidate by trial, above: re-walk the segment
constexpr uint32_t NORM_STRIDE = 16; // norm evaluation every NORM_STRIDE positions
constexpr float NORM_JUMP = float(0.35 * SC); // large log change between two samples: evaluate exactly
constexpr int THRESH_MARGIN = int(0.25 * SC); // extra room for rounding of prime logs
constexpr int LINES_PER_CHUNK = 8; // consecutive b per task
constexpr uint32_t NO_IDX = UINT32_MAX; // table entry is prime power
constexpr uint32_t SENT = UINT32_MAX; // position of prime not participating in this line

static_assert(SEG_LOG <= 16, "bucket offsets are uint16");

// one sieve entry
// prime power pk of prime p w/ root r (a = -b*r mod pk)
// index idx of p in the factor base
struct ent
{
	uint32_t pk, p, r, idx;
	int16_t lg; // SC*log(p) rounded
};

// one bucket hit
// offset inside segment
// log to subtract
// idx of prime (NO_IDX for prime powers)
struct bent
{
	uint32_t idx;
	uint16_t off;
	int16_t lg;
};

// entries of one side split by size
// p is non-decreasing inside each list (needed for p | b test)
struct side_tab
{
	vector<ent> small, large;
};

// round to the nearest int for |x| < 16384
static inline int16_t rnd16(float x)
{
	return (int16_t)((int32_t)(x + 16384.5f) - 16384);
}

// builds entry table from list of entries
// requires:
//   all entries in factor base order
// returns entries split into small + large
static side_tab make_side_tab(const vector<ent>& all)
{
	side_tab t;
	for (const ent& e : all) (e.pk < SEG ? t.small : t.large).push_back(e);
	return t;
}

// per thread state of one side of sieve
struct side_work
{
	const side_tab* tab = nullptr;
	vector<uint32_t> soff, spos, spos0; // small: first hit index at the start of the line, running position, position at the start of the segment
	vector<uint32_t> loff, lpos; // large: the same (no per-segment copy needed)
	vector<vector<bent>> bucket; // large hits of the current super segment, per segment
	vector<vector<uint32_t>> hits; // per candidate of the current segment: factor base indices dividing it
	vector<int16_t> arr; // sieve values of the current segment

	// sizes all arrays for the given table
	void init(const side_tab& t)
	{
		tab = &t;
		soff.assign(t.small.size(), 0);
		spos.assign(t.small.size(), 0);
		spos0.assign(t.small.size(), 0);
		loff.assign(t.large.size(), 0);
		lpos.assign(t.large.size(), 0);
		bucket.assign(SEGS_PER_SUPER, {});
		arr.assign(SEG, 0);
	}

	// first index of the line (a + A) hit by each entry, computed from scratch for line b
	void init_offsets(uint64_t b, int64_t A)
	{
		auto first = [&](const ent& e) -> uint32_t {
			const uint64_t a0 = (e.pk - (b % e.pk) * e.r % e.pk) % e.pk; // a = -b*r mod pk
			return (uint32_t)((a0 + (uint64_t)A) % e.pk);
		};
		for (size_t t = 0; t < soff.size(); ++t) soff[t] = first(tab->small[t]);
		for (size_t t = 0; t < loff.size(); ++t) loff[t] = first(tab->large[t]);
	}

	// moves the offsets from line b to line b+1: a0 changes by -r, so the index does too (mod pk)
	void advance_offsets()
	{
		for (size_t t = 0; t < soff.size(); ++t)
		{
			const ent& e = tab->small[t];
			const uint32_t o = soff[t];
			soff[t] = o >= e.r ? o - e.r : o + e.pk - e.r;
		}
		for (size_t t = 0; t < loff.size(); ++t)
		{
			const ent& e = tab->large[t];
			const uint32_t o = loff[t];
			loff[t] = o >= e.r ? o - e.r : o + e.pk - e.r;
		}
	}

	// sets the running positions for line b, primes dividing b do not sieve (a+b*r = a mod p and gcd(a,b)=1)
	void start_line(uint64_t b)
	{
		spos = soff;
		lpos = loff;
		for (size_t t = 0; t < spos.size(); ++t)
		{
			if (tab->small[t].p > b) break; // p non-decreasing: no later entry can divide b
			if (b % tab->small[t].p == 0) spos[t] = SENT;
		}
		for (size_t t = 0; t < lpos.size(); ++t)
		{
			if (tab->large[t].p > b) break;
			if (b % tab->large[t].p == 0) lpos[t] = SENT;
		}
	}

	// bucket sort the hits of all large entries falling in [sb, se) (sb is a multiple of SEG)
	void fill_buckets(uint32_t sb, uint32_t se)
	{
		for (auto& v : bucket) v.clear();
		for (size_t t = 0; t < lpos.size(); ++t)
		{
			uint32_t x = lpos[t];
			if (x >= se) continue; // also skips SENT
			const ent& e = tab->large[t];
			for (; x < se; x += e.pk) bucket[(x - sb) >> SEG_LOG].push_back({ e.idx, (uint16_t)(x & SEG_MASK), e.lg });
			lpos[t] = x;
		}
	}

	// subtracts the logs of all small entries over the segment [base, base+L)
	void sieve_small(uint32_t base, uint32_t L)
	{
		const uint32_t end = base + L;
		int16_t* a = arr.data() - base; // a[x] = value at index x
		for (size_t t = 0; t < spos.size(); ++t)
		{
			const ent& e = tab->small[t];
			uint32_t x = spos[t];
			for (; x < end; x += e.pk) a[x] -= e.lg;
			spos[t] = x; // continues in the next segment
		}
	}

	// subtracts the bucketed hits of segment number s of the current super segment
	void apply_bucket(uint32_t s)
	{
		int16_t* a = arr.data();
		for (const bent& h : bucket[s]) a[h.off] -= h.lg;
	}

	// appends to hits[k] the indices of the primes dividing candidate k, for candidates at segment offsets js
	// bm marks the candidate positions (segment offsets), slot maps a marked offset to its candidate number
	void collect_hits(uint32_t base, uint32_t L, uint32_t s, const vector<uint32_t>& js, const vector<uint64_t>& bm, const vector<uint32_t>& slot)
	{
		const size_t nc = js.size();
		for (size_t t = 0; t < spos0.size(); ++t)
		{
			const ent& e = tab->small[t];
			if (e.idx == NO_IDX || spos0[t] == SENT) continue;
			if (e.pk < TRIAL_LIMIT)
			{
				const uint32_t o = soff[t];
				for (size_t k = 0; k < nc; ++k)
					if ((base + js[k]) % e.pk == o) hits[k].push_back(e.idx); // candidate lies on the progression of this prime
			}
			else
			{
				const uint32_t end = base + L;
				for (uint32_t x = spos0[t]; x < end; x += e.pk)
				{
					const uint32_t j = x - base;
					if (bm[j >> 6] >> (j & 63) & 1) hits[slot[j]].push_back(e.idx);
				}
			}
		}
		for (const bent& h : bucket[s])
			if (h.idx != NO_IDX && (bm[h.off >> 6] >> (h.off & 63) & 1)) hits[slot[h.off]].push_back(h.idx);
	}
};

// read-only data shared by all threads
struct sieve_ctx
{
	const poly* f;
	const bigint* m;
	const factor_bases* fb;
	side_tab rt, at; // rational and algebraic entries
	struct proj_t
	{
		uint64_t p;
		float lg;
		uint32_t idx;
		unsigned k; // p^k || lc(f)
	};
	vector<proj_t> projs; // projective primes (p | lc(f))
	int64_t A;
	uint32_t len;
	uint64_t lpr, lpa; // large prime bounds per side (0 = none)
	int rthr, athr; // candidate thresholds in sieve units
	double logm;
	size_t d;
	vector<double> fc; // f as doubles
	size_t alg_off, quad_off;
	bool exact_norm;
};

// per phase seconds and counters, summed over the threads
struct phase_stats
{
	double fill = 0, sieve = 0, scan = 0, hits = 0, fact = 0;
	uint64_t ncand = 0, nlines = 0;
	void add(const phase_stats& o)
	{
		fill += o.fill; sieve += o.sieve; scan += o.scan; hits += o.hits; fact += o.fact;
		ncand += o.ncand; nlines += o.nlines;
	}
};

// fills out[0..L) with SC*(log|N(a+b*theta)| - proj_lg) for the indices base..base+L
// input: c[k] = f_k * (-b)^(d-k) (so the norm form is sum c[k]*(-a)^k), d = degree, A = half width
// the norm is evaluated at every NORM_STRIDE-th index and interpolated, except where it changes sign or jumps
static void fill_alg_norm(int16_t* out, uint32_t base, uint32_t L, int64_t A, const double* c, size_t d, double proj_lg, bool exact, vector<float>& lv, vector<char>& sg)
{
	auto eval = [&](uint32_t j) { // Horner in -a
		const double ap = -(double)((int64_t)(base + j) - A);
		double h = c[d];
		for (size_t k = d; k-- > 0;) h = h * ap + c[k];
		return h;
	};
	auto lnv = [&](double h) { return (float)(SC * (fast_log(std::max(fabs(h), 1.0)) - proj_lg)); }; // same estimate as before, in sieve units

	if (exact)
	{
		for (uint32_t j = 0; j < L; ++j) out[j] = rnd16(lnv(eval(j)));
		return;
	}

	const uint32_t ns = (L + NORM_STRIDE - 1) / NORM_STRIDE + 1; // sample points j = 0, S, 2S, ... (the last one may lie past L)
	lv.resize(ns);
	sg.resize(ns);
	for (uint32_t k = 0; k < ns; ++k)
	{
		const double h = eval(k * NORM_STRIDE);
		lv[k] = lnv(h);
		sg[k] = h < 0;
	}
	for (uint32_t k = 0; k + 1 < ns; ++k)
	{
		const uint32_t j0 = k * NORM_STRIDE, j1 = std::min<uint32_t>(L, j0 + NORM_STRIDE);
		if (j0 >= L) break;
		if (sg[k] == sg[k + 1] && fabsf(lv[k + 1] - lv[k]) <= NORM_JUMP)
		{
			const float dv = (lv[k + 1] - lv[k]) * (1.0f / NORM_STRIDE);
			for (uint32_t j = j0; j < j1; ++j) out[j] = rnd16(lv[k] + dv * (float)(j - j0));
		}
		else
			for (uint32_t j = j0; j < j1; ++j) out[j] = rnd16(lnv(eval(j))); // near a root of the norm form: exact
	}
}

// the sieve of one thread: owns every array that is written to
struct line_sieve
{
	const sieve_ctx& X;
	side_work R, Aw; // rational and algebraic
	vector<uint64_t> bm; // marks candidate offsets inside the segment
	vector<uint32_t> slot; // marked offset -> candidate number
	vector<uint32_t> js; // candidate offsets of the segment
	vector<float> lv;
	vector<char> sg;
	vector<double> c, bp;
	phase_stats st;

	explicit line_sieve(const sieve_ctx& x) : X(x), bm(SEG / 64, 0), slot(SEG, 0), c(x.d + 1), bp(x.d + 1)
	{
		R.init(X.rt);
		Aw.init(X.at);
	}

	// sieves b = b0 ... b1, appends full relations to out and relations with a large prime to partials
	void run_chunk(uint64_t b0, uint64_t b1, vector<nfs_relation>& out, vector<nfs_relation>& partials)
	{
		R.init_offsets(b0, X.A);
		Aw.init_offsets(b0, X.A);
		for (uint64_t b = b0; b <= b1; ++b)
		{
			if (b > b0)
			{
				R.advance_offsets();
				Aw.advance_offsets();
			}
			do_line(b, out, partials);
			++st.nlines;
		}
	}

	// one line: all segments, relations appended to out / partials
	void do_line(uint64_t b, vector<nfs_relation>& out, vector<nfs_relation>& partials)
	{
		using clk = chrono::steady_clock;
		auto sec = [](clk::time_point a, clk::time_point z) { return chrono::duration<double>(z - a).count(); };

		R.start_line(b);
		Aw.start_line(b);

		float proj_lg = 0; // projective primes dividing b divide every norm of the line
		vector<uint32_t> proj_idx;
		for (const auto& e : X.projs)
			if (b % e.p == 0)
			{
				unsigned vb = 0; // v_p(b)
				for (uint64_t t = b; t % e.p == 0; t /= e.p) ++vb;
				proj_lg += std::min(vb, e.k) * e.lg; // v_p(norm) is about min(v_p(b), k), can be larger when equal (THRESH_MARGIN covers a little)
				proj_idx.push_back(e.idx);
			}

		bp[0] = 1;
		for (size_t j = 1; j <= X.d; ++j) bp[j] = bp[j - 1] * -(double)b;
		for (size_t k = 0; k <= X.d; ++k) c[k] = X.fc[k] * bp[X.d - k];
		const int16_t R0 = rnd16((float)(SC * (log((double)b) + X.logm))); // rational side: log(b*m) for the whole line

		for (uint32_t sb = 0; sb < X.len; sb += SEGS_PER_SUPER * SEG)
		{
			const uint32_t se = (uint32_t)std::min<uint64_t>(X.len, (uint64_t)sb + (uint64_t)SEGS_PER_SUPER * SEG);
			const auto t0 = clk::now();
			R.fill_buckets(sb, se);
			Aw.fill_buckets(sb, se);
			st.sieve += sec(t0, clk::now());

			for (uint32_t base = sb; base < se; base += SEG)
			{
				const uint32_t L = std::min<uint32_t>(SEG, se - base);
				const uint32_t s = (base - sb) >> SEG_LOG;

				auto t1 = clk::now();
				std::fill(R.arr.begin(), R.arr.begin() + L, R0);
				fill_alg_norm(Aw.arr.data(), base, L, X.A, c.data(), X.d, proj_lg, X.exact_norm, lv, sg);
				auto t2 = clk::now();
				st.fill += sec(t1, t2);

				R.spos0 = R.spos;
				Aw.spos0 = Aw.spos; // start of the segment, for the divisor search
				R.sieve_small(base, L);
				Aw.sieve_small(base, L);
				R.apply_bucket(s);
				Aw.apply_bucket(s);
				auto t3 = clk::now();
				st.sieve += sec(t2, t3);

				// candidates: both arrays near 0 (up to the large prime allowance), gcd(a,b)=1
				js.clear();
				const int16_t* sr = R.arr.data();
				const int16_t* sa = Aw.arr.data();
				auto take = [&](uint32_t j) { // position j passed the thresholds: keep it if gcd(a,b)=1
					const int64_t a = (int64_t)(base + j) - X.A;
					if (gcd((uint64_t)(a < 0 ? -a : a), b) == 1) js.push_back(j);
				};
				uint32_t j = 0;
				for (; j + 8 <= L; j += 8) // 8 positions at a time, one rarely taken branch (no per position branch to mispredict)
				{
					int any = 0;
					for (uint32_t i = 0; i < 8; ++i) any |= (sr[j + i] <= X.rthr) & (sa[j + i] <= X.athr);
					if (!any) continue;
					for (uint32_t i = 0; i < 8; ++i)
						if (sr[j + i] <= X.rthr && sa[j + i] <= X.athr) take(j + i);
				}
				for (; j < L; ++j)
					if (sr[j] <= X.rthr && sa[j] <= X.athr) take(j);
				auto t4 = clk::now();
				st.scan += sec(t3, t4);
				if (js.empty()) continue;
				st.ncand += js.size();

				const size_t nc = js.size();
				for (auto* h : { &R.hits, &Aw.hits })
				{
					if (h->size() < nc) h->resize(nc);
					for (size_t k = 0; k < nc; ++k) (*h)[k].clear();
				}
				for (size_t k = 0; k < nc; ++k)
				{
					bm[js[k] >> 6] |= uint64_t(1) << (js[k] & 63);
					slot[js[k]] = (uint32_t)k;
				}
				R.collect_hits(base, L, s, js, bm, slot);
				Aw.collect_hits(base, L, s, js, bm, slot);
				if (!proj_idx.empty())
					for (size_t k = 0; k < nc; ++k)
					{
						for (uint32_t i : proj_idx) Aw.hits[k].push_back(i);
						sort(Aw.hits[k].begin(), Aw.hits[k].end()); // factor_*_hits wants ascending indices
					}
				for (uint32_t j : js) bm[j >> 6] = 0;
				auto t5 = clk::now();
				st.hits += sec(t4, t5);

				for (size_t k = 0; k < nc; ++k)
				{
					const int64_t a = (int64_t)(base + js[k]) - X.A;
					nfs_relation rel;
					rel.a = a;
					rel.b = (int64_t)b;
					uint64_t lr = 0, la = 0;
					if (!factor_rational_hits(a, rel.b, *X.m, *X.fb, R.hits[k], X.lpr, lr, rel.cols, rel.rfac)) continue;
					if (!factor_algebraic_hits(*X.f, a, rel.b, *X.fb, X.alg_off, Aw.hits[k], X.lpa, la, rel.cols, rel.afac)) continue;
					if (!quadratic_chars(a, rel.b, *X.fb, X.quad_off, rel.cols)) continue;
					if (X.fb->parity) rel.cols.push_back(uint32_t(X.quad_off + X.fb->quadratic.size()));

					rel.lp_r = lr;
					rel.lp_a = la;
					if (la && !alg_root(a, rel.b, la, rel.lp_a_root)) continue; // p | b: not a usable ideal

					if (lr || la) partials.push_back(move(rel));
					else out.push_back(move(rel));
				}
				st.fact += sec(t5, clk::now());
			}
		}
	}
};

} // namespace

// appends relations to out until out.size() + (cycles of the large prime graph) >= sp.target or b exceeds b_max
// relations with a large prime go to partials; next_b is the first b still to be processed (in/out)
// threads: GNFS_THREADS (default: all cores) take chunks of consecutive b; every b below the returned next_b is complete
// GNFS_EXACT_NORM=1 evaluates the algebraic norm at every position (for checking the interpolation)
void find_relations(const poly& f, const bigint& m, const factor_bases& fb, sieve_params sp, vector<nfs_relation>& out, vector<nfs_relation>& partials, uint64_t& next_b)
{
	if (fb.rational.empty()) throw domain_error("find_relations: empty rational base");

	const size_t ncol = num_columns(fb);
	const size_t d = f.size() - 1;
	if (sp.target == 0) sp.target = ncol + 32;
	if (sp.A <= 0) sp.A = std::max<int64_t>(10000, 4 * (int64_t)fb.rational.back().p);
	const uint64_t b_max = sp.b_max ? sp.b_max : DEFAULT_B_MAX;

	const int64_t A = sp.A;
	const uint64_t len64 = 2 * (uint64_t)A + 1;
	if (len64 >= (1ull << 31)) throw domain_error("find_relations: sieve line too long (positions are 32 bit)");

	// large prime bounds per side, clamped to B^2 so that any cofactor <= bound is automatically prime
	uint64_t Br = fb.rational.back().p, Ba = 0;
	for (const fb_pair& e : fb.algebraic)
		if (e.r != FB_PROJECTIVE) Ba = std::max<uint64_t>(Ba, e.p);
	if (Br >= (1ull << 31) || Ba >= (1ull << 31)) throw domain_error("find_relations: factor base primes must be below 2^31");
	uint64_t lpr = 0, lpa = 0;
	if (sp.lp_bound)
	{
		lpr = std::min<uint64_t>(sp.lp_bound, Br * Br);
		lpa = std::min<uint64_t>(sp.lp_bound, Ba * Ba);
		if (lpr <= Br) lpr = 0; // nothing above the base: no large primes on this side
		if (lpa <= Ba) lpa = 0;
	}
	const float rthr = sp.slack + (lpr ? (float)log((double)lpr) : 0.f); // candidate thresholds on the residual logs
	const float athr = sp.aslack + (lpa ? (float)log((double)lpa) : 0.f);

	sieve_ctx X;
	X.f = &f;
	X.m = &m;
	X.fb = &fb;
	X.A = A;
	X.len = (uint32_t)len64;
	X.lpr = lpr;
	X.lpa = lpa;
	X.rthr = (int)std::floor(SC * rthr) + THRESH_MARGIN;
	X.athr = (int)std::floor(SC * athr) + THRESH_MARGIN;
	X.logm = log_cpp(m);
	X.d = d;
	X.fc.resize(d + 1);
	for (size_t k = 0; k <= d; ++k) X.fc[k] = f[k].convert_to<double>(); // f as doubles, for the norm estimate
	X.alg_off = 1 + fb.rational.size(); // first column of each block
	X.quad_off = X.alg_off + fb.algebraic.size();
	{
		const char* e = getenv("GNFS_EXACT_NORM");
		X.exact_norm = e && atoi(e) != 0;
	}

	// entry tables: rational primes and the prime powers that fit in the interval, algebraic degree-1 primes
	{
		vector<ent> rv, av;
		for (size_t i = 0; i < fb.rational.size(); ++i)
		{
			const fb_pair& e = fb.rational[i];
			const int16_t lg = rnd16((float)(SC * log((double)e.p)));
			rv.push_back({ (uint32_t)e.p, (uint32_t)e.p, (uint32_t)e.r, uint32_t(i), lg });
			for (uint64_t pk = e.p * e.p; pk <= len64; pk *= e.p)
				rv.push_back({ (uint32_t)pk, (uint32_t)e.p, (uint32_t)(m % bigint(pk)).convert_to<uint64_t>(), NO_IDX, lg });
		}
		for (size_t i = 0; i < fb.algebraic.size(); ++i)
		{
			const fb_pair& e = fb.algebraic[i];
			if (e.r == FB_PROJECTIVE)
			{
				unsigned k = 0; // multiplicity of p in lc(f)
				for (bigint t = f.data().back(); t % bigint(e.p) == 0; t /= bigint(e.p)) ++k;
				X.projs.push_back({ e.p, (float)log((double)e.p), uint32_t(i), k }); // p | b makes p^min(v_p(b), k) divide the whole line
			}
			else av.push_back({ (uint32_t)e.p, (uint32_t)e.p, (uint32_t)e.r, uint32_t(i), rnd16((float)(SC * log((double)e.p))) });
		}
		X.rt = make_side_tab(rv);
		X.at = make_side_tab(av);
	}

	cycle_counter cc; // exact cycle count, updated with every partial relation
	for (const nfs_relation& x : partials) cc.add(x);

	unsigned nthreads = thread::hardware_concurrency();
	if (const char* e = getenv("GNFS_THREADS")) nthreads = (unsigned)std::max(1, atoi(e));
	if (nthreads == 0) nthreads = 1;

	mutex mu; // protects out, partials, cc, stats, progress
	atomic<uint64_t> claim(next_b); // next unclaimed b
	atomic<bool> stop(out.size() + cc.cycles >= sp.target);
	phase_stats total;

	auto worker = [&]() {
		line_sieve ls(X);
		vector<nfs_relation> lout, lpart;
		while (!stop.load())
		{
			const uint64_t c0 = claim.fetch_add(LINES_PER_CHUNK);
			if (c0 > b_max) break;
			const uint64_t c1 = std::min<uint64_t>(c0 + LINES_PER_CHUNK - 1, b_max);
			lout.clear();
			lpart.clear();
			ls.run_chunk(c0, c1, lout, lpart); // a claimed chunk is always finished, so next_b never skips a line

			lock_guard<mutex> lk(mu);
			for (nfs_relation& r : lpart)
			{
				cc.add(r); // before the move
				partials.push_back(move(r));
			}
			for (nfs_relation& r : lout) out.push_back(move(r));
			total.add(ls.st);
			ls.st = phase_stats();
			if (out.size() + cc.cycles >= sp.target) stop = true;
			const uint64_t c1floor = 500 * (c1 / 500);
			if (c0 < c1floor)
			{
				cerr << "b=" << c1 << "  full " << out.size() << ", partial " << partials.size() << " (" << cc.cycles << " cycles), target " << sp.target << ", "
				     << 100. * (out.size() + cc.cycles) / sp.target << " percent done; " << (double)total.ncand / std::max<uint64_t>(1, total.nlines) << " candidates/line; cpu s/line: fill "
				     << total.fill / std::max<uint64_t>(1, total.nlines) << " sieve " << total.sieve / std::max<uint64_t>(1, total.nlines) << " scan " << total.scan / std::max<uint64_t>(1, total.nlines)
				     << " hits " << total.hits / std::max<uint64_t>(1, total.nlines) << " factor " << total.fact / std::max<uint64_t>(1, total.nlines) << '\n';
			}
		}
	};

	if (!stop.load())
	{
		vector<thread> pool;
		for (unsigned t = 1; t < nthreads; ++t) pool.emplace_back(worker);
		worker();
		for (thread& t : pool) t.join();
	}

	next_b = std::max<uint64_t>(next_b, std::min<uint64_t>(claim.load(), b_max + 1)); // never moves back (next_b may already exceed b_max)
	if (total.nlines)
		cerr << "sieve phases, cpu s per line (" << total.nlines << " lines, " << nthreads << " threads): fill " << total.fill / total.nlines << " sieve " << total.sieve / total.nlines
		     << " scan " << total.scan / total.nlines << " hits " << total.hits / total.nlines << " factor " << total.fact / total.nlines << "; " << (double)total.ncand / total.nlines << " candidates/line\n";
	cerr << "final b=" << next_b << "\n";
}
