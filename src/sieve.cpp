// relation finding: line sieve on both sides, candidate factoring (trial or re-sieve), quadratic characters,
// single large prime per side (partial relations combined through cycles in the large-prime graph)

#include "sieve.hpp"
#include "numtheory.hpp"
#include <map>
#include <tuple>
#include <algorithm>
#include <cstring>

// ---------------------------------------------------------------------------------------------
// cofactor handling
// after dividing out every factor-base prime, v is either 1, or a product of primes > B.
// if v <= lp_bound <= B^2 it must be a single prime: accept it as a large prime.
static bool take_cofactor(const bigint& v, uint64_t lp_bound, uint64_t& lp)
{
	lp = 0;
	if (v == 1) return true;
	if (lp_bound && v <= bigint(lp_bound)) { lp = v.convert_to<uint64_t>(); return true; }
	return false;
}

// a+b*m smooth over rational fb (up to one large prime)?
// in-place update cols with sign(a+b*m) and prime factors
static bool factor_rational(int64_t a, int64_t b, const bigint& m, const factor_bases& fb, uint64_t lp_bound, uint64_t& lp, vector<uint32_t>& cols, vector<fac_entry>& fac)
{
	bigint v = bigint(a) + bigint(b) * m;
	if (v < 0) { cols.push_back(0); v = -v; }

	for (size_t i = 0; i < fb.rational.size(); ++i)
	{
		const uint64_t p = fb.rational[i].p;
		if (lin_mod(a, b, p, fb.rational[i].r) != 0) continue; // p does not divide a+bm

		const bigint P = p;
		unsigned e = 0;
		while (v % P == 0) { v /= P; ++e; }
		if (e) fac.push_back({ uint32_t(i), e });
		if (e & 1) cols.push_back(uint32_t(1 + i));
	}

	return take_cofactor(v, lp_bound, lp);
}

// N(a+b*theta) smooth over algebraic fb (up to one large prime)?
// a degree-1 prime ideal (p,r)|a+b*theta iff a+b*r=0 mod p
// projective ideal | a+b*theta iff p|b
static bool factor_algebraic(const poly& f, int64_t a, int64_t b, const factor_bases& fb, size_t offset, uint64_t lp_bound, uint64_t& lp, vector<uint32_t>& cols, vector<fac_entry>& fac)
{
	bigint N = abs(norm(f, a, b));
	if (N == 0) return false;

	for (size_t i = 0; i < fb.algebraic.size(); ++i)
	{
		const uint64_t p = fb.algebraic[i].p, r = fb.algebraic[i].r;
		const bool hit = (r == FB_PROJECTIVE) ? (uint64_t(b) % p == 0) : (lin_mod(a, b, p, r) == 0);
		if (!hit) continue;

		const bigint P = p;
		unsigned e = 0;
		while (N % P == 0) { N /= P; ++e; }
		if (e) fac.push_back({ uint32_t(i), e });
		if (e & 1) cols.push_back(uint32_t(offset + i));
	}

	return take_cofactor(N, lp_bound, lp);
}

// same, but the primes dividing the value are already known (hits = ascending indices into the base)
static bool factor_rational_hits(int64_t a, int64_t b, const bigint& m, const factor_bases& fb, const vector<uint32_t>& hits, uint64_t lp_bound, uint64_t& lp, vector<uint32_t>& cols, vector<fac_entry>& fac)
{
	bigint v = bigint(a) + bigint(b) * m;
	if (v < 0) { cols.push_back(0); v = -v; }

	for (uint32_t i : hits)
	{
		const bigint P = fb.rational[i].p;
		unsigned e = 0;
		while (v % P == 0) { v /= P; ++e; }
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
		while (N % P == 0) { N /= P; ++e; }
		if (e) fac.push_back({ i, e });
		if (e & 1) cols.push_back(uint32_t(offset + i));
	}

	return take_cofactor(N, lp_bound, lp);
}

// quadratic characters
// in-place update cols if ((a+b*s)/q)=-1
// returns false if a+b*s = 0 mod q
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

// ---------------------------------------------------------------------------------------------
// large-prime helpers

static uint64_t mulmod_u64(uint64_t a, uint64_t b, uint64_t m) { return (uint64_t)((unsigned __int128)a * b % m); }

// inverse of a mod prime p (a != 0 mod p), extended Euclid
static uint64_t inv_mod_u64(uint64_t a, uint64_t p)
{
	__int128 r0 = p, r1 = a % p, t0 = 0, t1 = 1;
	while (r1 != 0)
	{
		const __int128 q = r0 / r1;
		const __int128 r2 = r0 - q * r1; r0 = r1; r1 = r2;
		const __int128 t2 = t0 - q * t1; t0 = t1; t1 = t2;
	}
	if (t0 < 0) t0 += p;
	return (uint64_t)t0;
}

// root of the algebraic ideal of degree 1 above p dividing a+b*theta: r = -a/b mod p.  false if p | b.
static bool alg_root(int64_t a, int64_t b, uint64_t p, uint64_t& r)
{
	const uint64_t bm = uint64_t(b) % p;
	if (bm == 0) return false;
	int64_t am = a % (int64_t)p; if (am < 0) am += (int64_t)p;
	r = (p - mulmod_u64((uint64_t)am, inv_mod_u64(bm, p), p)) % p;
	return true;
}

// graph: vertex 0 = "no large prime"; every other vertex is a rational prime p or an algebraic ideal (p,r).
// a partial relation is an edge between its two vertices.
struct lp_graph
{
	std::map<std::tuple<int, uint64_t, uint64_t>, uint32_t> id;
	vector<std::array<uint32_t, 2>> edges;

	explicit lp_graph(const vector<nfs_relation>& rels)
	{
		auto vid = [&](int side, uint64_t p, uint64_t r) -> uint32_t {
			if (p == 0) return 0;
			auto it = id.find({ side, p, r });
			if (it != id.end()) return it->second;
			const uint32_t v = uint32_t(id.size() + 1);
			id.emplace(std::make_tuple(side, p, r), v);
			return v;
		};
		for (const nfs_relation& x : rels)
			edges.push_back({ vid(0, x.lp_r, 0), vid(1, x.lp_a, x.lp_a_root) });
	}
	size_t nvert() const { return id.size() + 1; }
};

static uint32_t uf_find(vector<uint32_t>& p, uint32_t x) { while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; } return x; }

// number of independent cycles = number of full relations obtainable by combining the partials
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
		if (u == v) ++cycles; else uf[u] = v;
	}
	return cycles;
}

// append one combined relation per fundamental cycle of the large-prime graph to out
void combine_partials(const vector<nfs_relation>& partials, vector<nfs_relation>& out)
{
	if (partials.empty()) return;
	lp_graph g(partials);
	const size_t V = g.nvert(), E = g.edges.size();

	vector<vector<uint32_t>> adj(V); // edge indices
	for (size_t e = 0; e < E; ++e) { adj[g.edges[e][0]].push_back((uint32_t)e); adj[g.edges[e][1]].push_back((uint32_t)e); }

	const uint32_t NONE = UINT32_MAX;
	vector<uint32_t> pv(V, NONE), pe(V, NONE), depth(V, 0);
	vector<char> seen(V, 0), tree(E, 0);

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
				seen[v] = 1; pv[v] = u; pe[v] = e; depth[v] = depth[u] + 1; tree[e] = 1;
				q.push_back(v);
			}
		}
	}

	for (size_t e0 = 0; e0 < E; ++e0)
	{
		if (tree[e0]) continue;
		vector<uint32_t> used{ (uint32_t)e0 };
		uint32_t u = g.edges[e0][0], v = g.edges[e0][1];
		while (u != v)
		{
			if (depth[u] >= depth[v]) { used.push_back(pe[u]); u = pv[u]; }
			else                      { used.push_back(pe[v]); v = pv[v]; }
		}

		nfs_relation c;
		c.a = partials[used[0]].a;
		c.b = partials[used[0]].b;
		vector<uint32_t> cols;
		for (uint32_t e : used)
			cols.insert(cols.end(), partials[e].cols.begin(), partials[e].cols.end());
		std::sort(cols.begin(), cols.end());
		for (size_t i = 0; i < cols.size();) // keep columns with odd multiplicity
		{
			size_t j = i;
			while (j < cols.size() && cols[j] == cols[i]) ++j;
			if ((j - i) & 1) c.cols.push_back(cols[i]);
			i = j;
		}
		c.members = std::move(used);
		out.push_back(std::move(c));
	}
}

// ---------------------------------------------------------------------------------------------
// line sieve over a in [-A,A] for b = 1,2,3...
// two arrays per b:
//  rational:  starts at log|a+bm|~log(b*m); a prime p | a+b*m subtracts log(p)
//  algebraic: starts at log|N(a+b*theta)| (evaluated in double); (p,r) with a+b*r=0 (mod p) subtracts log(p)
// candidates: residual log below slack (+ log of the large-prime bound, if enabled) on BOTH sides.
// those candidates are checked precisely afterwards:
// rational factorization, algebraic norm factorization, quadratic characters.
// a cofactor that is a single prime <= lp_bound on one or both sides gives a partial relation (-> partials);
// relations with no large prime go straight to out.  Call combine_partials(partials, out) when done.

// approximate ln(x) for x >= 1, no libm, no division, no branch:  x = 2^e * (1+t), t in [0,1),
// log2(1+t) ~ t (1.3465 - 0.3465 t)  (measured max error 0.005 in ln).
static inline double fast_ln(double x)
{
	uint64_t bits;
	std::memcpy(&bits, &x, sizeof bits);
	const int e = int((bits >> 52) & 0x7ff) - 1023;
	bits = (bits & 0x000fffffffffffffULL) | 0x3ff0000000000000ULL;
	double m;
	std::memcpy(&m, &bits, sizeof m);
	const double t = m - 1.0;
	return ((double)e + t * (1.3465 - 0.3465 * t)) * 0.6931471805599453;
}

void find_relations(const poly& f, const bigint& m, const factor_bases& fb, sieve_params sp, vector<nfs_relation>& out, vector<nfs_relation>& partials, uint64_t& next_b)
{
	if (fb.rational.empty()) throw domain_error("find_relations: empty rational base");

	const size_t ncol = num_columns(fb);
	const size_t d = f.size() - 1;
	if (sp.target == 0) sp.target = ncol + 32;
	if (sp.A <= 0) sp.A = std::max<int64_t>(10000, 4 * (int64_t)fb.rational.back().p);
	const uint64_t b_max = sp.b_max ? sp.b_max : DEFAULT_B_MAX;

	const int64_t A = sp.A;
	const uint64_t len = 2 * (uint64_t)A + 1;

	// large prime bounds per side, clamped to B^2 so that any cofactor <= bound is automatically prime
	uint64_t Br = fb.rational.back().p, Ba = 0;
	for (const fb_pair& e : fb.algebraic) if (e.r != FB_PROJECTIVE) Ba = std::max<uint64_t>(Ba, e.p);
	uint64_t lpr = 0, lpa = 0;
	if (sp.lp_bound)
	{
		lpr = std::min<uint64_t>(sp.lp_bound, Br * Br);
		lpa = std::min<uint64_t>(sp.lp_bound, Ba * Ba);
		if (lpr <= Br) lpr = 0;
		if (lpa <= Ba) lpa = 0;
	}
	const float rthr = sp.slack  + (lpr ? (float)std::log((double)lpr) : 0.f);
	const float athr = sp.aslack + (lpa ? (float)std::log((double)lpa) : 0.f);

	struct entry { uint64_t pk, p, r; float lg; uint32_t idx; };
	vector<entry> ents; // rational primes and prime powers that fit in the interval
	for (size_t i = 0; i < fb.rational.size(); ++i)
	{
		const fb_pair& e = fb.rational[i];
		const float lg = (float)std::log((double)e.p);
		ents.push_back({ e.p, e.p, e.r, lg, uint32_t(i) });
		for (uint64_t pk = e.p * e.p; pk <= len; pk *= e.p)
			ents.push_back({ pk, e.p, (m % bigint(pk)).convert_to<uint64_t>(), lg, uint32_t(i) });
	}

	struct aentry { uint64_t p, r; float lg; uint32_t idx; };
	vector<aentry> aents; // algebraic degree-1 primes
	struct pentry { uint64_t p; float lg; uint32_t idx; };
	vector<pentry> projs; // p | lc(f): the root at infinity.  p | F(a,b) for ALL a in the lines with p | b
	for (size_t i = 0; i < fb.algebraic.size(); ++i)
	{
		const fb_pair& e = fb.algebraic[i];
		if (e.r == FB_PROJECTIVE) projs.push_back({ e.p, (float)std::log((double)e.p), uint32_t(i) });
		else aents.push_back({ e.p, e.r, (float)std::log((double)e.p), uint32_t(i) });
	}

	// cost model (very rough, in "memory steps")
	double resieve_cost = 0;
	for (const entry& e : ents) if (e.pk == e.p) resieve_cost += (double)len / (double)e.p + 2;
	for (const aentry& e : aents) resieve_cost += (double)len / (double)e.p + 2;
	const double test_cost_per_cand = 4.0 * (double)(fb.rational.size() + fb.algebraic.size());

	vector<double> fc(d + 1);
	for (size_t k = 0; k <= d; ++k) fc[k] = f[k].convert_to<double>();

	const size_t alg_off  = 1 + fb.rational.size();
	const size_t quad_off = alg_off + fb.algebraic.size();
	const double logm = log_cpp(m);

	vector<float> s(len), sa(len);
	vector<double> bp(d + 1);
	size_t ncand = 0;
	vector<uint32_t> cands;
	vector<uint64_t> bitmap;
	vector<uint32_t> slot;
	vector<vector<uint32_t>> hit_r, hit_a;

	size_t cycles = count_lp_cycles(partials);

	uint64_t b = next_b;
	for (; b <= b_max && out.size() + cycles < sp.target; ++b)
	{
		std::fill(s.begin(), s.end(), (float)(std::log((double)b) + logm));
		float proj_lg = 0;
		for (const pentry& e : projs) if (b % e.p == 0) proj_lg += e.lg;

		// log|N(a+b*theta)| = log|sum f_k (-a)^k (-b)^(d-k)|
		bp[0] = 1;
		for (size_t j = 1; j <= d; ++j) bp[j] = bp[j - 1] * -(double)b;
		for (uint64_t i = 0; i < len; ++i)
		{
			const double ap = -(double)((int64_t)i - A);
			double apow = 1.0;
			double h = 0.0;
			for (size_t k = 0; k <= d; ++k)
			{
				h += fc[k] * apow * bp[d - k];
				apow *= ap;
			}
			sa[i] = (float)fast_ln(std::max(std::fabs(h), 1.0)) - proj_lg;
		}

		for (const entry& e : ents)
		{
			if (b % e.p == 0) continue; // a+bm = a mod p, and gcd(a,b)=1 excludes it
			const uint64_t a0 = (e.pk - (b % e.pk) * e.r % e.pk) % e.pk; // a = -b*m mod pk
			for (uint64_t i = (a0 + (uint64_t)A) % e.pk; i < len; i += e.pk) s[i] -= e.lg;
		}

		for (const aentry& e : aents)
		{
			if (b % e.p == 0) continue;
			const uint64_t a0 = (e.p - (b % e.p) * e.r % e.p) % e.p; // a = -b*r mod p
			for (uint64_t i = (a0 + (uint64_t)A) % e.p; i < len; i += e.p) sa[i] -= e.lg;
		}

		// candidates: both arrays near 0 (up to the large prime allowance), gcd(a,b)=1
		cands.clear();
		for (uint64_t i = 0; i < len; ++i)
		{
			if (s[i] > rthr || sa[i] > athr) continue;
			const int64_t a = (int64_t)i - A;
			if (std::gcd((uint64_t)(a < 0 ? -a : a), b) != 1) continue;
			cands.push_back((uint32_t)i);
		}
		ncand += cands.size();

		const bool resieve = !cands.empty() &&
			(sp.cand_mode == 2 || (sp.cand_mode == 0 && (double)cands.size() * test_cost_per_cand > resieve_cost));

		if (resieve)
		{
			if (bitmap.size() != (len + 63) / 64) { bitmap.assign((len + 63) / 64, 0); slot.assign(len, 0); }
			if (hit_r.size() < cands.size()) { hit_r.resize(cands.size()); hit_a.resize(cands.size()); }
			for (size_t k = 0; k < cands.size(); ++k)
			{
				bitmap[cands[k] >> 6] |= uint64_t(1) << (cands[k] & 63);
				slot[cands[k]] = (uint32_t)k;
				hit_r[k].clear();
				hit_a[k].clear();
			}

			for (const entry& e : ents)
			{
				if (e.pk != e.p || b % e.p == 0) continue; // prime powers divide only where the prime does
				const uint64_t a0 = (e.p - (b % e.p) * e.r % e.p) % e.p;
				for (uint64_t i = (a0 + (uint64_t)A) % e.p; i < len; i += e.p)
					if (bitmap[i >> 6] >> (i & 63) & 1) hit_r[slot[i]].push_back(e.idx);
			}
			for (const aentry& e : aents)
			{
				if (b % e.p == 0) continue;
				const uint64_t a0 = (e.p - (b % e.p) * e.r % e.p) % e.p;
				for (uint64_t i = (a0 + (uint64_t)A) % e.p; i < len; i += e.p)
					if (bitmap[i >> 6] >> (i & 63) & 1) hit_a[slot[i]].push_back(e.idx);
			}

			if (proj_lg > 0) // the projective primes dividing b divide every candidate of the line
				for (size_t k = 0; k < cands.size(); ++k)
				{
					for (const pentry& e : projs) if (b % e.p == 0) hit_a[k].push_back(e.idx);
					std::sort(hit_a[k].begin(), hit_a[k].end());
				}

			for (uint32_t pos : cands) bitmap[pos >> 6] = 0;
		}

		for (size_t k = 0; k < cands.size(); ++k)
		{
			const int64_t a = (int64_t)cands[k] - A;

			nfs_relation rel;
			rel.a = a;
			rel.b = (int64_t)b;
			uint64_t lr = 0, la = 0;
			if (resieve)
			{
				if (!factor_rational_hits(a, rel.b, m, fb, hit_r[k], lpr, lr, rel.cols, rel.rfac)) continue;
				if (!factor_algebraic_hits(f, a, rel.b, fb, alg_off, hit_a[k], lpa, la, rel.cols, rel.afac)) continue;
			}
			else
			{
				if (!factor_rational(a, rel.b, m, fb, lpr, lr, rel.cols, rel.rfac)) continue;
				if (!factor_algebraic(f, a, rel.b, fb, alg_off, lpa, la, rel.cols, rel.afac)) continue;
			}
			if (!quadratic_chars(a, rel.b, fb, quad_off, rel.cols)) continue;
			if (fb.parity) rel.cols.push_back(uint32_t(quad_off + fb.quadratic.size()));

			rel.lp_r = lr;
			rel.lp_a = la;
			if (la && !alg_root(a, rel.b, la, rel.lp_a_root)) continue; // p | b: not a usable ideal

			if (lr || la) partials.push_back(std::move(rel));
			else          out.push_back(std::move(rel));
		}

		if (b % 100 == 0 && !partials.empty()) cycles = count_lp_cycles(partials);

		if (b % 500 == 0)
			cerr << "b=" << b << "  full " << out.size() << ", partial " << partials.size()
			     << " (" << cycles << " cycles), target " << sp.target << " (" << ncand << " candidates checked)\n";
	}

	next_b = b;
}
