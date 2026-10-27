// product and remainder trees + CRT reconstruction

#include "crt_tree.hpp"
#include "numtheory.hpp"

struct prod_tree
{
	vector<vector<bigint>> level; // level[0]=leaves, level.back()={product of all}

	explicit prod_tree(vector<bigint> leaves)
	{
		if (leaves.empty()) leaves.push_back(1);
		level.push_back(std::move(leaves));
		while (level.back().size() > 1)
		{
			const auto& cur = level.back();
			vector<bigint> nxt;
			for (size_t i = 0; i + 1 < cur.size(); i += 2) nxt.push_back(cur[i] * cur[i + 1]); // node value = product of children
			if (cur.size() & 1) nxt.push_back(cur.back());
			level.push_back(std::move(nxt));
		}
	}

	const bigint& root() const { return level.back()[0]; }
};

// return product of entries in v
bigint balanced_product(vector<bigint> v) { return prod_tree(std::move(v)).root(); }

// returns value mod (leaf_i)^2 for every leaf
static vector<bigint> rem_tree_sq(const bigint& value, const prod_tree& T)
{
	vector<bigint> cur{ value % (T.root() * T.root()) };
	for (size_t lv = T.level.size() - 1; lv-- > 0;)
	{
		const auto& nodes = T.level[lv];
		vector<bigint> nxt(nodes.size());
		// split each value into two branches
		for (size_t i = 0; i < nodes.size(); ++i) nxt[i] = cur[i / 2] % (nodes[i] * nodes[i]);
		cur = std::move(nxt);
	}
	return cur;
}

// value mod p for each prime
vector<u64> reduce_all(const bigint& value, const vector<u64>& primes)
{
	vector<u64> out;
	if (primes.empty()) return out;
	// construct product tree to minimize # of large mods computed
	const prod_tree T(vector<bigint>(primes.begin(), primes.end()));
	const vector<bigint> r = rem_tree_sq(value, T);
	for (size_t i = 0; i < primes.size(); ++i) out.push_back(mod_u64(r[i], primes[i]));
	return out;
}

// returns vector X, |X[c]|<P/2, s.t. X[c] = res[i][c] (mod primes[i])
// for all i and c < d
vector<bigint> crt_tree(const vector<u64>& primes, const vector<vector<u64>>& res, size_t d)
{
	const size_t k = primes.size();
	const prod_tree T(vector<bigint>(primes.begin(), primes.end()));
	const bigint& P = T.root();
	const vector<bigint> Pmod = rem_tree_sq(P, T);

	vector<u64> u(k);
	for (size_t i = 0; i < k; ++i)
	{
		const u64 q = mod_u64(Pmod[i] / bigint(primes[i]), primes[i]); // (P/p_i) mod p_i = ((P mod (p_i^2)) / p_i) mod p_i
		u[i] = powmod64(q, primes[i] - 2, primes[i]); // (P/p_i)^-1 mod p_i
	}

	vector<bigint> X(d);
	for (size_t c = 0; c < d; ++c)
	{
		vector<bigint> cur(k);
		for (size_t i = 0; i < k; ++i) cur[i] = mulmod64(res[i][c], u[i], primes[i]); // res_i * (P/p_i)^-1 mod p_i

		for (size_t lv = 0; lv + 1 < T.level.size(); ++lv)
		{
			const auto& nodes = T.level[lv];
			vector<bigint> nxt;
			for (size_t i = 0; i + 1 < cur.size(); i += 2) nxt.push_back(cur[i] * nodes[i + 1] + cur[i + 1] * nodes[i]); // tree recurrence, sum = sum_l*c_r + sum_r*c_l
			if (cur.size() & 1) nxt.push_back(cur.back());
			cur = std::move(nxt);
		}

		// put x in [-P/2,P/2]
		X[c] = cur[0] % P;
		if (X[c] > P / 2) X[c] -= P;
	}
	return X;
}
