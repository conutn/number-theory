// product and remainder trees, and the Chinese remainder reconstruction built on them

#include "crt_tree.hpp"
#include "numtheory.hpp"

struct prod_tree
{
	vector<vector<bigint>> level; // level[0] = leaves, level.back() = { product of all }

	explicit prod_tree(vector<bigint> leaves)
	{
		if (leaves.empty()) leaves.push_back(1);
		level.push_back(std::move(leaves));
		while (level.back().size() > 1)
		{
			const auto& cur = level.back();
			vector<bigint> nxt;
			for (size_t i = 0; i + 1 < cur.size(); i += 2) nxt.push_back(cur[i] * cur[i + 1]);
			if (cur.size() & 1) nxt.push_back(cur.back());
			level.push_back(std::move(nxt));
		}
	}

	const bigint& root() const { return level.back()[0]; }
};

bigint balanced_product(vector<bigint> v) { return prod_tree(std::move(v)).root(); }

// value mod (leaf_i)^2 for every leaf, top-down: node modulus = (product of the leaves below)^2
static vector<bigint> rem_tree_sq(const bigint& value, const prod_tree& T)
{
	vector<bigint> cur{ value % (T.root() * T.root()) };
	for (size_t lv = T.level.size() - 1; lv-- > 0;)
	{
		const auto& nodes = T.level[lv];
		vector<bigint> nxt(nodes.size());
		for (size_t i = 0; i < nodes.size(); ++i) nxt[i] = cur[i / 2] % (nodes[i] * nodes[i]);
		cur = std::move(nxt);
	}
	return cur;
}

vector<u64> reduce_all(const bigint& value, const vector<u64>& primes)
{
	vector<u64> out;
	if (primes.empty()) return out;
	const prod_tree T(vector<bigint>(primes.begin(), primes.end()));
	const vector<bigint> r = rem_tree_sq(value, T);
	for (size_t i = 0; i < primes.size(); ++i) out.push_back(mod_u64(r[i], primes[i]));
	return out;
}

vector<bigint> crt_tree(const vector<u64>& primes, const vector<vector<u64>>& res, size_t d)
{
	const size_t k = primes.size();
	const prod_tree T(vector<bigint>(primes.begin(), primes.end()));
	const bigint& P = T.root();
	const vector<bigint> Pmod = rem_tree_sq(P, T);

	vector<u64> u(k);
	for (size_t i = 0; i < k; ++i)
	{
		const u64 q = mod_u64(Pmod[i] / bigint(primes[i]), primes[i]); // (P/p_i) mod p_i
		u[i] = powmod64(q, primes[i] - 2, primes[i]);
	}

	vector<bigint> X(d);
	for (size_t c = 0; c < d; ++c)
	{
		vector<bigint> cur(k);
		for (size_t i = 0; i < k; ++i) cur[i] = mulmod64(res[i][c], u[i], primes[i]);

		for (size_t lv = 0; lv + 1 < T.level.size(); ++lv)
		{
			const auto& nodes = T.level[lv];
			vector<bigint> nxt;
			for (size_t i = 0; i + 1 < cur.size(); i += 2) nxt.push_back(cur[i] * nodes[i + 1] + cur[i + 1] * nodes[i]);
			if (cur.size() & 1) nxt.push_back(cur.back());
			cur = std::move(nxt);
		}

		X[c] = cur[0] % P;
		if (X[c] > P / 2) X[c] -= P;
	}
	return X;
}
