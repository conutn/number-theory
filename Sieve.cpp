#include <algorithm>
#include <iostream>
#include <unordered_set>
#include <vector>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/math/tools/polynomial.hpp>

using namespace boost::multiprecision;
using namespace boost::math;
using namespace boost::math::tools;
using namespace std;

typedef cpp_bin_float_oct float256_t;
typedef number<cpp_bin_float<512, digit_base_2>> float4096_t;
typedef polynomial<cpp_int> poly;

float256_t log10_2 = log10(float256_t(2));

#define O(A) std::cout << (A) << '\n';

template <typename T>
ostream& operator << (ostream& o, const vector<T>& v)
{
	cout << "[ ";
	for (T t : v)
	{
		cout << t << " ";
	}
	cout << "]";
	return o;
}

static vector<cpp_int> gen_prime(cpp_int lim)
{
	unordered_set<cpp_int> sieve;

	if (lim >= 2)
		sieve.insert(2);
	if (lim >= 3)
		sieve.insert(3);

	for (cpp_int x = 1; x * x < lim; ++x)
	{
		for (cpp_int y = 1; y * y < lim; ++y)
		{
			cpp_int n = 4 * x * x + y * y;
			if (n <= lim && (n % 12 == 1 || n % 12 == 5))
			{
				if (sieve.find(n) != sieve.end()) sieve.erase(n);
				else sieve.insert(n);
			}

			n = 3 * x * x + y * y;
			if (n <= lim && n % 12 == 7)
			{
				if (sieve.find(n) != sieve.end()) sieve.erase(n);
				else sieve.insert(n);
			}

			n = 3 * x * x - y * y;
			if (n <= lim && x > y && n % 12 == 11)
			{
				if (sieve.find(n) != sieve.end()) sieve.erase(n);
				else sieve.insert(n);
			}
		}
	}

	for (cpp_int r = 5; r * r <= lim; ++r)
		if (sieve.find(r) != sieve.end())
			for (cpp_int i = r * r; i <= lim; i += r * r)
				sieve.erase(i);

	vector<cpp_int> list;

	for (cpp_int i : sieve)
		list.push_back(i);

	sort(list.begin(), list.end());

	return list;
}

static cpp_int mod_pow(cpp_int base, cpp_int exp, cpp_int mod)
{
	base %= mod;
	cpp_int res = 1;

	while (exp > 0)
	{
		if (exp & 1) res = res * base % mod;
		base = base * base % mod;
		exp >>= 1;
	}

	return res;
}

static size_t _sdcnt(const cpp_int& n)
{
	return (n.backend().size() << 5) - (size_t)_lzcnt_u32(n.backend().limbs()[n.backend().size() - 1]);
}

static poly mod(const poly& a, const cpp_int& p) // coefficients of a in Z/pZ
{
	poly res(a);

	for (size_t i = 0; i < a.size(); ++i)
	{
		res[i] = a[i] % p;
	}

	return res + poly();
}

static poly mod(const poly& a, const poly& b, const cpp_int& p) // a % b in Z/pZ
{
	poly ar(a), br(b);

	while (ar.size() >= br.size())
	{
		cpp_int mult = ar[ar.size() - 1] * mod_pow(br[br.size() - 1], p - 2, p) % p;

		br = mod(br * mult, p);
		
		size_t s = ar.size() - br.size();

		br.data().insert(br.data().begin(), s, 0);
		ar -= br;
		br.data().erase(br.data().begin(), br.data().begin() + s);

		ar = mod(ar, p);
	}

	return ar;
}

static poly gcd(const poly& a, const poly& b, const cpp_int& p) // gcd(a,b) in Z/pZ
{
	poly ar = mod(a, p), br = mod(b, p);
	
	if (a.size() == 0 && b.size() == 0)
	{
		return poly();
	}

	while (ar.size() > 0 && br.size() > 0)
	{
		if (ar.size() > br.size()) ar = mod(ar, br, p);
		else br = mod(br, ar, p);
	}
	
	if (std::max(ar.size(), br.size()) <= 1) return poly();

	if (ar.size() > br.size()) return ar;
	if (br.size() > ar.size()) return br;

	return poly();
}

static poly synth_div(const poly& a, const cpp_int& n, const cpp_int& p)
{
	poly res(vector<cpp_int>(a.size() - 1));
	cpp_int nr = n % p, carry = 0;

	for (size_t j = a.size() - 1; j > 0; --j)
	{
		res[j - 1] = (a[j] + carry) % p;
		carry = (res[j - 1] * nr) % p;
	}

	return res;
}

class mod_p {
public:
	vector<cpp_int> fac;
	cpp_int p;
	mod_p(const size_t& lim, const cpp_int& _p) : fac(lim), p(_p)
	{
		fac[0] = 1;
		for (size_t i = 1; i < fac.size(); ++i)
		{
			fac[i] = fac[i - 1] * i % p;
		}
	};
	~mod_p() { fac.clear(); };

	cpp_int factorial(const size_t& n) const
	{
		if (n >= fac.size()) return 0;
		return fac[n];
	}

	cpp_int nCr(const size_t& a, const size_t& b) const
	{
		return factorial(a) * mod_pow(factorial(b) * factorial(a - b), p - 2, p);
	}
};

static poly shift(const poly& f, const cpp_int& s, const cpp_int& p) // f(x-s)
{
	mod_p m(f.size(), p);

	vector<cpp_int> pows(f.size());
	pows[0] = 1;

	for (size_t i = 1; i < pows.size(); ++i)
	{
		pows[i] = pows[i - 1] * -s % p;
	}

	poly res(f);

	for (size_t i = 0; i < res.size(); ++i)
	{
		for (size_t j = 0; j < i; ++j)
		{
			res[j] = pows[i - j] * m.nCr(i, j) * res[i] + res[j];
		}
	}

	return mod(res, p);
}

template <typename T>
static void concat(vector<T>& a, vector<T>& b)
{
	a.insert(a.begin(), b.begin(), b.end());
}

template <typename T>
static void concat(vector<T>& a, vector<T>&& b)
{
	a.insert(a.begin(), make_move_iterator(b.begin()), make_move_iterator(b.end()));
}

static vector<cpp_int> _rtfndr(const poly& a, const cpp_int& p)
{
	if (a.size() <= 1) return {};
	if (a.size() == 2) return {a[0] * mod_pow(-a[1], p - 2, p) % p};

	cpp_int power = (p - 1) / 2, offset = -1;
	vector<cpp_int> sols;
	poly res, shifted;

	do
	{
		shifted = shift(a, ++offset, p);

		if (shifted(0) % p == 0)
		{
			sols.push_back(0);
			shifted.data().erase(shifted.data().begin());
			concat(sols, _rtfndr(shifted, p));

			for (size_t i = 0; i < sols.size(); ++i)
			{
				sols[i] = (sols[i] - offset) % p;
			}

			return sols;
		}
		
		vector<poly> pows(_sdcnt(power));
		pows[0] = poly{ 0, 1 };

		for (size_t i = 1; i < pows.size(); ++i)
		{
			pows[i] = mod(pows[i - 1] * pows[i - 1], shifted, p);
		}

		res = poly{ 1 };
		cpp_int mask = 1;
		size_t sh = 0;

		while (mask <= power)
		{
			if (power & mask) res = mod(res * pows[sh], shifted, p);
			mask <<= 1, sh += 1;
		}
	} while (res.size() == 1 && ((res[0] - 1) % p == 0 || (res[0] + 1) % p == 0));
	
	res[0] += 1;

	for (size_t i = 0; i < 2; ++i)
	{
		concat(sols, _rtfndr(gcd(shifted, res, p), p));
		res[0] -= 2;
	}

	for (size_t i = 0; i < sols.size(); ++i)
	{
		sols[i] = (sols[i] - offset) % p;
	}

	return sols;
}

static vector<cpp_int> roots(const poly& f, const cpp_int& p) // f(x)=0 (mod p)
{
	if (p == 2)
	{
		vector<cpp_int> r;

		if (f(0) % 2 == 0) r.push_back(0);
		if (f(1) % 2 == 0) r.push_back(1);

		return r;
	}

	poly a = mod(f, p), res{ 1 };

	vector<poly> pows(_sdcnt(p));
	pows[0] = poly{ 0, 1 };

	for (size_t i = 1; i < pows.size(); ++i)
	{
		pows[i] = mod(pows[i - 1] * pows[i - 1], a, p);
	}

	cpp_int mask = 1;
	size_t sh = 0;

	while (mask <= p)
	{
		if (p & mask) res = mod(res * pows[sh], a, p);
		mask <<= 1, sh += 1;
	}

	res[1] -= 1;

	return _rtfndr(gcd(a, res, p), p);
}

static cpp_int norm(poly f, cpp_int a, cpp_int b) // N(a+b*theta)
{
	vector<cpp_int> a_pow(f.size()), b_pow(f.size());
	a_pow[0] = 1, b_pow[b_pow.size() - 1] = 1;

	for (size_t i = 1; i < a_pow.size(); ++i)
	{
		a_pow[i] = a_pow[i - 1] * a;
		b_pow[b_pow.size() - i - 1] = b_pow[b_pow.size() - i] * b;
	}

	cpp_int res = 0;

	for (size_t n = 0; n < f.size(); ++n)
	{
		cpp_int x = f[n] * a_pow[n] * b_pow[n];

		if ((f.size() + n) & 1)	res += x;
		else res -= x;
	}

	return res;
}

int main()
{
	cout << setprecision(std::numeric_limits<float256_t>::max_digits10);
	poly p{ 1, 1, 1, 1, 1, 1 };
	cout << roots(p, 2) << "\n";
}
