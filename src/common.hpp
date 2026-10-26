#pragma once
// common includes, namespaces and the basic types used by every module

#include <algorithm>
#include <iomanip>
#include <cmath>
#include <string>
#include <cstring>
#include <complex>
#include <map>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <unordered_set>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#ifndef GNFS_NO_GMP
#include <boost/multiprecision/gmp.hpp> // link with -lgmp; build with -DGNFS_NO_GMP to use Boost's own cpp_int instead
#endif
#include <boost/multiprecision/miller_rabin.hpp>
#include <boost/math/tools/polynomial.hpp>

// the code base uses these namespaces unqualified
using namespace boost::multiprecision;
using namespace boost::math;
using namespace boost::math::tools;
using namespace std;

// big integer type used everywhere.  GMP's mpz is much faster than cpp_int for large numbers:
// cpp_int division is quadratic (measured: 0.9s at 10^6 bits), which defeats the product/remainder trees
#ifdef GNFS_NO_GMP
typedef boost::multiprecision::cpp_int bigint;
#else
typedef boost::multiprecision::mpz_int bigint;
#endif
typedef polynomial<bigint> poly;

typedef uint64_t u64;
typedef unsigned __int128 u128;

// debugging aid: print a value
#define O(A) std::cout << (A) << '\n';

// print a vector as [ a b c ]
template <typename T>
ostream& operator << (ostream& o, const vector<T>& v)
{
	cout << "[ ";
	for (auto& t : v)
	{
		cout << t << " ";
	}
	cout << "]";
	return o;
}
