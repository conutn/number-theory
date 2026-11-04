#pragma once
// common includes, namespaces, and basic types

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <mutex>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <unordered_set>
#include <unordered_map>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#ifndef GNFS_NO_GMP
#include <boost/multiprecision/gmp.hpp> 
#endif
#include <boost/multiprecision/miller_rabin.hpp>
#include <boost/math/tools/polynomial.hpp>

// code base uses these namespaces unqualified
using namespace boost::multiprecision;
using namespace boost::math;
using namespace boost::math::tools;
using namespace std;

using std::swap;
using std::min;
using std::max;

// mpz_int is much faster than cpp_int
#ifdef GNFS_NO_GMP
typedef cpp_int bigint;
#else
typedef mpz_int bigint;
#endif
typedef polynomial<bigint> poly;

typedef uint32_t u32;
typedef uint64_t u64;
typedef unsigned __int128 u128;

// debugging
#define O(A) std::cout << (A) << '\n';

// print vector as [ a b c ]
template <typename T>
ostream& operator << (ostream& o, const vector<T>& v)
{
	o << "[ ";
	for (auto& t : v)
	{
		o << t << " ";
	}
	o << "]";
	return o;
}
