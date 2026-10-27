#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest.h>
#include "common.hpp"
#include "crt_tree.hpp"

TEST_CASE("reduce_all")
{
	bigint value{ 123454321 };
	vector<u64> primes{ 2, 3, 5, 7, 11, 13, 17, 19, 23, 29 };
	
	vector<u64> result = reduce_all(value, primes);

	REQUIRE(result == vector<u64>{ 1, 1, 1, 4, 1, 3, 15, 16, 4, 16 });

	bigint large_val{ "58944179072512134454730901998550607675832106253317606582218767307351440902373" };
	vector<u64> large_ps{ 8132933252884780463ull, 10141842613149767987ull, 2088896335063434733ull, 13315064649787534103ull };

	vector<u64> large_res = reduce_all(large_val, large_ps);

	REQUIRE(large_res == vector<u64>{ 3186221815267225710ull, 8057698696427857757ull, 431670929934640352ull, 3769940679750126006ull });
}

TEST_CASE("crt_tree")
{
	vector<u64> primes{ 2, 3, 5, 7, 11, 13, 17, 19, 23, 29 };
	vector<vector<u64>> res{
		{  0,  0,  1, 50 },
		{  0,  1,  2, 50 }, 
		{  0,  2,  4, 50 }, 
		{  0,  3,  6, 50 }, 
		{  0,  4, 10, 50 }, 
		{  0,  5, 12, 50 }, 
		{  0,  6, 16, 50 }, 
		{  0,  7, 18, 50 }, 
		{  0,  8, 22, 50 }, 
		{  0,  9, 28, 50 }
	};

	vector<bigint> result = crt_tree(primes, res, 4);

	REQUIRE(result == vector<bigint>{ 0, -703693778, -1, 50 });

	vector<u64> large_ps{ 2419808266336676179ull, 16261357598916606971ull, 16377197336239428979ull, 8464028951218056779ull };
	vector<vector<u64>> large_res{ 
		{  2127447350312611634ull, 0 }, 
		{ 13313188221646406652ull, 0 }, 
		{  7474479101033451421ull, 0 }, 
		{  5564766933607153098ull, 0 }
	};

	vector<bigint> large_result = crt_tree(large_ps, large_res, 1);

	REQUIRE(large_result == vector<bigint>{ bigint{ "-1902811939361092227462764977724395823895334874420494285437402368649519228466" } });
}
