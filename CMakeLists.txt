cmake_minimum_required(VERSION 3.10)
project(gnfs CXX)

# C++17 is enough; pass -DCMAKE_CXX_STANDARD=20 (or 23) to build with a newer standard
if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# GMP is used for all big integers by default (Boost's cpp_int has quadratic division, which
# is far too slow for the large numbers in the square root).  Configure with -DGNFS_USE_GMP=OFF
# to use cpp_int instead (no GMP needed, but much slower for big inputs).
option(GNFS_USE_GMP "Use GMP (via Boost.Multiprecision) for big integers" ON)

find_package(Boost REQUIRED)

add_executable(sieve
  src/main.cpp
  src/numtheory.cpp
  src/polymod.cpp
  src/poly_select.cpp
  src/planner.cpp
  src/factor_base.cpp
  src/sieve.cpp
  src/matrix.cpp
  src/lanczos.cpp
  src/crt_tree.cpp
  src/sqrt_newton.cpp
  src/sqrt_couveignes.cpp
  src/sqrt.cpp
)
target_include_directories(sieve PRIVATE src ${Boost_INCLUDE_DIRS})

# every file includes the same heavy Boost headers (via common.hpp): precompile them once
if(NOT CMAKE_VERSION VERSION_LESS 3.16)
  target_precompile_headers(sieve PRIVATE src/common.hpp)
endif()

if(GNFS_USE_GMP)
  find_library(GMP_LIBRARY gmp)
  find_path(GMP_INCLUDE_DIR gmp.h)
  if(NOT GMP_LIBRARY OR NOT GMP_INCLUDE_DIR)
    message(FATAL_ERROR "GMP not found (install libgmp-dev) or configure with -DGNFS_USE_GMP=OFF")
  endif()
  target_include_directories(sieve PRIVATE ${GMP_INCLUDE_DIR})
  target_link_libraries(sieve PRIVATE ${GMP_LIBRARY})
else()
  target_compile_definitions(sieve PRIVATE GNFS_NO_GMP)
endif()
