cmake_minimum_required(VERSION 3.10)
project(gnfs CXX)

if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()

if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# option -DGNFS_USE_GMP=OFF to use cpp_int (much slower)
option(GNFS_USE_GMP "Use GMP for big integers" ON)

find_package(Boost REQUIRED)

# suppress false warning
add_compile_options(
    $<$<CXX_COMPILER_ID:GNU>:-Wno-free-nonheap-object>
)

add_library(gnfs_core STATIC
  src/colreduce.cpp
  src/crt_tree.cpp
  src/dump.cpp
  src/factor_base.cpp
  src/lanczos.cpp
  src/matrix.cpp
  src/numtheory.cpp
  src/polymod.cpp
  src/sieve.cpp
  src/sqrt.cpp
  src/sqrt_couveignes.cpp
  src/sqrt_newton.cpp
  src/state.cpp
)

target_include_directories(gnfs_core PUBLIC src ${Boost_INCLUDE_DIRS})

if(NOT CMAKE_VERSION VERSION_LESS 3.16)
  target_precompile_headers(gnfs_core PRIVATE src/common.hpp)
endif()

if(GNFS_USE_GMP)
  find_library(GMP_LIBRARY gmp)
  find_path(GMP_INCLUDE_DIR gmp.h)
  if(NOT GMP_LIBRARY OR NOT GMP_INCLUDE_DIR)
    message(FATAL_ERROR "GMP not found")
  endif()
  target_include_directories(gnfs_core PUBLIC ${GMP_INCLUDE_DIR})
  target_link_libraries(gnfs_core PUBLIC ${GMP_LIBRARY})
else()
  target_compile_definitions(gnfs_core PUBLIC GNFS_NO_GMP)
endif()

add_executable(sieve
  src/main.cpp
  src/planner.cpp
  src/poly_select.cpp
)
target_link_libraries(sieve PRIVATE gnfs_core)

add_executable(postprocess
  src/postprocess.cpp
)
target_link_libraries(postprocess PRIVATE gnfs_core)

add_executable(finish
  src/finish.cpp
)
target_link_libraries(finish PRIVATE gnfs_core)

enable_testing()

add_executable(gnfs_tests
  tests/test_crt_tree.cpp
)
target_include_directories(gnfs_tests PRIVATE doctest)
target_link_libraries(gnfs_tests PRIVATE gnfs_core)
add_test(NAME gnfs_tests COMMAND gnfs_tests)
