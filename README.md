# GNFS factoring program

A General Number Field Sieve implementation: planner, polynomial selection, line sieve, block Lanczos and two square-root methods.
Reads `n` on standard input, prints `n = p * q`; progress and the parameter plan go to standard error.

## Build

    cmake -S . -B build
    cmake --build build

Options: `-DGNFS\_USE\_GMP=OFF` uses Boost's `cpp\_int` (no GMP needed)

## Use

    echo 100000980001501 | ./build/sieve      # factor a number
    ./build/sieve plan < n.txt                # predict plan and run time
    ./build/sieve [bound] [A] [max b] < n.txt # fix factor base bound, sieve half-width A, b limit

## Layout (`src/`)

| module | role |
|---|---|
| `common.hpp` | includes, namespaces, `bigint` (GMP or `cpp\_int`), `poly`, `u64` |
| `numtheory` | modular arithmetic, 64-bit mulmod/powmod, Jacobi symbol, integer roots, prime generation |
| `polymod` | polynomials over Z and F\_p: gcd, roots, factorization degrees, irreducibility test, `monic\_transform` |
| `poly\_select` | base-m polynomials; lattice (LLL) search for better, non-monic polynomials: `select\_polynomial` |
| `planner` | model of the number of relations; chooses bound, sieve width and line count; predicts the running time |
| `factor\_base` | rational / algebraic (incl. projective) / quadratic-character bases, bound heuristics, norms |
| `sieve` | relation finding: line sieve on both sides, candidate factoring, relation storage |
| `matrix` | sparse matrix over F\_2, singleton and surplus-row filtering, products with B and B^T |
| `lanczos` | block Lanczos on B\*B^T, dense fallback; both return dependencies among the rows |
| `crt\_tree` | product and remainder trees, CRT reconstruction |
| `sqrt\_newton` | square root in Z[alpha] by p-adic Newton lifting (any degree) |
| `sqrt\_couveignes` | Couveignes' square root (odd degree) |
| `sqrt` | picks the method and turns a dependency into a factor |
| `main.cpp` | the driver: polynomial, plan, sieve, filter, linear algebra, square roots |

Workflow in `main`: `select\_polynomial` -> `planner::plan\_parameters` -> `make\_factor\_bases` -> `find\_relations`
-> `lanczos\_dependencies` (or `dense\_dependencies`) -> `try\_dependency`
