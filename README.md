# GNFS factoring program

A General Number Field Sieve implementation: planner, polynomial selection, line sieve, block Lanczos and two square-root methods.
Reads `n` on standard input, prints `n = p * q`; progress and the parameter plan go to standard error.

## Build

    cmake -S . -B build              # Release build, GMP for big integers (needs libgmp-dev and Boost headers)
    cmake --build build              # produces build/sieve

Options: `-DGNFS_USE_GMP=OFF` uses Boost's `cpp_int` (no GMP needed, but much slower on large inputs),
`-DCMAKE_CXX_STANDARD=20` builds with a newer standard (17 is the default).

## Use

    echo 100000980001501 | ./build/sieve        # everything automatic
    ./build/sieve plan < n.txt                  # only print the parameter plan and the predicted running time
    ./build/sieve [bound] [A] [max b] < n.txt   # fix the factor base bound, the sieve half-width A, the b limit

Environment variables (mostly for testing):
`GNFS_POLY=base` plain base-m polynomial, `GNFS_POLY=select` always search; `GNFS_DEGREE=d` forces the degree;
`GNFS_NONMONIC=1` (with `select`) leaves the monic candidates out.

## Layout (`src/`)

| module | role |
|---|---|
| `common.hpp` | includes, namespaces, `bigint` (GMP or `cpp_int`), `poly`, `u64` |
| `numtheory` | modular arithmetic, 64-bit mulmod/powmod, Jacobi symbol, integer roots, prime generation |
| `polymod` | polynomials over Z and F\_p: gcd, roots, factorization degrees, irreducibility test, `monic_transform` |
| `poly_select` | base-m polynomials; lattice (LLL) search for better, non-monic polynomials: `select_polynomial` |
| `planner` | model of the number of relations; chooses bound, sieve width and line count; predicts the running time |
| `factor_base` | rational / algebraic (incl. projective) / quadratic-character bases, bound heuristics, norms |
| `sieve` | relation finding: line sieve on both sides, candidate factoring, relation storage |
| `matrix` | sparse matrix over F_2, singleton and surplus-row filtering, products with B and B^T |
| `lanczos` | block Lanczos on B*B^T, dense fallback; both return dependencies among the rows |
| `crt_tree` | product and remainder trees, CRT reconstruction |
| `sqrt_newton` | square root in Z[alpha] by p-adic Newton lifting (any degree) |
| `sqrt_couveignes` | Couveignes' square root (odd degree) |
| `sqrt` | picks the method and turns a dependency into a factor |
| `main.cpp` | the driver: polynomial, plan, sieve, filter, linear algebra, square roots |

Each module is a `.hpp` (types, documentation, prototypes) and a `.cpp`. Only names used by other modules are in a header;
everything else is `static` inside its `.cpp`.

Header dependencies (no cycles):

    common
    ├── numtheory, polymod, crt\_tree, planner, poly\_select, factor\_base
    ├── sieve        -> factor\_base
    ├── matrix       -> sieve
    ├── lanczos      -> matrix
    ├── sqrt\_newton  -> sieve
    ├── sqrt\_couveignes -> factor\_base, sieve      (also uses numtheory, polymod, crt\_tree)
    └── sqrt         -> factor\_base, sieve         (calls sqrt\_newton, sqrt\_couveignes)

Data flow in `main`: `select_polynomial` -> `planner::plan_parameters` -> `make_factor_bases` -> `find_relations`
(repeated with a larger target until `filter_matrix` leaves enough rows) -> `lanczos_dependencies` (or `dense_dependencies`)
-> `try_dependency` on each dependency until one gives a factor.
