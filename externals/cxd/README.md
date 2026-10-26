# CXD

CXD is a C++20 library for two-dimensional and three-dimensional convex
collision detection. It implements the Gilbert-Johnson-Keerthi (GJK)
algorithm for intersection tests and the Expanding Polytope Algorithm
(EPA) for penetration information. A Separating Axis Theorem (SAT)
narrowphase is also provided for small convex meshes (n <= 8 vertices).

The library is adapter-based and supports floating-point vectors from the
C++ standard library (`std::vector`, `std::array`), Boost.uBLAS, and Eigen.
The standalone benchmarks additionally exercise the Siconos vector type.

## Quick start

```cpp
#include "cxd.hpp"
#include "cxd_std.hpp"  // std::vector adapter

// Two shapes as collections of vertices
std::vector<std::vector<double>> a = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};
std::vector<std::vector<double>> b = {{2, 2}, {4, 2}, {4, 4}, {2, 4}};

// GJK intersection test
auto [intersects, simplex0, simplex1, simplex2] =
    cxd::gjk_inter<2>(a, b, {1.0, 0.0});

// EPA penetration (depth, normal, contact vertices)
auto [depth, normal, v0, v1] =
    cxd::epa<2>(a, b, simplex0, simplex1, simplex2, simplex2, 1e-6);
```

## Headers

| Header                      | Purpose                                           |
|-----------------------------|---------------------------------------------------|
| `src/cxd.hpp`               | Core scalar GJK, EPA, SAT APIs and types          |
| `src/cxd_batch.hpp`         | Batched GJK/EPA/SAT over multiple pairs           |
| `src/cxd_std.hpp`           | `std::vector` adapter (also used for Broad phase) |
| `src/cxd_std_array_vec.hpp` | `std::array` fixed-size vector adapter            |
| `src/cxd_ublas.hpp`         | Boost.uBLAS adapter                               |

Every adapter provides `env_t<Point>` with `scalar`, `vector`, `support`,
and `inner_prod` traits. Include the adapter header *before* `cxd.hpp`
(or `cxd_batch.hpp`) so the point type is visible at definition time.

## Collision API (scalar)

### `gjk_inter<Dim>`

```cpp
template<std::size_t Dim, typename ShapesA, typename ShapesB, typename Vector>
constexpr auto gjk_inter(const ShapesA& sa, const ShapesB& sb,
                         const Vector& init,
                         const gjk_options options = {})
```

Tests two shapes for intersection using GJK.

- **`Dim`** must be 2 or 3.
- **`sa`, `sb`**: non-empty collections of floating-point vectors (2D or 3D).
- **`init`**: a non-zero initial search direction.
- **`options`**: optional `gjk_options` (`tolerance = 1e-6`,
  `max_iterations = 64`).

Returns `std::tuple<bool, support, support, support>` (2D) or
`stdtuple<bool, support, support, support, support>` (3D):

- **bool**: `true` if the origin is enclosed (shapes intersect),
  `false` if separated.
- **support**: the GJK simplex vertices (3 for 2D, up to 4 for 3D).

Throws `std::invalid_argument` for empty shapes, wrong dimension, zero
direction, or non-positive tolerance/iterations. Throws
`std::runtime_error` if the iteration limit is exceeded. A repeated
support point reports separation (no further progress possible).

### `epa`

```cpp
template<std::size_t Dim, typename ShapesA, typename ShapesB, typename Support>
constexpr auto epa(const ShapesA& sa, const ShapesB& sb,
                   const Support& supporta, const Support& supportb,
                   const Support& supportc, const Support& supportd,
                   const any_scalar auto& tol,
                   const std::size_t max_iterations = 64)
```

Expands the complete Minkowski-difference polytope from a non-degenerate
GJK simplex to find the deepest penetration.

- **`supporta`–`supportd`**: the simplex vertices from `gjk_inter`.
- **`tol`**: positive tolerance for convergence.
- **`max_iterations`**: positive iteration budget.

Returns `std::tuple<scalar, vector, support, support>` (2D) or
`std::tuple<scalar, vector, support, support, support>` (3D):

- **depth**: penetration depth (positive when enclosed).
- **normal**: contact normal oriented A -> B.
- **support**: the vertices of the deepest face (2 for 2D, 3 for 3D).

### `broad_phase`

```cpp
auto broad_phase(const std::vector<Shapes>& shapes, double padding = 0.0)
```

Sweep-and-prune over axis-aligned bounding boxes. Returns index pairs whose
boxes overlap (including edge-touching). Pairs not in the result are
definitely separated. The dimension is parameterized by the shape vertex
count (2D for `std::array<double, 2>`, 3D for `std::array<double, 3>`).

### SAT (scalar)

```cpp
template<std::size_t Dim, typename ShapesA, typename ShapesB>
constexpr auto sat_inter(const ShapesA& sa, const ShapesB& sb,
                         const gjk_options options = {})

template<std::size_t Dim, typename ShapesA, typename ShapesB>
constexpr auto sat(const ShapesA& sa, const ShapesB& sb,
                   const gjk_options options = {})

template<std::size_t Dim, typename ShapesA, typename ShapesB>
constexpr auto sat_distance(const ShapesA& sa, const ShapesB& sb,
                            const gjk_options options = {})
```

SAT-based narrowphase for small convex meshes. `sat_inter` returns a
tuple `(bool intersects, normal, h)`. `sat` returns `(depth, normal,
vertices)`. `sat_distance` returns the signed distance (positive when
separated, negative when enclosed).

SAT is standalone: it is not automatically wired into `gjk_batch` because
SAT and GJK can disagree on degenerate rungs (where GJK returns
`max_iterations`).

## Batch API

The batch API processes multiple pairs in a single call using packed
vertex pools. It is designed for the per-step narrowphase of a
constraint solver (e.g., Siconos vkernel).

### `packed_pool<Point>`

A single buffer holding every shape's vertices, plus one
`shape_entry{vertex_offset, vertex_count}` per shape:

```cpp
template<typename Point>
struct packed_pool {
    std::vector<Point> vertex_buffer;     // all vertices concatenated
    std::vector<shape_entry> entries;     // one per shape
    using vector = Point;
    using scalar = typename traits::env_t<Point>::scalar;
    using support = typename traits::env_t<Point>::support;
};
```

Build one with `pack_shapes<Dim>(shapes)` or from two shape lists with
the 3D convenience overload.

### `pair_record<Scalar, Dim>`

```cpp
template<typename Scalar, std::size_t Dim>
struct pair_record {
    std::uint32_t shape_id_A;
    std::uint32_t shape_id_B;
    Scalar init_dir[Dim];   // initial GJK search direction
};
```

### `gjk_batch`

```cpp
template<std::size_t Dim, typename Pool, typename PairRecord, typename Support>
void gjk_batch(const Pool& pool,
               const std::vector<PairRecord>& pairs,
               const gjk_options options,
               std::vector<batch_result<Support, Dim + 1>>& results,
               const pool_bounds<typename Pool::scalar, Dim>* bounds = nullptr,
               batch_timing* timing = nullptr)
```

- **`results`**: caller-owned buffer, reused across calls (`assign` refills
  to the pair count, no stale state).
- **`bounds`**: optional circumradius prefilter (see below).
- **`timing`**: optional per-pair timing (null = zero cost).

An allocating convenience overload returns the result vector directly.

A 3D convenience overload accepts two shape lists and a shared initial
direction, packing them into one pool with a zipped pair list.

### `epa_batch`

```cpp
template<std::size_t Dim, typename Pool, typename PairRecord, typename Support>
void epa_batch(const Pool& pool,
               const std::vector<PairRecord>& pairs,
               const std::vector<std::uint32_t>& hit_pairs,
               const std::vector<batch_result<Support, Dim + 1>>& gjk_results,
               const gjk_options options,
               std::vector<epa_result<typename Pool::scalar, Dim>>& answers,
               epa_arena<Pool, Dim>& arena,
               batch_timing* timing = nullptr)
```

Runs EPA only over the enclosed pairs identified by `hit_pairs` (indices
into `gjk_results`). The `epa_arena` provides fixed-capacity storage reused
across calls. Returns `epa_result<Scalar, Dim>` rows:

```cpp
template<typename Scalar, std::size_t Dim>
struct epa_result {
    Scalar depth;
    std::array<Scalar, Dim> normal;     // A -> B
    std::array<std::pair<std::size_t, std::size_t>, Dim> vertices{};
    // vertex index pairs (idx in A, idx in B) of the deepest face
};
```

### `sat_batch` / `sat_batch_pair`

```cpp
template<std::size_t Dim, typename Pool, typename PairRecord>
void sat_batch(const Pool& pool,
               const std::vector<PairRecord>& pairs,
               const gjk_options options,
               std::vector<sat_result<typename Pool::scalar, Dim>>& answers,
               const pool_bounds<typename Pool::scalar, Dim>* bounds = nullptr,
               batch_timing* timing = nullptr)

template<std::size_t Dim, typename Pool, typename PairRecord>
auto sat_batch_pair(const Pool& pool, const PairRecord& pair,
                    const gjk_options options,
                    const pool_bounds<typename Pool::scalar, Dim>* bounds = nullptr)
    -> sat_result<typename Pool::scalar, Dim>
```

Batched SAT narrowphase returning `sat_result` rows:

```cpp
template<typename Scalar, std::size_t Dim>
struct sat_result {
    status st;                        // separated, enclosed, invalid_input
    std::array<Scalar, Dim> normal{}; // oriented A -> B
    Scalar h{};                       // +gap (separated), -depth (enclosed)
};
```

A 3D shapes convenience overload packs two shape lists into a pool.

### `sat_pair`

```cpp
template<std::size_t Dim, typename ShapeA, typename ShapeB>
auto sat_pair(const ShapeA& sa, const ShapeB& sb,
              const gjk_options options = {})
    -> sat_result<scalar, Dim>
```

Scalar SAT answer for a single pair. Wraps `sat_inter` via `sat_from_tuple`.
The 2D scan requires vertices in convex cyclic order; unordered input
falls back to hull-building.

### `pool_bounds` / circumradius prefilter

```cpp
template<typename Scalar, std::size_t Dim>
struct pool_bounds {
    std::vector<Scalar> centre;        // Dim per shape
    std::vector<Scalar> radius;        // circumradius per shape
    std::vector<std::uint8_t> convex;  // convexity flag per shape
};

template<std::size_t Dim, typename Pool>
void compute_pool_bounds(const Pool& pool,
                         pool_bounds<typename Pool::scalar, Dim>& out)
```

Computes per-shape centroid and circumradius in O(total_vertices), once
per batch step. When passed to `gjk_batch` / `sat_batch` / `gjk_batch_pair`
/ `sat_batch_pair`, pairs whose centre distance exceeds `rA + rB + tolerance`
are answered as `separated` in O(1) per pair — before any simplex or SAT
work. Separated pairs carry a real normal (centre-to-centre direction) and
`h = distance - (rA + rB)` (a lower bound on the true gap).

For 2D shapes, `compute_pool_bounds` also marks convexity via a cyclic
winding test: full 2pi total turn with consistent cross-sign. Sign
consistency alone is insufficient (a golden-angle star walk turns
uniformly while crossing itself), so the total-turn closure is enforced.

### `batch_timing`

```cpp
struct batch_timing {
    std::vector<double> seconds;  // one row per pair (per hit for EPA)
};
```

Optional per-pair wall-clock breakdown, zero cost when the pointer is null.

### `status`

```cpp
enum class status { separated, enclosed, max_iterations, invalid_input };
```

GJK batch results use all four. SAT never returns `max_iterations` (SAT
decides per axis, not iteratively).

## Device backend

The batch cores are designed for accelerator porting. `assert_device_contract`
statically verifies that pool vertices and `shape_entry` are trivially
copyable. `packed_pool` buffers travel flat (one `Pool::vector` per vertex,
one `shape_entry` per shape); only the pool object itself is a host handle.

A CUDA backend (`src/device/cxd_batch_device.cuh`) provides thread-per-pair
kernels with golden parity against the host cores.

## Building

### Requirements

- C++20 compiler (GCC 14+ or Clang 21+)
- CMake 3.10+
- CppUnit (for tests)
- Boost (headers)
- Eigen 3.3+ (optional, for some adapters and benchmarks)
- Google Benchmark (for benchmarks only)
- Siconos (for the full benchmark suite; not needed for library or tests)

### Guix (recommended)

```sh
guix shell gcc-toolchain cppunit boost cmake make -- \
  cmake -B build -DCMAKE_BUILD_TYPE=Debug -DCXD_BUILD_BENCHMARKS=OFF .
make -C build -j$(nproc)
ctest --test-dir build --output-on-failure
```

### Quick CMake

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
ctest --test-dir build --output-on-failure
```

## Testing

The test suite (68 tests) covers:

- GJK intersection (2D/3D, scalar and batch, multiple adapters)
- EPA penetration (2D/3D, including degenerate seeds and adversarial cases)
- SAT narrowphase (2D/3D, pair and batch, prefilter)
- Broad phase (sweep and kd-tree)
- Batch parity (golden parity against committed scalar references)
- Buffer reuse (zero-allocation across steps)

Run with: `ctest --test-dir build --output-on-failure`

## Benchmarks

`benchmark/bench_narrowphase` times the batch cores on the Siconos vkernel
report distribution. Measurements (Release, gcc-15.2 -O2, N bodies = 907 / 3811 /
15622 pairs):

Batched SAT vs batched GJK per pair (N=250 / 1000 / 4000):

| benchmark                            | N=250 | N=1000 | N=4000      |
|--------------------------------------|-------|--------|-------------|
| `gjk_batch`                          | 68.3  | 78.1   | 79.2        |
| `gjk_batch` + circumradius prefilter | 61.9  | 61.8   | 63.5 (-20%) |
| `sat_batch` (no prefilter)           | 268.8 | 269.7  | 217.9       |
| `sat_batch` + circumradius prefilter | 90.9  | 87.2   | 87.5 (-60%) |
| `pairloop_gjk_pair`                  | 75.0  | 74.4   | 74.6        |
| `chain_gjk+epa`                      | 93.7  | 94.5   | 92.8        |

At N=4000, `sat_batch` (217.9 ns/pair) is ~2.8x slower than `gjk_batch`
(79.2 ns/pair): SAT checks all edge normals (8 projections at 4-vertex
squares) vs GJK's 2-3 support queries. The circumradius prefilter answers
~86.7% separated pairs in O(1) per pair, dropping GJK cost by 20% and SAT
cost by 60%.

Vertex-count sweep (regular n-gons, same circumradius, N=4000, 15622 pairs):

| n (verts) | encl % | `gjk_batch` ns/pair | SAT reference ns/pair | ratio              |
|-----------|--------|---------------------|-----------------------|--------------------|
| 4         | 13.4%  | 79.4                | 66.9                  | SAT 0.84x (faster) |
| 8         | 39.0%  | 102.7               | 336.5                 | GJK 3.3x faster    |
| 16        | 42.7%  | 164.5               | 813.7                 | GJK 5.0x faster    |

SAT checks 2n edge normals per pair (O(n) per pair), so its cost grows
quadratically with vertex count. GJK's cost grows ~linearly (2-3 support
queries per iteration, each O(n)). The crossover lands between n=4 and n=8
on bare cores: SAT wins at n=4, GJK wins at n>=8. SAT is best for small n
(n <= 8); at n=16 it is 5x more expensive than batch GJK.

## License

Apache 2.0. CXD is part of the Siconos Software. Copyright 2026 INRIA.
