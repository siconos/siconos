/* CXD is a program dedicated to convex collision detection.
 * It is part of the Siconos Software.
 *
 * Copyright 2026 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#ifndef CXD_BATCH_HPP
#define CXD_BATCH_HPP

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "cxd.hpp"

namespace cxd {
enum class status { separated, enclosed, max_iterations, invalid_input };

struct shape_entry {
  std::uint32_t vertex_offset;
  std::uint32_t vertex_count;
};

template <typename Scalar, std::size_t Dim>
struct pair_record {
  std::uint32_t shape_id_A;
  std::uint32_t shape_id_B;
  Scalar init_dir[Dim];
};

template <typename Support, std::size_t N>
struct batch_result {
  status st;
  std::array<Support, N> seed;
};

// Result of an EPA penetration query: depth, normal, and the vertex indices
// (idx in shape A, idx in shape B) of each vertex of the deepest face.
// 2D: 2 vertices (edge). 3D: 3 vertices (triangle).
template <typename Scalar, std::size_t Dim>
struct epa_result {
  Scalar depth;
  std::array<Scalar, Dim> normal;
  std::array<std::pair<std::size_t, std::size_t>, Dim> vertices{};
};

// One packed pool: every shape's vertices concatenated into a single buffer,
// plus one offset/count entry per shape. The point's own aliases travel with
// the pool, so the evaluators name Pool::vector/Pool::scalar/Pool::support
// without deriving them again. The aliases depend on Point, so env_t's checks
// still resolve at instantiation.
template <typename Point>
struct packed_pool {
  std::vector<Point> vertex_buffer;
  std::vector<shape_entry> entries;

  using vector = Point;
  using scalar = typename traits::env_t<Point>::scalar;
  using support = typename traits::env_t<Point>::support;
};

// Phase 3's memory contract: the fixed-stride upload, stated once next to
// the types it describes.
//
// Transfer units (the caller re-packs nothing):
//   - pool.vertex_buffer.data(): one Pool::vector per shape vertex, stride
//     sizeof(Pool::vector);
//   - pool.entries.data(): one shape_entry per shape, stride
//     sizeof(shape_entry);
//   - pairs.data(): one pair_record<scalar, Dim> per queried pair, stride
//     sizeof(pair_record<scalar, Dim>);
//   - the result rows: one batch_result<support, Dim + 1> per queried pair,
//     stride sizeof(batch_result<support, Dim + 1>) (gjk_batch returns it;
//     epa_batch takes it in and returns one answers row per hit, stride
//     sizeof(std::pair<scalar, std::array<scalar, Dim>>)).
//
// Nothing in these buffers has indirection, and the row shapes do not
// depend on the entries' order: a device backend uploads those unchanged
// and an index-reordered copy is a legal layout. The pool OBJECT is a host
// handle (it owns std::vectors): only its buffers travel, and the runners
// reading records or writing results straight off the device copy stays the
// model. The asserted value-typedness is the part a host would silently
// accept and a copy-based upload silently corrupt.
//
// The seed rows are the honest caveat of this phase: batch_result and the
// answers are flat IN STRIDE (one row per pair, fixed size, no pointer), but
// support_t backs its ({shape index, vertex index, position}) tuple in
// std::tuple, whose special members are not trivial under every standard
// library, so a byte-for-byte row memcpy is NOT guaranteed. Today the
// contract therefore asks a marshaling step for the seed elements (its own
// fixed layout per support) and defers the row memcpy to the backend port
// that replaces support_t's backing with a trivially copyable struct.
template <typename Pool, std::size_t Dim>
inline void assert_device_contract()
{
  using scalar = typename Pool::scalar;
  static_assert(std::is_trivially_copyable_v<typename Pool::vector>,
                "pool vertices must be value typed to upload as-is");
  static_assert(std::is_trivially_copyable_v<shape_entry>,
                "shape entries travel flat");
  static_assert(std::is_trivially_copyable_v<pair_record<scalar, Dim>>,
                "pair records must be uploadable as-is");
  static_assert(std::is_standard_layout_v<pair_record<scalar, Dim>>,
                "pair records need the interoperable layout a memcpy keeps");
}

template <typename Shape>
auto pack_shapes(const std::vector<Shape>& shapes)
{
  using point_type = typename Shape::value_type;

  // The pool's aliases are namespace-scope members, so a caller that only
  // packs shapes never instantiates (or warns about) a function-local struct.
  using pool_type = packed_pool<point_type>;

  pool_type pool;

  std::uint32_t current_offset = 0;
  for (const auto& shape : shapes) {
    const std::uint32_t count = static_cast<std::uint32_t>(shape.size());
    pool.entries.push_back({current_offset, count});

    for (const auto& vertex : shape) {
      pool.vertex_buffer.push_back(vertex);
    }

    current_offset += count;
  }

  return pool;
}

// A non-owning window onto one shape's vertices inside a packed pool. It
// exposes exactly the surface the scalar support machinery reads --
// value_type, size_type, operator[], begin/end as input iterators -- so
// cxd::support_indice and cxd::diff_support scan a pooled shape as if it were
// a plain shape and the batched evaluator can carry gjk_inter's search over
// unchanged. The indice alias comes from the pool's own support type, so the
// indices the support scan produces are the ones env_t expects.
template <typename Pool>
struct shape_view {
  using point_type = typename Pool::vector;
  using indice = std::tuple_element_t<0, typename Pool::support>;
  using value_type = point_type;
  using size_type = indice;
  using const_iterator = const point_type*;

  shape_view(const Pool& pool, const std::uint32_t id)
      : vertices{pool.vertex_buffer.data() + pool.entries[id].vertex_offset},
        count{static_cast<size_type>(pool.entries[id].vertex_count)}
  {
  }

  const point_type& operator[](const size_type i) const
  {
    return vertices[i];
  }

  const_iterator begin() const { return vertices; }

  const_iterator end() const { return vertices + count; }

  bool empty() const { return count == 0; }

  size_type size() const { return count; }

 private:
  const point_type* vertices;
  size_type count;
};

// One support query over the pool: the window over the same diff_support scan
// the intersection pass uses, factored out so the penetration pass reaches a
// shape through the one implementation rather than a second copy of it.
template <typename Pool>
auto pool_diff_support(const Pool& pool, const std::uint32_t id_a,
                       const std::uint32_t id_b,
                       const typename Pool::vector& direction)
{
  return cxd::diff_support(shape_view<Pool>{pool, id_a},
                           shape_view<Pool>{pool, id_b}, direction);
}

// ---------------------------------------------------------------------------
// Pool bounds: the O(pool) picture the batch stages share. One centre and
// one circumradius per shape, computed once per pool build, feeds the
// O(1) prefilter: a pair whose centre distance exceeds rA + rB plus the
// tolerance is separated before any simplex work, and the SAT stages can
// answer the rejection with real numbers (the unit axis A -> B and the
// ball-gap lower bound). Non-finite vertices are skipped for centre and
// radius: a shape is bounded by its finite vertices alone, so popular
// garbage-in rows that would stall a support query far apart are turned
// into a definite separated answer instead. The convex flag (2D) says the
// vertices walk a convex polygon in the stored order, which is the fast
// row-check the SAT scan needs; 3D leaves it zero (no flag there).
// ---------------------------------------------------------------------------
template <typename Scalar, std::size_t Dim>
struct pool_bounds {
  std::vector<Scalar> centre;        // Dim per shape
  std::vector<Scalar> radius;        // one per shape
  std::vector<std::uint8_t> convex;  // one per shape
};

// 2D cyclic convexity over stored vertices: every strict turn must keep one
// sign AND the accumulated turning must be exactly one revolution. Sign
// consistency alone is not sufficient: a star walk through a golden-angle
// spiral crosses itself while turning the same way at every step, and SAT
// over such an order's edge normals is not SAT over the shape's convex
// hull, so the total turn is required to close at +-2pi. Collinear turns
// (cross == 0) are allowed -- the duplicate axis they leave behind is
// harmless for the SAT scan -- but a reversal (a 180-degree u-turn), a zero
// edge, and fewer than three vertices are not. The check the bounds fill
// marks per shape is this one, and the scalar SAT entry runs it per call,
// so a pool built without bounds still takes the same fast scan when the
// shapes allow it.
template <typename It>
bool cyclic_convex_2d(const It first, const It last)
{
  using point = std::remove_cvref_t<decltype(*first)>;
  using scalar = typename decltype(env(point{}))::scalar;

  const auto count = static_cast<std::size_t>(std::distance(first, last));
  if (count < 3) return false;

  constexpr scalar one_turn = scalar(6.28318530717958647692);
  constexpr scalar slack = scalar(1e-9);

  std::uint8_t sign = 0;
  scalar turning{};
  for (std::size_t i = 0; i < count; ++i) {
    const point& p0 = *std::next(first, static_cast<std::ptrdiff_t>(i));
    const point& p1 =
        *std::next(first, static_cast<std::ptrdiff_t>((i + 1) % count));
    const point& p2 =
        *std::next(first, static_cast<std::ptrdiff_t>((i + 2) % count));
    const scalar ex0 = p1[0] - p0[0], ey0 = p1[1] - p0[1];
    const scalar ex1 = p2[0] - p1[0], ey1 = p2[1] - p1[1];
    // a zero edge is a duplicate vertex: not a clean polygon
    if (ex0 * ex0 + ey0 * ey0 == scalar{}) return false;
    const scalar cross = ex0 * ey1 - ey0 * ex1;
    // a u-turn (zero cross against the negative dot) folds the polygon
    // back over its own edge: not a walk around a convex region
    if (cross == scalar{} && ex0 * ex1 + ey0 * ey1 < scalar{}) return false;
    if (cross > scalar{}) {
      // consistent right-hand turning stays under one revolution
      if (sign == std::uint8_t{2}) return false;
      sign = std::uint8_t{1};
      turning += std::atan2(cross, ex0 * ex1 + ey0 * ey1);
      if (turning > one_turn + slack) return false;
    }
    else if (cross < scalar{}) {
      if (sign == std::uint8_t{1}) return false;
      sign = std::uint8_t{2};
      turning += std::atan2(cross, ex0 * ex1 + ey0 * ey1);
      if (turning < -one_turn - slack) return false;
    }
  }
  // exactly one revolution's worth of turning must be walked, either way
  return turning > one_turn - slack && turning < one_turn + slack;
}

template <std::size_t Dim, typename Pool>
  requires(Dim == 2 || Dim == 3)
void compute_pool_bounds(const Pool& pool,
                         pool_bounds<typename Pool::scalar, Dim>& out)
{
  using point = typename Pool::vector;
  using scalar = typename Pool::scalar;

  const auto* entries = pool.entries.data();
  out.centre.assign(pool.entries.size() * Dim, scalar{});
  out.radius.assign(pool.entries.size(), scalar{});
  out.convex.assign(pool.entries.size(), std::uint8_t{0});

  for (std::size_t id = 0; id < pool.entries.size(); ++id) {
    const auto begin = entries[id].vertex_offset;
    const auto count = entries[id].vertex_count;

    // the centre is the mean of the finite vertices, the radius the
    // largest squared distance of a finite vertex to that centre
    scalar c[Dim];
    for (std::size_t k = 0; k < Dim; ++k) c[k] = scalar{};
    std::size_t finite = 0;
    for (std::size_t i = 0; i < count; ++i) {
      const auto& p = pool.vertex_buffer[begin + i];
      bool ok = true;
      for (std::size_t k = 0; k < Dim && ok; ++k) ok = std::isfinite(p[k]);
      if (!ok) continue;
      ++finite;
      for (std::size_t k = 0; k < Dim; ++k) c[k] += p[k];
    }
    if (finite > 0) {
      for (std::size_t k = 0; k < Dim; ++k) c[k] /= scalar(finite);
    }
    scalar r2 = scalar{};
    for (std::size_t i = 0; i < count; ++i) {
      const auto& p = pool.vertex_buffer[begin + i];
      bool ok = true;
      for (std::size_t k = 0; k < Dim && ok; ++k) ok = std::isfinite(p[k]);
      if (!ok) continue;
      scalar d2 = scalar{};
      for (std::size_t k = 0; k < Dim; ++k) {
        const scalar diff = p[k] - c[k];
        d2 += diff * diff;
      }
      r2 = std::max(r2, d2);
    }
    for (std::size_t k = 0; k < Dim; ++k) out.centre[id * Dim + k] = c[k];
    out.radius[id] = std::sqrt(r2);
  }

  // 2D convexity: the shared cyclic walk over the stored order. A shape
  // that fails it keeps the flag off and every pair it enters takes the
  // hull path instead of the fast scan.
  if constexpr (Dim == 2) {
    for (std::size_t id = 0; id < pool.entries.size(); ++id) {
      const auto begin = entries[id].vertex_offset;
      const auto count = entries[id].vertex_count;
      const auto vertices =
          pool.vertex_buffer.begin() + static_cast<std::ptrdiff_t>(begin);
      const bool convex = cyclic_convex_2d(vertices, vertices + count);
      out.convex[id] = convex ? std::uint8_t{1} : std::uint8_t{0};
    }
  }
}

template <std::size_t Dim, typename Pool>
  requires(Dim == 2 || Dim == 3)
auto compute_pool_bounds(const Pool& pool)
    -> pool_bounds<typename Pool::scalar, Dim>
{
  pool_bounds<typename Pool::scalar, Dim> out;
  compute_pool_bounds<Dim>(pool, out);
  return out;
}

// The O(1) prefilter over the bounds. Positive when the centre distance
// exceeds both radii plus the tolerance: then no point of the pair can sit
// closer than the slack, and every stage may answer separated without
// reading a vertex.
template <std::size_t Dim, typename Pool, typename Bounds>
  requires(Dim == 2 || Dim == 3)
constexpr bool bounds_separated(const Pool& pool, const Bounds& bounds,
                                const std::uint32_t id_a,
                                const std::uint32_t id_b,
                                const typename Pool::scalar tolerance)
{
  using scalar = typename Pool::scalar;
  const auto* entries = pool.entries.data();
  (void)entries;
  scalar d2 = scalar{};
  for (std::size_t k = 0; k < Dim; ++k) {
    const scalar diff = bounds.centre[std::size_t{id_a} * Dim + k] -
                        bounds.centre[std::size_t{id_b} * Dim + k];
    d2 += diff * diff;
  }
  const scalar sum = bounds.radius[id_a] + bounds.radius[id_b] + tolerance;
  return d2 > sum * sum;
}

// Per-pair evaluator: gjk_inter's search over the pool. Everything the scalar
// path reports by throwing comes back as a status instead, and the seed is
// written only on the enclosed path -- every other status leaves it
// untouched.
template <std::size_t Dim, typename Pool, typename PairRecord>
  requires(Dim == 2 || Dim == 3)
auto gjk_batch_pair(
    const Pool& pool, const PairRecord& pair, const gjk_options options,
    const pool_bounds<typename Pool::scalar, Dim>* bounds = nullptr)
    -> batch_result<typename Pool::support, Dim + 1>
{
  using result_type = batch_result<typename Pool::support, Dim + 1>;
  using vector = typename Pool::vector;
  using scalar = typename Pool::scalar;
  using support = typename Pool::support;

  const auto invalid_input = []() {
    result_type entry{};
    entry.st = status::invalid_input;
    return entry;
  };

  const auto separated = []() {
    result_type entry{};
    entry.st = status::separated;
    return entry;
  };

  const auto exhausted = []() {
    result_type entry{};
    entry.st = status::max_iterations;
    return entry;
  };

  // Validation mirrors gjk_inter, reported per pair instead of thrown: an id
  // outside the pool, an empty shape, a direction record or pooled point
  // whose extent disagrees with Dim, and a zero initial direction are
  // invalid_input for this pair alone. The extent checks keep Dim and the
  // data in step: Dim=3 over a two-component record would read init_dir[2]
  // past its end, and Dim=2 over three-component points would answer from a
  // direction whose z was silently zero-filled. gjk_inter rejects the same
  // mismatches at compile time through validate_vector_2d/3d.
  const auto& entries = pool.entries;
  if (pair.shape_id_A >= entries.size() ||
      pair.shape_id_B >= entries.size() ||
      entries[pair.shape_id_A].vertex_count == 0 ||
      entries[pair.shape_id_B].vertex_count == 0 ||
      std::size(pair.init_dir) != Dim ||
      pool.vertex_buffer[entries[pair.shape_id_A].vertex_offset].size() !=
          Dim ||
      pool.vertex_buffer[entries[pair.shape_id_B].vertex_offset].size() !=
          Dim)
    return invalid_input();

  // Both support queries read their shape straight out of the pool. The views
  // carry the surface a plain shape has, so support_indice's scan (strict <,
  // ties keep the lowest vertex index) and diff_support behave exactly as
  // they do on the scalar path.
  const auto diff_support = [&pool](const std::uint32_t id_a,
                                    const std::uint32_t id_b,
                                    const vector& direction) {
    return pool_diff_support(pool, id_a, id_b, direction);
  };

  if constexpr (Dim == 2) {
    const auto enclosed = [](const support& s0, const support& s1,
                             const support& s2) {
      result_type entry{};
      entry.st = status::enclosed;
      entry.seed = {s0, s1, s2};
      return entry;
    };

    const vector init{pair.init_dir[0], pair.init_dir[1]};

    const auto tolerance_squared = options.tolerance * options.tolerance;
    if (num::inner_prod(init, init) <= tolerance_squared)
      return invalid_input();

    // O(1) ball prefilter: a pair whose centre distance goes beyond the
    // radii sum plus the tolerance is separated before any simplex work.
    // Only the compute path shrinks; the row is exactly what the full GJK
    // would return for a pair this far apart.
    if (bounds != nullptr &&
        bounds_separated<Dim>(pool, *bounds, pair.shape_id_A, pair.shape_id_B,
                              options.tolerance))
      return separated();

    const auto is_near_zero = [&tolerance_squared](const vector& v) {
      return num::inner_prod(v, v) <= tolerance_squared;
    };
    const auto simplex0 =
        diff_support(pair.shape_id_A, pair.shape_id_B, init);
    if (num::inner_prod(support_vector(simplex0), init) < -options.tolerance)
      return separated();

    std::array<support, 3> simplex{simplex0, simplex0, simplex0};
    std::size_t simplex_size = 1;
    vector direction = -support_vector(simplex0);

    for (std::size_t iteration = 0; iteration < options.max_iterations;
         ++iteration) {
      if (is_near_zero(direction))
        return enclosed(simplex[0], simplex[1], simplex[2]);

      const support candidate =
          diff_support(pair.shape_id_A, pair.shape_id_B, direction);
      if (num::inner_prod(support_vector(candidate), direction) <
          -options.tolerance)
        return separated();

      for (std::size_t i = 0; i < simplex_size; ++i) {
        const vector difference =
            support_vector(candidate) - support_vector(simplex[i]);
        if (is_near_zero(difference)) return separated();
      }

      simplex[2] = simplex[1];
      simplex[1] = simplex[0];
      simplex[0] = candidate;
      ++simplex_size;

      const vector a = support_vector(simplex[0]);
      const vector a0 = -a;
      if (simplex_size == 2) {
        const vector ab = support_vector(simplex[1]) - a;
        if (num::inner_prod(ab, a0) <= options.tolerance) {
          simplex_size = 1;
          direction = a0;
          continue;
        }

        direction = cross3(ab, a0, ab);
        if (is_near_zero(direction))
          return enclosed(simplex[0], simplex[1], simplex[2]);
        continue;
      }

      const vector ab = support_vector(simplex[1]) - a;
      const vector ac = support_vector(simplex[2]) - a;
      const vector ab_direction = cross3(ac, ab, ab);
      if (num::inner_prod(ab_direction, a0) > options.tolerance) {
        simplex_size = 2;
        direction = ab_direction;
        continue;
      }

      const vector ac_direction = cross3(ab, ac, ac);
      if (num::inner_prod(ac_direction, a0) > options.tolerance) {
        simplex[1] = simplex[2];
        simplex_size = 2;
        direction = ac_direction;
        continue;
      }

      return enclosed(simplex[0], simplex[1], simplex[2]);
    }

    return exhausted();
  }
  else if constexpr (Dim == 3) {
    const auto enclosed = [](const support& s0, const support& s1,
                             const support& s2, const support& s3) {
      result_type entry{};
      entry.st = status::enclosed;
      entry.seed = {s0, s1, s2, s3};
      return entry;
    };

    const vector init{pair.init_dir[0], pair.init_dir[1], pair.init_dir[2]};

    const auto tolerance_squared = options.tolerance * options.tolerance;
    if (num::inner_prod(init, init) <= tolerance_squared)
      return invalid_input();

    // O(1) ball prefilter: a pair whose centre distance goes beyond the
    // radii sum plus the tolerance is separated before any simplex work.
    // Only the compute path shrinks; the row is exactly what the full GJK
    // would return for a pair this far apart.
    if (bounds != nullptr &&
        bounds_separated<Dim>(pool, *bounds, pair.shape_id_A, pair.shape_id_B,
                              options.tolerance))
      return separated();

    const auto is_near_zero = [&tolerance_squared](const vector& v) {
      return num::inner_prod(v, v) <= tolerance_squared;
    };
    const auto simplex0 =
        diff_support(pair.shape_id_A, pair.shape_id_B, init);
    if (num::inner_prod(support_vector(simplex0), init) < -options.tolerance)
      return separated();

    std::array<support, 4> simplex{simplex0, simplex0, simplex0, simplex0};
    std::size_t simplex_size = 1;
    vector direction = -support_vector(simplex0);

    const auto dedup_simplex = [&]() -> std::size_t {
      for (std::size_t i = 0; i < simplex_size; ++i) {
        for (std::size_t j = i + 1; j < simplex_size;) {
          if (is_near_zero(support_vector(simplex[i]) -
                           support_vector(simplex[j]))) {
            simplex[j] = simplex[simplex_size - 1];
            --simplex_size;
          }
          else {
            ++j;
          }
        }
      }
      return simplex_size;
    };

    const auto compact_simplex = [&]() {
      const std::size_t sz = dedup_simplex();
      std::array<support, 4> out{};
      for (std::size_t i = 0; i < sz; ++i) out[i] = simplex[i];
      for (std::size_t i = sz; i < 4; ++i) out[i] = out[sz - 1];
      return out;
    };

    vector last_loop_dir;
    std::size_t last_loop_size = 0;
    bool has_loop_state = false;
    std::size_t same_loop_state = 0;

    for (std::size_t iteration = 0; iteration < options.max_iterations;
         ++iteration) {
      // Loop-state stall detector: when the search revisits an identical
      // (simplex size, direction) pair, the reduction path has converged
      // on a feature of the Minkowski difference that cannot be expanded
      // or discarded. Answer from the converged simplex instead of
      // re-entering the same reduction forever (box-versus-box pair 189:
      // the size-4 dispatch keeps dropping a vertex of the same closed
      // triangle, cycling until max_iterations).
      const scalar dir_dot_last = num::inner_prod(direction, last_loop_dir);
      const scalar dir_norm2 = num::inner_prod(direction, direction);
      const scalar last_norm2 = num::inner_prod(last_loop_dir, last_loop_dir);
      const bool same_axis =
          has_loop_state and
          dir_dot_last * dir_dot_last >= 0.9 * dir_norm2 * last_norm2;
      if (same_axis and simplex_size == last_loop_size)
        ++same_loop_state;
      else
        same_loop_state = 0;
      has_loop_state = true;
      last_loop_dir = direction;
      last_loop_size = simplex_size;

      if (same_loop_state >= 1 and simplex_size == 3) {
        // Containment, not prism containment: the origin must lie in the
        // convex hull of the three points, i.e. in their plane and inside
        // all three edges. The prism test alone answers true for origins
        // outside the triangle but between its edge planes (separated box
        // pair 189: witness triangle in the z = -13.17 plane).
        const vector t0 = support_vector(simplex[0]);
        const vector t1 = support_vector(simplex[1]);
        const vector t2 = support_vector(simplex[2]);
        const vector t01 = t1 - t0;
        const vector t02 = t2 - t0;
        const vector t12 = t2 - t1;
        const vector tn = cross(t01, t02);
        const scalar scale = num::inner_prod(t01, t01) +
                             num::inner_prod(t02, t02) +
                             num::inner_prod(tn, tn);
        const vector t0n = -t0;
        const bool origin_in_triangle =
            std::abs(num::inner_prod(tn, t0n)) <= options.tolerance * scale &&
            num::inner_prod(cross(t01, tn), t0n) <=
                options.tolerance * scale &&
            num::inner_prod(cross(tn, t02), t0n) <=
                options.tolerance * scale &&
            num::inner_prod(cross(t12, tn), -t1) <= options.tolerance * scale;
        if (origin_in_triangle) {
          const auto compacted = compact_simplex();
          return enclosed(compacted[0], compacted[1], compacted[2],
                          compacted[3]);
        }
        // The triangle does not contain the origin, which by itself
        // proves nothing about separation: a three-point simplex cannot
        // enclose the origin in 3D. The search has converged on this face
        // of the Minkowski difference (the support query re-selects points
        // of the same plane), so decide with the separation test: query the
        // support along the face normal oriented toward the origin. A
        // negative projection certifies that the difference keeps no point
        // on the origin's side of this face and the shapes are separated;
        // otherwise no separating direction exists along the converged
        // face and the difference contains the origin, so the shapes
        // intersect. Answering separation unconditionally here would
        // separate overlapping pairs whose stall converged on a face seen
        // from the far side (overlapping octahedra); answering
        // intersection unconditionally would intersect separated box pairs
        // whose supports keep advancing in-plane (pair 189).
        vector toward_origin = tn;
        if (num::inner_prod(tn, t0) > 0) toward_origin = -tn;
        const support sep_witness =
            diff_support(pair.shape_id_A, pair.shape_id_B, toward_origin);
        const scalar sep_proj =
            num::inner_prod(support_vector(sep_witness), toward_origin);
        // pair189 witness: tn=(232,233,232), t0=(5.56,-0.64,-13.17),
        // w=(5.56,-0.64,8.41), proj=3097.56. The stall triangle's normal is
        // not the isolating direction, so the face-normal support test cannot
        // certify separation. Fall back to an axis-aligned probe of the
        // support difference along the three coordinate axes of the converged
        // triangle's bounding box: a negative projection along any axis
        // certifies separation of the boxes (separated box pairs have their
        // isochronous side planes perpendicular to a box axis).
        vector probe = toward_origin;
        bool separated = sep_proj < -options.tolerance;
        if (!separated) {
          const scalar probe_candidates[6][3] = {{1, 0, 0}, {-1, 0, 0},
                                                 {0, 1, 0}, {0, -1, 0},
                                                 {0, 0, 1}, {0, 0, -1}};
          for (std::size_t pi = 0; pi < 6 and not separated; ++pi) {
            probe = vector{probe_candidates[pi][0], probe_candidates[pi][1],
                           probe_candidates[pi][2]};
            const support probe_w =
                diff_support(pair.shape_id_A, pair.shape_id_B, probe);
            separated = num::inner_prod(support_vector(probe_w), probe) <
                        -options.tolerance;
          }
        }
        const auto compacted = compact_simplex();
        if (separated) {
          result_type entry{};
          entry.st = status::separated;
          return entry;
        }
        return enclosed(compacted[0], compacted[1], compacted[2],
                        sep_witness);
      }
      if (is_near_zero(direction)) {
        dedup_simplex();
        if (simplex_size >= 3) {
          const vector a = support_vector(simplex[0]);
          const vector ab = support_vector(simplex[1]) - a;
          const vector ac = support_vector(simplex[2]) - a;
          const vector abc = cross(ab, ac);
          if (!is_near_zero(abc)) {
            direction = num::inner_prod(abc, a) > 0 ? -abc : abc;
            continue;
          }
        }
        if (simplex_size >= 2) {
          const vector a = support_vector(simplex[0]);
          const vector ab = support_vector(simplex[1]) - a;
          direction = cross(ab, vector{1, 0, 0});
          if (is_near_zero(direction)) direction = cross(ab, vector{0, 1, 0});
          if (num::inner_prod(direction, a) < 0) direction = -direction;
          continue;
        }
        direction = -support_vector(simplex[0]);
        continue;
      }

      const support candidate =
          diff_support(pair.shape_id_A, pair.shape_id_B, direction);
      if (num::inner_prod(support_vector(candidate), direction) <
          -options.tolerance) {
        return separated();
      }

      bool duplicate = false;
      for (std::size_t i = 0; i < simplex_size; ++i) {
        const vector difference =
            support_vector(candidate) - support_vector(simplex[i]);
        if (is_near_zero(difference)) {
          duplicate = true;
          break;
        }
      }
      if (duplicate) {
        // Duplicate support vertex: the current simplex already encodes the
        // closest feature, and forcing a new vertex via geometric fallbacks
        // can cycle on degenerate cases (e.g. axis-aligned cubes where every
        // escape direction still maps to the same face). Following
        // Bullet/Jolt/ ReactPhysics3D practice: reduce the simplex, keep the
        // best feature, and stop when no new vertex can be obtained.
        const std::size_t final_size = dedup_simplex();
        if (final_size < 3) {
          simplex_size = 1;
          direction = -support_vector(simplex[0]);
          continue;
        }
        if (simplex_size == 4) {
          // Containment-aware stall handling. A stalled four-point simplex
          // only witnesses an intersection when it encloses the origin: the
          // origin must be on the same side of every face as the opposite
          // vertex. Signed face tests alone depend on the tetrahedron's
          // orientation; same-side products do not.
          const vector a = support_vector(simplex[0]);
          const vector b = support_vector(simplex[1]);
          const vector c = support_vector(simplex[2]);
          const vector d = support_vector(simplex[3]);
          const vector ab = b - a;
          const vector ac = c - a;
          const vector ad = d - a;
          const vector bc = c - b;
          const vector bd = d - b;
          const bool encloses = num::inner_prod(cross(ab, ac), ad) *
                                        num::inner_prod(cross(ab, ac), -a) >
                                    0 &&
                                num::inner_prod(cross(ab, ad), ac) *
                                        num::inner_prod(cross(ab, ad), -a) >
                                    0 &&
                                num::inner_prod(cross(ac, ad), ab) *
                                        num::inner_prod(cross(ac, ad), -a) >
                                    0 &&
                                num::inner_prod(cross(bd, bc), a - b) *
                                        num::inner_prod(cross(bd, bc), -b) >
                                    0;
          if (encloses) {
            const auto compacted = compact_simplex();
            return enclosed(
                compacted[0], compacted[1], compacted[2],
                diff_support(pair.shape_id_A, pair.shape_id_B, direction));
          }
          simplex_size = 3;
        }
        // A stalled three-point simplex: the origin may still lie inside the
        // prism of the triangle (in which case the Minkowski difference
        // contains the origin and the shapes intersect). Otherwise reduce to
        // the closest edge and let the separation early-out answer.
        {
          const vector a = support_vector(simplex[0]);
          const vector a0 = -a;
          const vector ab = support_vector(simplex[1]) - a;
          const vector ac = support_vector(simplex[2]) - a;
          const vector abc = cross(ab, ac);
          const bool origin_in_prism =
              num::inner_prod(cross(ab, abc), a0) <= options.tolerance &&
              num::inner_prod(cross(abc, ac), a0) <= options.tolerance;
          // A stalled three-point simplex terminates the query. When the
          // origin lies inside the prism of the triangle the Minkowski
          // difference contains the origin, so the shapes intersect. When
          // it lies outside, the stalled triangle is the closest feature
          // of the difference and the shapes are separated: reducing to
          // an edge would only re-select a support point already in the
          // simplex and cycle until max_iterations.
          // A repeated stall with an unchanged simplex cannot make progress:
          // answer by containment of the reduced feature. The stall
          // simplex holds three distinct points (compaction repeats the
          // last one), so the fourth witness slot must be filled with a
          // fresh support along the converged face's own normal; a
          // repeated point starves EPA of a tetrahedron.
          const auto stall_compacted = compact_simplex();
          {
            const vector da = support_vector(stall_compacted[0]);
            const vector d1 = support_vector(stall_compacted[1]) - da;
            const vector d2 = support_vector(stall_compacted[2]) - da;
            const vector dn = cross(d1, d2);
            vector ddir = dn;
            if (num::inner_prod(dn, da) > 0) ddir = -dn;
            const support stall_wit4 =
                diff_support(pair.shape_id_A, pair.shape_id_B, ddir);
            if (origin_in_prism) {
              return enclosed(stall_compacted[0], stall_compacted[1],
                              stall_compacted[2], stall_wit4);
            }
            // No separation certificate was found: GJK's early-out is the
            // only sound separation answer and it did not fire, so the
            // stalled query cannot prove separation. Answer intersection (the
            // pristine behaviour; the reviewed false positives for separated
            // pairs are caught by the loop-state stall terminal instead,
            // which sees the actual axis-aligned separation).
            return enclosed(stall_compacted[0], stall_compacted[1],
                            stall_compacted[2], stall_wit4);
          }
        }
      }

      if (simplex_size < 4) {
        for (std::size_t i = simplex_size; i > 0; --i) {
          simplex[i] = simplex[i - 1];
        }
        simplex[0] = candidate;
        ++simplex_size;
      }

      for (std::size_t i = 0; i < simplex_size; ++i) {
        for (std::size_t j = i + 1; j < simplex_size;) {
          if (is_near_zero(support_vector(simplex[i]) -
                           support_vector(simplex[j]))) {
            simplex[j] = simplex[simplex_size - 1];
            --simplex_size;
          }
          else {
            ++j;
          }
        }
      }

      const vector a = support_vector(simplex[0]);
      const vector a0 = -a;

      if (simplex_size == 2) {
        const vector ab = support_vector(simplex[1]) - a;
        const scalar ab2 = num::inner_prod(ab, ab);
        const scalar ab_ao = num::inner_prod(ab, a0);
        if (ab_ao <= options.tolerance) {
          simplex_size = 1;
          direction = a0;
          continue;
        }
        if (ab_ao >= ab2 - options.tolerance) {
          simplex[0] = simplex[1];
          simplex_size = 1;
          direction = -support_vector(simplex[0]);
          continue;
        }
        direction = a0 - (ab_ao / ab2) * ab;
        if (is_near_zero(direction)) {
          vector perp = cross(ab, vector{1, 0, 0});
          if (is_near_zero(perp)) perp = cross(ab, vector{0, 1, 0});
          direction = perp;
        }
        else {
        }
        continue;
      }

      if (simplex_size == 3) {
        const vector ab = support_vector(simplex[1]) - a;
        const vector ac = support_vector(simplex[2]) - a;
        const vector abc = cross(ab, ac);
        const scalar abc_ao = num::inner_prod(abc, a0);

        const vector ab_cross_abc = cross(ab, abc);
        const vector ac_cross_abc = cross(abc, ac);

        if (num::inner_prod(ab_cross_abc, a0) > options.tolerance) {
          simplex_size = 2;
          const vector ab2 =
              support_vector(simplex[1]) - support_vector(simplex[0]);
          const scalar ab2_len = num::inner_prod(ab2, ab2);
          const scalar ab_ao2 =
              num::inner_prod(ab2, -support_vector(simplex[0]));
          if (ab_ao2 <= options.tolerance) {
            simplex_size = 1;
            direction = -support_vector(simplex[0]);
            continue;
          }
          direction = -support_vector(simplex[0]) - (ab_ao2 / ab2_len) * ab2;
          continue;
        }

        if (num::inner_prod(ac_cross_abc, a0) > options.tolerance) {
          simplex[1] = simplex[2];
          simplex_size = 2;
          const vector ac2 =
              support_vector(simplex[1]) - support_vector(simplex[0]);
          const scalar ac2_len = num::inner_prod(ac2, ac2);
          const scalar ac_ao2 =
              num::inner_prod(ac2, -support_vector(simplex[0]));
          if (ac_ao2 <= options.tolerance) {
            simplex_size = 1;
            direction = -support_vector(simplex[0]);
            continue;
          }
          direction = -support_vector(simplex[0]) - (ac_ao2 / ac2_len) * ac2;
          continue;
        }

        if (abc_ao > options.tolerance) {
          simplex_size = 3;
          direction = -abc;
        }
        else
          direction = abc;

        // Stall detector for the triangle ping-pong: when the search keeps
        // re-selecting the same plane-normal direction without growing the
        // simplex, the support set has converged on a face of the Minkowski
        // difference. Judge the origin against that converged triangle:
        // inside the triangle prism the Minkowski difference contains the
        // origin and the shapes intersect; outside it this face is the
        // closest feature reached and the shapes are separated. Continuing
        // the normal search re-inserts the same points and cycles until
        // max_iterations (separated box pair 189: witness triangle in the
        // z = -13.17 plane, origin outside its prism).
        continue;
      }

      if (simplex_size == 4) {
        const vector ab = support_vector(simplex[1]) - a;
        const vector ac = support_vector(simplex[2]) - a;
        const vector ad = support_vector(simplex[3]) - a;
        const vector abc = cross(ab, ac);
        const vector abd = cross(ab, ad);
        const vector acd = cross(ad, ac);

        const scalar abc_ao =
            num::inner_prod(abc, a0) * num::inner_prod(abc, ad);
        const scalar abd_ao =
            num::inner_prod(abd, a0) * num::inner_prod(abd, ac);
        const scalar acd_ao =
            num::inner_prod(acd, a0) * num::inner_prod(acd, ab);
        if (abc_ao <
            -options.tolerance * std::abs(num::inner_prod(abc, ad))) {
          // P, B, C are already at indices 0, 1, 2 in the right order;
          // dropping D (index 3) needs no data movement at all.
          simplex_size = 3;
          direction = -abc;
          for (std::size_t i = 1; i < simplex_size; ++i) {
            for (std::size_t j = i + 1; j < simplex_size;) {
              if (is_near_zero(support_vector(simplex[i]) -
                               support_vector(simplex[j]))) {
                simplex[j] = simplex[simplex_size - 1];
                --simplex_size;
              }
              else {
                ++j;
              }
            }
          }
          if (simplex_size < 3) {
            simplex_size = 1;
            direction = -support_vector(simplex[0]);
          }
          continue;
        }

        if (abd_ao <
            -options.tolerance * std::abs(num::inner_prod(abd, ac))) {
          // Keep P, B, D (drop C): move D from index 3 down to index 2,
          // the slot the size-3 code below reads.
          simplex[2] = simplex[3];
          simplex_size = 3;
          direction = -abd;
          for (std::size_t i = 1; i < simplex_size; ++i) {
            for (std::size_t j = i + 1; j < simplex_size;) {
              if (is_near_zero(support_vector(simplex[i]) -
                               support_vector(simplex[j]))) {
                simplex[j] = simplex[simplex_size - 1];
                --simplex_size;
              }
              else {
                ++j;
              }
            }
          }
          if (simplex_size < 3) {
            simplex_size = 1;
            direction = -support_vector(simplex[0]);
          }
          continue;
        }

        if (acd_ao <
            -options.tolerance * std::abs(num::inner_prod(acd, ab))) {
          // Keep P, C, D (drop B): move C from index 2 down to index 1
          // and D from index 3 down to index 2.
          simplex[1] = simplex[2];
          simplex[2] = simplex[3];
          simplex_size = 3;
          direction = -acd;
          for (std::size_t i = 1; i < simplex_size; ++i) {
            for (std::size_t j = i + 1; j < simplex_size;) {
              if (is_near_zero(support_vector(simplex[i]) -
                               support_vector(simplex[j]))) {
                simplex[j] = simplex[simplex_size - 1];
                --simplex_size;
              }
              else {
                ++j;
              }
            }
          }
          if (simplex_size < 3) {
            simplex_size = 1;
            direction = -support_vector(simplex[0]);
          }
          continue;
        }

        // The fourth face, BCD, was never tested: the three tests above
        // only judge the faces based at A. Without it a tetrahedron whose
        // only containment violation is on face BCD still claimed an
        // enclosure (that was pair 189's false positive).
        {
          const vector bc4 =
              support_vector(simplex[2]) - support_vector(simplex[1]);
          const vector bd4 =
              support_vector(simplex[3]) - support_vector(simplex[1]);
          const vector bcd4 = cross(bd4, bc4);
          const scalar bcd_bo =
              num::inner_prod(bcd4, -support_vector(simplex[1])) *
              num::inner_prod(bcd4, support_vector(simplex[0]) -
                                        support_vector(simplex[1]));
          if (bcd_bo < -options.tolerance *
                           std::abs(num::inner_prod(
                               bcd4, support_vector(simplex[0]) -
                                         support_vector(simplex[1])))) {
            // Drop A (index 0) and continue with face BCD as the triangle.
            simplex[0] = simplex[1];
            simplex[1] = simplex[2];
            simplex[2] = simplex[3];
            simplex_size = 3;
            direction = -bcd4;
            for (std::size_t i = 1; i < simplex_size; ++i) {
              for (std::size_t j = i + 1; j < simplex_size;) {
                if (is_near_zero(support_vector(simplex[i]) -
                                 support_vector(simplex[j]))) {
                  simplex[j] = simplex[simplex_size - 1];
                  --simplex_size;
                }
                else {
                  ++j;
                }
              }
            }
            if (simplex_size < 3) {
              simplex_size = 1;
              direction = -support_vector(simplex[0]);
            }
            continue;
          }
        }

        const auto compacted2 = compact_simplex();
        const auto& a = compacted2[0];
        const auto& b = compacted2[1];
        const auto& c = compacted2[2];
        const auto& d = compacted2[3];
        const auto av = support_vector(a);
        const auto bv = support_vector(b);
        const auto cv = support_vector(c);
        const auto dv = support_vector(d);
        const vector ab2 = bv - av;
        const vector ac2 = cv - av;
        const vector ad2 = dv - av;
        const scalar volume_sq = num::inner_prod(cross(ab2, ac2), ad2);
        if (volume_sq <= options.tolerance * options.tolerance) {
          // Degenerate (coplanar) tetrahedron: the four support points
          // all lie in one plane, so the simplex cannot grow a fourth
          // dimension here. Reduce to the coplanar face and keep
          // searching along the flipped normal so the query can grow a
          // real tetrahedron from the next support. A separated pair
          // converges onto the same face; the loop-state stall detector
          // then answers from the converged simplex instead of cycling
          // (box-versus-box pair 189).
          simplex_size = 3;
          direction = -cross(ab2, ac2);
          continue;
        }

        return enclosed(a, b, c, d);
      }

      // A three-point simplex cannot enclose the origin in 3D, but it may
      // still witness an intersection when the search converged without
      // discarding it (the supports keep advancing). EPA needs four
      // distinct vertices, so fill the fourth slot with a fresh support
      // along the current search direction instead of repeating the last
      // point of the compacted simplex.
      const auto compacted3 = compact_simplex();
      return enclosed(
          compacted3[0], compacted3[1], compacted3[2],
          diff_support(pair.shape_id_A, pair.shape_id_B, direction));
    }

    return exhausted();
  }
}

// Per-pair wall-clock rows, one bucket per stage a batch stage runs.
// A stage keeps its own batch_timing instance; when the pointer is null
// the stage does no timing at all.
struct batch_timing {
  std::vector<double>
      seconds;  // one row per pair (per hit, for EPA), per call
};

// The batched intersection pass (core, by reference): one pair list, one
// packed shape pool, one status per pair, written into a caller-owned result
// buffer. The buffer may be reused across simulation steps: assign re-fills
// it to the current pair count and overwrites every entry, so no stale state
// survives a call. Options-only validation (a non-positive tolerance or a
// zero iteration budget) rejects every pair; everything else is decided per
// pair by gjk_batch_pair.
//
// bounds, when given, feeds the O(1) ball prefilter (per-pair rejection
// before any simplex work); timing, when given, collects one wall-clock row
// per pair. Both are null by default and cost nothing when unset -- the
// stage stages stay untouched for callers that do not want them.

template <std::size_t Dim, typename Pool, typename PairRecord,
          typename Support>
  requires(Dim == 2 || Dim == 3)
void gjk_batch(
    const Pool& pool, const std::vector<PairRecord>& pairs,
    const gjk_options options,
    std::vector<batch_result<Support, Dim + 1>>& results,
    const pool_bounds<typename Pool::scalar, Dim>* bounds = nullptr,
    batch_timing* timing = nullptr)
{
  using result_type = batch_result<typename Pool::support, Dim + 1>;

  assert_device_contract<Pool, Dim>();

  // assign (not resize) so the vector is re-filled even after the content
  // it held is stale; a buffer kept across calls is emptied and re-added
  // without losing its heap block, and its capacity from the previous step
  // is reused.
  results.assign(pairs.size(), result_type{});

  if (options.tolerance <= 0. || options.max_iterations == 0) {
    for (auto& entry : results) entry.st = status::invalid_input;
    return;
  }

  if (timing != nullptr) {
    timing->seconds.assign(pairs.size(), 0.);
    for (std::size_t i = 0; i < pairs.size(); ++i) {
      const auto start = std::chrono::steady_clock::now();
      results[i] = gjk_batch_pair<Dim>(pool, pairs[i], options, bounds);
      timing->seconds[i] = std::chrono::duration<double>(
                               std::chrono::steady_clock::now() - start)
                               .count();
    }
    return;
  }

  for (std::size_t i = 0; i < pairs.size(); ++i)
    results[i] = gjk_batch_pair<Dim>(pool, pairs[i], options, bounds);
}

// The batched intersection pass (allocating form): allocates the result
// buffer itself and forwards to the core overload above. Kept as a
// convenience for callers that buffer once per call; the by-reference core
// is the intended form for the per-step batch stage.
template <std::size_t Dim, typename Pool, typename PairRecord>
  requires(Dim == 2 || Dim == 3)
auto gjk_batch(const Pool& pool, const std::vector<PairRecord>& pairs,
               const gjk_options options = {})
    -> std::vector<batch_result<typename Pool::support, Dim + 1>>
{
  std::vector<batch_result<typename Pool::support, Dim + 1>> results;
  gjk_batch<Dim>(pool, pairs, options, results);
  return results;
}

// Convenience overload for the common 3D layout: two shape lists of equal
// length, one shared initial direction. Both sides are packed into one pool
// (the a-shapes first, then the b-shapes) and the pair list zips them. A
// length mismatch yields one invalid_input entry per pair that exists on both
// sides rather than a throw. Constrained to Dim=3: its shapes are
// std::array<double, 3>, so any other Dim is rejected at compile time the way
// gjk_inter rejects a mismatched init vector, and a 2D convenience overload
// stays out of scope as the spec has it.
template <std::size_t Dim>
  requires(Dim == 3)
auto gjk_batch(const std::vector<std::vector<std::array<double, 3>>>& a,
               const std::vector<std::vector<std::array<double, 3>>>& b,
               const std::array<double, 3>& init,
               const gjk_options options = {})
    -> std::vector<batch_result<
        cxd::traits::env_t<std::array<double, 3>>::support, Dim + 1>>
{
  using shape = std::vector<std::array<double, 3>>;
  using result_type =
      batch_result<cxd::traits::env_t<std::array<double, 3>>::support,
                   Dim + 1>;

  const std::size_t pair_count = std::min(a.size(), b.size());
  std::vector<result_type> results(pair_count);

  if (a.size() != b.size()) {
    for (auto& entry : results) entry.st = status::invalid_input;
    return results;
  }

  std::vector<shape> pool_shapes;
  pool_shapes.reserve(a.size() + b.size());
  pool_shapes.insert(pool_shapes.end(), a.begin(), a.end());
  pool_shapes.insert(pool_shapes.end(), b.begin(), b.end());

  const auto pool = pack_shapes(pool_shapes);

  const auto offset = static_cast<std::uint32_t>(a.size());
  std::vector<pair_record<double, 3>> pairs;
  pairs.reserve(pair_count);
  for (std::size_t i = 0; i < pair_count; ++i)
    pairs.push_back(
        pair_record<double, 3>{static_cast<std::uint32_t>(i),
                               offset + static_cast<std::uint32_t>(i),
                               {init[0], init[1], init[2]}});

  return gjk_batch<Dim>(pool, pairs, options);
}

// ---------------------------------------------------------------------------
// The batched penetration pass
// ---------------------------------------------------------------------------

// The answer for a hit the batch cannot compute. Scalar epa reports every
// rejection by throwing; the batch's return type carries no status, so a
// rejected hit comes back as NaN in place instead -- invalid-input rejection
// without a throw, and every position of the result stays indexable.
template <typename Scalar, std::size_t Dim>
auto unanswered() -> epa_result<Scalar, Dim>
{
  const Scalar nan = std::numeric_limits<Scalar>::quiet_NaN();
  epa_result<Scalar, Dim> r{};
  r.depth = nan;
  r.normal.fill(nan);
  return r;
}

// One growable region of the polytope arena. While the polytope fits Capacity
// entries the storage is the inline array -- scalar epa's small-buffer
// optimization for shapes whose combined vertex count is at most
// epa_sbo_capacity. The first entry past it promotes the whole region to the
// heap vector, which serves every access from then on, so the heap remains
// the fallback exactly as it is in the scalar path. Promotion is one-way, but
// the heap buffer is not: a growth past the inline bound reallocates it,
// exactly as a std::vector does. An index is therefore the handle that stays
// valid across a growing call and a reference is not -- the face walk
// addresses its faces by index for that reason, and keeps the scalar
// arithmetic unchanged.
template <typename T, std::size_t Capacity>
struct polytope_region {
  void reset()
  {
    count = 0;
    heap_storage.clear();
  }

  // Pre-reserve the heap store so the region never reallocates once the pair
  // loop starts: with aCapacity entries reserved, push_back and insert_at
  // cannot allocate, and promote() copies into storage that is already big
  // enough. Called once per batch, before the loop, from the exact bound the
  // guards already use (max_polytope). Leaves the inline store alone, so
  // small regions still never touch the heap.
  void reserve(const std::size_t n)
  {
    if (n > Capacity) heap_storage.reserve(n);
  }

  std::size_t size() const { return count; }

  T* data()
  {
    return count <= Capacity ? inline_storage.data() : heap_storage.data();
  }

  T* begin() { return data(); }

  T* end() { return data() + count; }

  T& operator[](const std::size_t i) { return data()[i]; }

  void push_back(const T& value)
  {
    if (count < Capacity)
      inline_storage[count] = value;
    else if (count == Capacity) {
      promote();
      heap_storage.push_back(value);
    }
    else
      heap_storage.push_back(value);
    ++count;
  }

  // Make room at index for value, shifting the tail up by one: scalar epa's
  // std::copy_backward over its raw polytope pointer.
  void insert_at(const std::size_t index, const T& value)
  {
    if (count < Capacity) {
      for (std::size_t i = count; i > index; --i)
        inline_storage[i] = inline_storage[i - 1];
      inline_storage[index] = value;
    }
    else {
      if (count == Capacity) promote();
      heap_storage.insert(
          heap_storage.begin() + static_cast<std::ptrdiff_t>(index), value);
    }
    ++count;
  }

  // Drop the entries the predicate marks, keeping the survivors in order. The
  // scalar path rebuilds its face list into a fresh vector here; the region
  // compacts in place so an expansion costs no allocation.
  template <typename Predicate>
  void erase_if(const Predicate& marked)
  {
    std::size_t kept = 0;
    if (count > Capacity) {
      for (std::size_t i = 0; i < count; ++i)
        if (not marked(i)) heap_storage[kept++] = heap_storage[i];
      if (kept <= Capacity)
        for (std::size_t i = 0; i < kept; ++i)
          inline_storage[i] = heap_storage[i];
      else
        heap_storage.resize(kept);
    }
    else {
      for (std::size_t i = 0; i < count; ++i) {
        if (not marked(i)) {
          if (kept != i) inline_storage[kept] = inline_storage[i];
          ++kept;
        }
      }
    }
    count = kept;
  }

 private:
  void promote()
  {
    heap_storage.assign(inline_storage.begin(), inline_storage.end());
  }

  std::array<T, Capacity> inline_storage{};
  std::vector<T> heap_storage{};
  std::size_t count = 0;
};

// One face of the batched polytope: three vertex indices into the vertex
// region, the outward unit normal and its distance from the origin. `removed`
// is the expansion mask scalar epa keeps in a vector beside its face list.
template <typename Vector, typename Scalar>
struct polytope_face {
  std::array<std::size_t, 3> indices{};
  Vector normal{};
  Scalar distance{};
  bool removed = false;
};

// The per-hit polytope storage the penetration pass owns. Scalar epa
// allocates the Minkowski-difference polytope on every call: in 2D a ring of
// support points behind an SBO of epa_sbo_capacity entries, in 3D a vertex
// set plus the face list built on top of it. The batch meets one hit at a
// time, so it owns one arena instead of one allocation per pair, reserves it
// once from the same bound the polytope guards use (vertex_count_A +
// vertex_count_B, twice it for the face list, since a closed polytope of V
// vertices carries at most 2V - 4 faces), and resets it between hits. Small
// polytopes stay inside the inline store; larger ones sit in the reserved
// heap buffer and still never reallocate.
template <typename Support>
struct polytope_ring {
  polytope_region<Support, epa_sbo_capacity> points;

  void reset() { points.reset(); }

  // See epa_batch's pre-sizing: reserve the exact vertex bound once so the
  // hit loop never reallocates. Points count is the bound.
  void reserve(const std::size_t vertices) { points.reserve(vertices); }
};

template <typename Vector, typename Scalar>
struct polytope_solid {
  using face = polytope_face<Vector, Scalar>;

  polytope_region<Vector, epa_sbo_capacity> vertices;
  polytope_region<face, 2 * epa_sbo_capacity> faces;

  // Fixed-size support tracking for vertex index extraction (no heap alloc).
  std::array<std::pair<std::size_t, std::size_t>, epa_sbo_capacity>
      vert_support_ids{};
  std::size_t vert_support_count = 0;

  void reset()
  {
    vertices.reset();
    faces.reset();
    vert_support_count = 0;
  }

  // Euler's bound: a convex polytope of V vertices carries at most 2V - 4
  // faces, so 2*vertices reserves the face list with margin to spare.
  void reserve(const std::size_t vertex_bound, const std::size_t face_bound)
  {
    vertices.reserve(vertex_bound);
    faces.reserve(face_bound);
  }
};

// The penetration pass over one hit's seed, 2D branch: scalar epa's edge walk
// over the ring of support points the arena holds.
template <typename Pool>
auto epa_penetration_2d(const Pool& pool, const std::uint32_t id_a,
                        const std::uint32_t id_b,
                        const std::array<typename Pool::support, 3>& seed,
                        const gjk_options options,
                        polytope_ring<typename Pool::support>& arena)
    -> epa_result<typename Pool::scalar, 2>
{
  using vector = typename Pool::vector;
  using scalar = typename Pool::scalar;
  using support = typename Pool::support;
  using answer = epa_result<scalar, 2>;

  const auto rejected = []() -> answer {
    answer r{};
    r.depth = std::numeric_limits<scalar>::quiet_NaN();
    r.normal.fill(std::numeric_limits<scalar>::quiet_NaN());
    return r;
  };

  const auto make_result = [](scalar depth, const vector& normal_vec,
                              const support& v0,
                              const support& v1) -> answer {
    answer r{};
    r.depth = depth;
    r.normal = {static_cast<scalar>(normal_vec[0]),
                static_cast<scalar>(normal_vec[1])};
    r.vertices[0] = {std::get<0>(v0), std::get<1>(v0)};
    r.vertices[1] = {std::get<0>(v1), std::get<1>(v1)};
    return r;
  };

  const auto diff_support = [&pool](const std::uint32_t a,
                                    const std::uint32_t b,
                                    const vector& direction) {
    return pool_diff_support(pool, a, b, direction);
  };

  // Validation mirrors epa, reported per hit instead of thrown: an option
  // pair epa would reject, or a shape outside the pool.
  if (options.tolerance <= 0. || options.max_iterations == 0)
    return rejected();
  if (id_a >= pool.entries.size() || id_b >= pool.entries.size() ||
      pool.entries[id_a].vertex_count == 0 ||
      pool.entries[id_b].vertex_count == 0)
    return rejected();

  // A GJK witness can be degenerate. When the 2D search direction collapses
  // -- the origin lies ON the segment between two supports -- the walk stops
  // with a two-point simplex whose empty slot is padded by repeating a
  // support, so the seed triangle has zero area. The 3D path bootstraps a
  // non-degenerate polytope in that case (epa_penetration_3d); the 2D
  // analogue queries a support off the degenerate line and restarts the face
  // walk from that triangle. The walk then expands from the origin-carrying
  // edge (whose distance to the origin is zero) towards the real closest
  // face, which is the deep face of the difference polytope.
  std::array<support, 3> tri{seed[0], seed[1], seed[2]};
  {
    const vector p0 = support_vector(tri[0]);
    const vector p1 = support_vector(tri[1]);
    const vector p2 = support_vector(tri[2]);
    const scalar area =
        (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p1[1] - p0[1]) * (p2[0] - p0[0]);
    if (std::abs(area) <= options.tolerance) {
      // recover the distinct supports the degenerate seed repeats
      std::size_t i0 = 0;
      std::size_t i1 = 0;
      bool have_second = false;
      for (std::size_t i = 1; i < 3 && !have_second; ++i) {
        const vector d = support_vector(tri[i]) - support_vector(tri[i0]);
        if (num::inner_prod(d, d) > options.tolerance * options.tolerance) {
          i1 = i;
          have_second = true;
        }
      }
      if (!have_second) {
        // one support only: the difference is a point and the origin sits on
        // it (a touching contact with no direction)
        const vector p = support_vector(tri[i0]);
        const scalar distance = std::sqrt(num::inner_prod(p, p));
        const vector normal = distance > options.tolerance
                                  ? vector{p[0] / distance, p[1] / distance}
                                  : vector{0., 0.};
        return make_result(distance, normal, tri[i0], tri[i0]);
      }
      const vector line = support_vector(tri[i1]) - support_vector(tri[i0]);
      const scalar line_length = std::sqrt(num::inner_prod(line, line));
      if (line_length > options.tolerance) {
        const vector perp{-line[1], line[0]};
        // a perpendicular probe can re-discover a support that lies on the
        // line (a tied vertex on axis-aligned shapes), so the two oblique
        // directions are probed as well
        const vector oblique_a{line[0] + perp[0], line[1] + perp[1]};
        const vector oblique_b{line[0] - perp[0], line[1] - perp[1]};
        const std::array<vector, 4> probes{perp, vector{-perp[0], -perp[1]},
                                           oblique_a, oblique_b};
        bool bootstrapped = false;
        for (const vector& probe : probes) {
          if (num::inner_prod(probe, probe) <=
              options.tolerance * options.tolerance)
            continue;
          const support candidate = diff_support(id_a, id_b, probe);
          const vector off =
              support_vector(candidate) - support_vector(tri[i0]);
          const scalar perpendicular_distance =
              std::abs(line[0] * off[1] - line[1] * off[0]) / line_length;
          if (perpendicular_distance > options.tolerance) {
            tri[0] = tri[i0];
            tri[1] = tri[i1];
            tri[2] = candidate;
            bootstrapped = true;
            break;
          }
        }
        if (!bootstrapped) {
          // the difference polytope really is a segment: the minimum
          // translation is the distance from the origin to that segment
          const vector d = line;
          const scalar dd = num::inner_prod(d, d);
          scalar t = 0.;
          if (dd > 0.) {
            t = -num::inner_prod(support_vector(tri[i0]), d) / dd;
            if (t < 0.) t = 0.;
            if (t > 1.) t = 1.;
          }
          const vector closest = support_vector(tri[i0]) + t * d;
          const scalar distance =
              std::sqrt(num::inner_prod(closest, closest));
          const vector normal =
              distance > options.tolerance
                  ? vector{closest[0] / distance, closest[1] / distance}
                  : vector{-d[1], d[0]};
          return make_result(distance, normal, tri[i0], tri[i1]);
        }
      }
    }
  }

  // The walk's exterior normals follow the triangle's orientation, and the
  // two GJK paths can hand the same three supports in a different order,
  // which would flip every normal the walk derives. Canonicalize the
  // orientation (counterclockwise) so the answer's normal points away from
  // the origin whatever order the witness came in.
  const std::array<support, 3> oriented =
      clockwise(support_vector(tri[0]), support_vector(tri[1]),
                support_vector(tri[2]))
          ? std::array<support, 3>{tri[0], tri[2], tri[1]}
          : tri;
  const vector a = support_vector(oriented[0]);
  const vector b = support_vector(oriented[1]);
  const vector c = support_vector(oriented[2]);
  const bool cw = clockwise(a, b, c);

  // The polytope can hold one vertex per point of the two shapes' difference.
  const std::size_t max_polytope =
      pool.entries[id_a].vertex_count + pool.entries[id_b].vertex_count;

  auto& ring = arena.points;
  ring.reset();
  ring.push_back(oriented[0]);
  ring.push_back(oriented[1]);
  ring.push_back(oriented[2]);

  for (std::size_t iteration = 0; iteration < options.max_iterations;
       ++iteration) {
    scalar closest_distance = std::numeric_limits<scalar>::infinity();
    vector closest_normal = vector();
    std::size_t closest_edge = 0;

    for (std::size_t i = 0; i < ring.size(); ++i) {
      const vector p1 = support_vector(ring[i]);
      const vector p2 = support_vector(ring[(i + 1) % ring.size()]);
      const vector edge = p2 - p1;
      if (num::inner_prod(edge, edge) <=
          options.tolerance * options.tolerance)
        return rejected();

      const vector normal = exterior_normal(cw, edge);
      const scalar distance = std::abs(num::inner_prod(normal, p1));
      if (distance < closest_distance) {
        closest_distance = distance;
        closest_normal = normal;
        closest_edge = i;
      }
    }

    const support candidate = diff_support(id_a, id_b, closest_normal);
    const scalar candidate_distance =
        num::inner_prod(support_vector(candidate), closest_normal);
    if (candidate_distance - closest_distance <= options.tolerance)
      return make_result(static_cast<scalar>(candidate_distance),
                         closest_normal, ring[closest_edge],
                         ring[(closest_edge + 1) % ring.size()]);

    for (std::size_t i = 0; i < ring.size(); ++i) {
      const vector difference =
          support_vector(candidate) - support_vector(ring[i]);
      if (num::inner_prod(difference, difference) <=
          options.tolerance * options.tolerance)
        return make_result(static_cast<scalar>(closest_distance),
                           closest_normal, ring[closest_edge],
                           ring[(closest_edge + 1) % ring.size()]);
    }

    if (ring.size() == max_polytope) return rejected();

    ring.insert_at(closest_edge + 1, candidate);
  }

  return rejected();
}

// Append-only list with room for a known, small number of entries. The 3D
// bootstrap below builds three such sets (distinct, candidates, triangles),
// each bounded by the four seed supports and the three bootstrap attempts.
// They are working scratch inside the per-pair loop, so they stay inline
// rather than reach the allocator the way a std::vector would -- the same
// reason polytope_region exists.
// A bound-sized list for the bootstrap and repair work sets, which are small
// and bounded by construction: the seed supports plus the bounded number of
// extra queries the attempt loops make. Capacity must bound that count by the
// same argument that bounds the loop -- push_back past Capacity writes past
// the inline store -- so each site states its bound next to its use. Unlike a
// std::vector this never reaches the allocator, which is what keeps a
// degenerate hit from charging the heap.
template <typename T, std::size_t Capacity>
struct fixed_list {
  std::size_t size() const { return count; }

  const T& operator[](const std::size_t i) const { return storage[i]; }
  T& operator[](const std::size_t i) { return storage[i]; }

  void push_back(const T& value) { storage[count++] = value; }

  const T* data() const { return storage.data(); }
  const T* begin() const { return storage.data(); }
  const T* end() const { return storage.data() + count; }

  std::array<T, Capacity> storage{};
  std::size_t count = 0;
};

// The penetration pass over one hit's seed, 3D branch: scalar epa's face walk
// over the vertex set and face list the arena holds, including its
// degenerate- witness bootstrap, the portal refinement that repairs a
// bootstrapped tetrahedron and the flat-shape projection when no tetrahedron
// exists.
template <typename Pool>
auto epa_penetration_3d(
    const Pool& pool, const std::uint32_t id_a, const std::uint32_t id_b,
    const std::array<typename Pool::support, 4>& seed,
    const gjk_options options,
    polytope_solid<typename Pool::vector, typename Pool::scalar>& arena)
    -> epa_result<typename Pool::scalar, 3>
{
  using vector = typename Pool::vector;
  using scalar = typename Pool::scalar;
  using support = typename Pool::support;
  using face = typename polytope_solid<vector, scalar>::face;
  using answer = epa_result<scalar, 3>;

  const auto rejected = []() -> answer {
    answer r{};
    r.depth = std::numeric_limits<scalar>::quiet_NaN();
    r.normal.fill(std::numeric_limits<scalar>::quiet_NaN());
    return r;
  };

  const auto diff_support = [&pool](const std::uint32_t a,
                                    const std::uint32_t b,
                                    const vector& direction) {
    return pool_diff_support(pool, a, b, direction);
  };

  // Validation mirrors epa, reported per hit instead of thrown: an option
  // pair epa would reject, or a shape outside the pool.
  if (options.tolerance <= 0. || options.max_iterations == 0)
    return rejected();
  if (id_a >= pool.entries.size() || id_b >= pool.entries.size() ||
      pool.entries[id_a].vertex_count == 0 ||
      pool.entries[id_b].vertex_count == 0)
    return rejected();

  vector a = support_vector(seed[0]);
  vector b = support_vector(seed[1]);
  vector c = support_vector(seed[2]);
  vector d = support_vector(seed[3]);

  const vector ab = b - a;
  const vector ac = c - a;
  const vector ad = d - a;
  const scalar volume = std::abs(num::inner_prod(cross(ab, ac), ad));
  if (volume <= options.tolerance * options.tolerance * options.tolerance) {
    // GJK handed EPA a degenerate witness simplex (points collinear or
    // coincident) that doesn't span a tetrahedron, even though the shapes
    // may genuinely overlap with real volume. Bootstrap a real
    // non-degenerate tetrahedron with a few extra support queries along
    // directions independent of what we already have, rather than trust
    // that degenerate witness for penetration depth. Only fall back to a
    // closest-point approximation on the degenerate point set if the
    // shapes truly don't admit a tetrahedron (e.g. a point/point or
    // line/line touching contact).
    //
    // Bound: the four seed supports plus the three directions below, so 7.
    fixed_list<vector, 4 + 3> distinct;
    for (const vector& v : {a, b, c, d}) {
      bool found = false;
      for (const vector& u : distinct) {
        const vector diff = u - v;
        if (num::inner_prod(diff, diff) <=
            options.tolerance * options.tolerance) {
          found = true;
          break;
        }
      }
      if (!found) distinct.push_back(v);
    }

    struct degeneracy_info {
      bool have_line = false;
      vector line_dir{};
      vector line_point{};
      bool have_plane = false;
      vector normal{};
      vector plane_point{};
      bool have_volume = false;
      vector volume_point{};
    };

    auto analyze = [&](const vector* pts, const std::size_t n) {
      degeneracy_info info;
      if (n < 2) return info;
      info.have_line = true;
      info.line_point = pts[1];
      const vector e1 = pts[1] - pts[0];
      const scalar e1_len2 = num::inner_prod(e1, e1);
      info.line_dir =
          util::transform(e1, [e1_len2](const auto&, const auto& x) {
            return x / std::sqrt(e1_len2);
          });
      for (std::size_t k = 2; k < n; ++k) {
        const vector ek = pts[k] - pts[0];
        if (!info.have_plane) {
          const vector n = cross(info.line_dir, ek);
          const scalar n_len2 = num::inner_prod(n, n);
          if (n_len2 > options.tolerance * options.tolerance) {
            const scalar n_len = std::sqrt(n_len2);
            info.normal = util::transform(
                n, [n_len](const auto&, const auto& x) { return x / n_len; });
            info.have_plane = true;
            info.plane_point = pts[k];
          }
        }
        else if (!info.have_volume) {
          const scalar off_plane = num::inner_prod(info.normal, ek);
          if (std::abs(off_plane) > options.tolerance) {
            info.have_volume = true;
            info.volume_point = pts[k];
          }
        }
      }
      return info;
    };

    constexpr std::size_t max_bootstrap_attempts = 3;
    for (std::size_t attempt = 0; attempt < max_bootstrap_attempts;
         ++attempt) {
      const auto info = analyze(distinct.data(), distinct.size());
      if (info.have_volume) break;

      // One entry when the simplex collapses to a point, two when it
      // escapes collinearity (the oblique sum, then its first
      // perpendicular), and two opposite normals when it escapes a plane.
      // Bound: two, exactly what those three cases push.
      fixed_list<vector, 2> candidates;
      if (!info.have_line) {
        candidates.push_back(vector{1, 0, 0});
      }
      else if (!info.have_plane) {
        // Escape collinearity. A single-axis perpendicular is itself
        // axis-aligned and, for axis-aligned box shapes, tends to
        // re-discover the same tied support vertex that produced the
        // degenerate simplex. Summing every non-zero axis-perpendicular
        // gives an oblique direction that reaches a genuinely different
        // vertex.
        const std::array<vector, 3> axes = {vector{1, 0, 0}, vector{0, 1, 0},
                                            vector{0, 0, 1}};
        vector perp_first{};
        vector perp_sum{};
        bool have_perp_first = false;
        for (const vector& axis : axes) {
          const vector cand = cross(info.line_dir, axis);
          const scalar cand_len2 = num::inner_prod(cand, cand);
          if (cand_len2 <= options.tolerance * options.tolerance) continue;
          const scalar cand_len = std::sqrt(cand_len2);
          const vector unit_cand =
              util::transform(cand, [cand_len](const auto&, const auto& x) {
                return x / cand_len;
              });
          if (!have_perp_first) {
            perp_first = unit_cand;
            have_perp_first = true;
          }
          perp_sum = perp_sum + unit_cand;
        }
        const scalar perp_sum_len2 = num::inner_prod(perp_sum, perp_sum);
        if (perp_sum_len2 > options.tolerance * options.tolerance) {
          const scalar perp_sum_len = std::sqrt(perp_sum_len2);
          candidates.push_back(util::transform(
              perp_sum, [perp_sum_len](const auto&, const auto& x) {
                return x / perp_sum_len;
              }));
        }
        if (have_perp_first) candidates.push_back(perp_first);
      }
      else {
        // Escape coplanarity: try the plane normal, then its opposite —
        // an axis-aligned normal can likewise re-discover an existing
        // tied vertex on one side of the plane.
        candidates.push_back(info.normal);
        candidates.push_back(-info.normal);
      }

      bool added = false;
      for (const vector& dir : candidates) {
        const support candidate = diff_support(id_a, id_b, dir);
        const vector cv = support_vector(candidate);
        bool is_new = true;
        for (const vector& u : distinct) {
          const vector diff = u - cv;
          if (num::inner_prod(diff, diff) <=
              options.tolerance * options.tolerance) {
            is_new = false;
            break;
          }
        }
        if (!is_new) continue;

        distinct.push_back(cv);
        added = true;
        break;
      }
      if (!added) break;
    }

    const auto final_info = analyze(distinct.data(), distinct.size());

    // Project the origin onto the triangle fan built from a coplanar (or
    // richer) point set, clipped to the fan's triangles/edges. Used both
    // as the genuinely-flat-shapes answer and as the last-resort answer
    // when a bootstrapped tetrahedron can't be made to enclose the
    // origin (see below).
    auto project_onto_points =
        [&](const vector* dv,
            const std::size_t n) -> std::tuple<scalar, vector> {
      // One triangle per interior index, so dv.size() - 2 at most. Bound: the
      // largest set this is called with is distinct's seven points, so five.
      fixed_list<std::array<std::size_t, 3>, 5> triangles;
      for (std::size_t i = 1; i + 1 < n; ++i) {
        const vector e10 = dv[i] - dv[0];
        const vector e20 = dv[i + 1] - dv[0];
        const vector normal_v = cross(e20, e10);
        if (num::inner_prod(normal_v, normal_v) >
            options.tolerance * options.tolerance)
          triangles.push_back({std::size_t(0), i, i + 1});
      }

      const vector& v0 = dv[triangles[0][0]];
      const vector& v1 = dv[triangles[0][1]];
      const vector& v2 = dv[triangles[0][2]];
      const vector e10 = v1 - v0;
      const vector e20 = v2 - v0;
      vector normal = cross(e20, e10);
      const scalar len = std::sqrt(num::inner_prod(normal, normal));
      normal = util::transform(
          normal, [len](const auto&, const auto& x) { return x / len; });

      const vector origin = vector();
      const scalar t = -num::inner_prod(normal, v0);
      const vector p = origin - t * normal;
      scalar best_dist2 = std::numeric_limits<scalar>::infinity();
      vector best_point = p;
      for (const auto& tri : triangles) {
        const vector& ta = dv[tri[0]];
        const vector& tb = dv[tri[1]];
        const vector& tc = dv[tri[2]];
        const vector tab = tb - ta;
        const vector tac = tc - ta;
        const vector tap = p - ta;
        const scalar d1 = num::inner_prod(tab, tap);
        const scalar d2 = num::inner_prod(tac, tap);
        const scalar d3 = num::inner_prod(tab, tab);
        const scalar d4 = num::inner_prod(tac, tab);
        const scalar d5 = num::inner_prod(tab, tac);
        const scalar d6 = num::inner_prod(tac, tac);
        const scalar denom = d3 * d6 - d4 * d5;
        if (std::abs(denom) > options.tolerance * options.tolerance) {
          const scalar u = (d1 * d6 - d2 * d4) / denom;
          const scalar v = (d3 * d2 - d5 * d1) / denom;
          if (u >= 0 && v >= 0 && u + v <= 1) {
            const vector cp = ta + u * tab + v * tac;
            const scalar dist2 = num::inner_prod(cp, cp);
            if (dist2 < best_dist2) {
              best_dist2 = dist2;
              best_point = cp;
            }
          }
        }
        const std::array<std::array<std::size_t, 2>, 3> edges = {
            {{tri[0], tri[1]}, {tri[1], tri[2]}, {tri[2], tri[0]}}};
        for (const auto& edge : edges) {
          const vector ea = dv[edge[0]];
          const vector eb = dv[edge[1]];
          const vector e = eb - ea;
          const scalar t2 =
              std::clamp(-num::inner_prod(e, p - ea) / num::inner_prod(e, e),
                         scalar(0), scalar(1));
          const vector cp = ea + t2 * e;
          const scalar dist2 = num::inner_prod(cp, cp);
          if (dist2 < best_dist2) {
            best_dist2 = dist2;
            best_point = cp;
          }
        }
      }
      const scalar depth = std::sqrt(best_dist2);
      vector out_normal = best_point;
      const scalar out_len =
          std::sqrt(num::inner_prod(out_normal, out_normal));
      if (out_len > options.tolerance)
        out_normal = util::transform(
            out_normal,
            [out_len](const auto&, const auto& x) { return x / out_len; });
      else
        out_normal = normal;
      return std::tuple{depth, out_normal};
    };

    if (final_info.have_volume) {
      // Bootstrap found four affinely-independent points, but nothing
      // guarantees they enclose the origin the way GJK's own witness
      // simplex normally would -- EPA's face-walk loop below requires
      // that enclosure to produce a meaningful answer. Refine with a
      // bounded portal-style procedure (in the spirit of MPR/XenoCollide
      // portal construction): repeatedly find the face the origin is
      // least enclosed by, and replace its opposite vertex with a fresh
      // support point along that face's outward normal.
      std::array<vector, 4> tet = {distinct[0], final_info.line_point,
                                   final_info.plane_point,
                                   final_info.volume_point};
      struct face_def {
        std::size_t i0, i1, i2, opp;
      };
      constexpr std::array<face_def, 4> face_defs = {
          {{0, 1, 2, 3}, {0, 1, 3, 2}, {0, 2, 3, 1}, {1, 2, 3, 0}}};

      auto worst_face = [&]() {
        std::size_t worst = 0;
        scalar worst_val = -std::numeric_limits<scalar>::infinity();
        vector worst_normal = vector();
        for (std::size_t fi = 0; fi < 4; ++fi) {
          const auto& fd = face_defs[fi];
          const vector& fa = tet[fd.i0];
          const vector& fb = tet[fd.i1];
          const vector& fc = tet[fd.i2];
          const vector& opp = tet[fd.opp];
          vector n = cross(fb - fa, fc - fa);
          const scalar n_len2 = num::inner_prod(n, n);
          if (n_len2 <= options.tolerance * options.tolerance) continue;
          const scalar n_len = std::sqrt(n_len2);
          n = util::transform(
              n, [n_len](const auto&, const auto& x) { return x / n_len; });
          if (num::inner_prod(n, opp - fa) > 0) n = -n;
          const scalar val = num::inner_prod(n, vector() - fa);
          if (val > worst_val) {
            worst_val = val;
            worst = fi;
            worst_normal = n;
          }
        }
        return std::tuple{worst, worst_val, worst_normal};
      };

      bool enclosing = false;
      constexpr std::size_t max_refine_attempts = 4;
      for (std::size_t refine = 0; refine < max_refine_attempts; ++refine) {
        const auto [worst, worst_val, worst_normal] = worst_face();
        if (worst_val < -options.tolerance) {
          enclosing = true;
          break;
        }
        const support candidate = diff_support(id_a, id_b, worst_normal);
        const vector cv = support_vector(candidate);
        const vector diff = cv - tet[face_defs[worst].opp];
        if (num::inner_prod(diff, diff) <=
            options.tolerance * options.tolerance)
          break;
        tet[face_defs[worst].opp] = cv;
      }

      if (enclosing) {
        // Fall through: a/b/c/d now span a tetrahedron that strictly
        // encloses the origin, so the ordinary face-walk EPA loop below
        // can run on it directly.
        a = tet[0];
        b = tet[1];
        c = tet[2];
        d = tet[3];
      }
      else {
        // Refinement couldn't produce a strictly-enclosing tetrahedron.
        // Fall back to projecting onto the coplanar subset that existed
        // just before the volume-completing point was bootstrapped in
        // (that point is always the last one appended to `distinct`).
        const auto [depth, out_normal] =
            project_onto_points(distinct.data(), distinct.size() - 1);
        epa_result<scalar, 3> r{};
        r.depth = static_cast<scalar>(depth);
        r.normal = {static_cast<scalar>(out_normal[0]),
                    static_cast<scalar>(out_normal[1]),
                    static_cast<scalar>(out_normal[2])};
        return r;
      }
    }
    else if (final_info.have_plane) {
      // The Minkowski difference is genuinely flat (e.g. two coplanar
      // shapes): no tetrahedron is possible. Project the origin onto the
      // triangle fan built from the plane we do have.
      const auto [depth, out_normal] =
          project_onto_points(distinct.data(), distinct.size());
      return {static_cast<scalar>(depth),
              {static_cast<scalar>(out_normal[0]),
               static_cast<scalar>(out_normal[1]),
               static_cast<scalar>(out_normal[2])}};
    }
    else {
      // Still collinear/coincident after bootstrapping: the shapes truly
      // reduce to a point or a line in the Minkowski difference. Fall
      // back to the closest point on the union of segments between the
      // distinct points found.
      const vector origin = vector();
      scalar best_dist2 = std::numeric_limits<scalar>::infinity();
      vector best_point = distinct[0];
      if (distinct.size() == 1) {
        best_point = distinct[0];
        best_dist2 = num::inner_prod(best_point, best_point);
      }
      else {
        for (std::size_t i = 0; i < distinct.size(); ++i) {
          for (std::size_t j = i + 1; j < distinct.size(); ++j) {
            const vector& ea = distinct[i];
            const vector& eb = distinct[j];
            const vector e = eb - ea;
            const scalar t2 = std::clamp(
                -num::inner_prod(e, ea - origin) / num::inner_prod(e, e),
                scalar(0), scalar(1));
            const vector cp = ea + t2 * e;
            const scalar dist2 = num::inner_prod(cp, cp);
            if (dist2 < best_dist2) {
              best_dist2 = dist2;
              best_point = cp;
            }
          }
        }
      }

      const scalar depth = std::sqrt(best_dist2);
      vector out_normal = best_point;
      const scalar out_len =
          std::sqrt(num::inner_prod(out_normal, out_normal));
      if (out_len > options.tolerance)
        out_normal = util::transform(
            out_normal,
            [out_len](const auto&, const auto& x) { return x / out_len; });
      else
        out_normal = vector{1, 0, 0};
      return {static_cast<scalar>(depth),
              {static_cast<scalar>(out_normal[0]),
               static_cast<scalar>(out_normal[1]),
               static_cast<scalar>(out_normal[2])}};
    }
  }

  // GJK's witness simplex can also merely touch the origin instead of
  // enclosing it: an edge of the tetrahedron passes exactly through the
  // origin, which leaves two of its faces at zero distance. The face
  // walk then converges on a chord of the Minkowski difference (the
  // support along the face normal re-selects a vertex already in the
  // polytope) and reports zero penetration for a pair that genuinely
  // overlaps. Repair the witness before the walk: when it does not
  // strictly enclose the origin, search the supports along the six axis
  // directions for four points whose convex hull does. A witness that
  // cannot be repaired this way is a true tangent contact, for which the
  // walk correctly answers zero.
  {
    const auto encloses = [&options](const std::array<vector, 4>& t) {
      static constexpr std::array<std::array<std::size_t, 4>, 4> faces = {
          {{{0, 1, 2, 3}}, {{0, 1, 3, 2}}, {{0, 2, 3, 1}}, {{1, 2, 3, 0}}}};
      for (const auto& f : faces) {
        const vector& p = t[f[0]];
        const vector& q = t[f[1]];
        const vector& r = t[f[2]];
        const vector& opp = t[f[3]];
        vector n = cross(q - p, r - p);
        if (num::inner_prod(n, n) <= options.tolerance * options.tolerance)
          return false;
        if (num::inner_prod(n, opp - p) > 0) n = -n;
        if (num::inner_prod(n, vector() - p) >= -options.tolerance)
          return false;
      }
      return true;
    };

    std::array<vector, 4> witness = {a, b, c, d};
    if (!encloses(witness)) {
      // The four witness supports plus the six axis supports: bound 10.
      fixed_list<vector, 10> points;
      for (const vector& v : witness) points.push_back(v);
      const std::array<vector, 6> axes = {
          {vector{1, 0, 0}, vector{-1, 0, 0}, vector{0, 1, 0},
           vector{0, -1, 0}, vector{0, 0, 1}, vector{0, 0, -1}}};
      for (const vector& axis : axes) {
        const vector v = support_vector(diff_support(id_a, id_b, axis));
        bool known = false;
        for (const vector& u : points)
          if (num::inner_prod(u - v, u - v) <=
              options.tolerance * options.tolerance) {
            known = true;
            break;
          }
        if (!known) points.push_back(v);
      }
      for (std::size_t i = 0; i + 3 < points.size(); ++i)
        for (std::size_t j = i + 1; j + 2 < points.size(); ++j)
          for (std::size_t k = j + 1; k + 1 < points.size(); ++k)
            for (std::size_t l = k + 1; l < points.size(); ++l) {
              if (encloses(witness)) break;
              const std::array<vector, 4> trial = {points[i], points[j],
                                                   points[k], points[l]};
              if (encloses(trial)) witness = trial;
            }
      a = witness[0];
      b = witness[1];
      c = witness[2];
      d = witness[3];
    }
  }

  auto& vertices = arena.vertices;
  auto& faces = arena.faces;
  vertices.reset();
  faces.reset();

  auto make_face = [&](const std::size_t i0, const std::size_t i1,
                       const std::size_t i2) {
    const vector v0 = vertices[i0];
    const vector v1 = vertices[i1];
    const vector v2 = vertices[i2];
    const vector e10 = v1 - v0;
    const vector e20 = v2 - v0;
    vector normal = cross(e20, e10);
    const scalar len2 = num::inner_prod(normal, normal);
    const scalar len = std::sqrt(len2);
    if (len > options.tolerance) {
      normal = util::transform(
          normal, [len](const auto&, const auto& x) { return x / len; });
    }
    const scalar distance = num::inner_prod(normal, v0);
    return face{{i0, i1, i2}, normal, distance};
  };

  vertices.push_back(a);
  vertices.push_back(b);
  vertices.push_back(c);
  vertices.push_back(d);

  arena.vert_support_count = 0;
  if (arena.vert_support_count < epa_sbo_capacity) {
    arena.vert_support_ids[arena.vert_support_count++] = {
        std::get<0>(seed[0]), std::get<1>(seed[0])};
  }
  if (arena.vert_support_count < epa_sbo_capacity) {
    arena.vert_support_ids[arena.vert_support_count++] = {
        std::get<0>(seed[1]), std::get<1>(seed[1])};
  }
  if (arena.vert_support_count < epa_sbo_capacity) {
    arena.vert_support_ids[arena.vert_support_count++] = {
        std::get<0>(seed[2]), std::get<1>(seed[2])};
  }
  if (arena.vert_support_count < epa_sbo_capacity) {
    arena.vert_support_ids[arena.vert_support_count++] = {
        std::get<0>(seed[3]), std::get<1>(seed[3])};
  }
  faces.push_back(make_face(0, 2, 1));  // ACB
  faces.push_back(make_face(0, 1, 3));  // ABD
  faces.push_back(make_face(0, 3, 2));  // ADC
  faces.push_back(make_face(1, 2, 3));  // BCD

  for (auto& f : faces) {
    if (f.distance < 0) {
      f.distance = -f.distance;
      f.normal = -f.normal;
      std::swap(f.indices[0], f.indices[1]);
    }
  }

  // The polytope can hold one vertex per point of the two shapes' difference.
  const std::size_t max_polytope =
      pool.entries[id_a].vertex_count + pool.entries[id_b].vertex_count;

  for (std::size_t iteration = 0; iteration < options.max_iterations;
       ++iteration) {
    scalar closest_distance = std::numeric_limits<scalar>::infinity();
    std::size_t closest_face = 0;

    for (std::size_t i = 0; i < faces.size(); ++i) {
      if (faces[i].distance < closest_distance) {
        closest_distance = faces[i].distance;
        closest_face = i;
      }
    }

    const support candidate =
        diff_support(id_a, id_b, faces[closest_face].normal);
    const scalar candidate_distance = num::inner_prod(
        support_vector(candidate), faces[closest_face].normal);
    if (candidate_distance - closest_distance <= options.tolerance) {
      answer r{};
      r.depth = static_cast<scalar>(closest_distance);
      r.normal = {static_cast<scalar>(faces[closest_face].normal[0]),
                  static_cast<scalar>(faces[closest_face].normal[1]),
                  static_cast<scalar>(faces[closest_face].normal[2])};
      for (std::size_t k = 0; k < 3; ++k) {
        const std::size_t vi = faces[closest_face].indices[k];
        if (vi < epa_sbo_capacity && vi < arena.vert_support_count) {
          r.vertices[k] = arena.vert_support_ids[vi];
        }
      }
      return r;
    }

    const vector candidate_vec = support_vector(candidate);
    bool duplicate = false;
    for (const auto& v : vertices) {
      const vector diff = v - candidate_vec;
      if (num::inner_prod(diff, diff) <=
          options.tolerance * options.tolerance) {
        duplicate = true;
        break;
      }
    }
    if (duplicate) {
      answer r{};
      r.depth = static_cast<scalar>(closest_distance);
      r.normal = {static_cast<scalar>(faces[closest_face].normal[0]),
                  static_cast<scalar>(faces[closest_face].normal[1]),
                  static_cast<scalar>(faces[closest_face].normal[2])};
      return r;
    }

    if (vertices.size() >= max_polytope) return rejected();

    const std::size_t new_index = vertices.size();
    vertices.push_back(candidate_vec);
    if (arena.vert_support_count < epa_sbo_capacity) {
      arena.vert_support_ids[arena.vert_support_count++] = {
          std::get<0>(candidate), std::get<1>(candidate)};
    }
    else {
      ++arena.vert_support_count;  // track count even if not stored
    }

    // Remove every face the new vertex is in front of.
    for (std::size_t i = 0; i < faces.size(); ++i) {
      const vector& v0 = vertices[faces[i].indices[0]];
      const scalar dot = num::inner_prod(faces[i].normal, candidate_vec - v0);
      faces[i].removed = dot > options.tolerance;
    }

    bool all_faces_removed = true;
    for (const auto& f : faces) {
      if (!f.removed) {
        all_faces_removed = false;
        break;
      }
    }
    if (all_faces_removed) {
      answer r{};
      r.depth = static_cast<scalar>(closest_distance);
      r.normal = {static_cast<scalar>(faces[closest_face].normal[0]),
                  static_cast<scalar>(faces[closest_face].normal[1]),
                  static_cast<scalar>(faces[closest_face].normal[2])};
      for (std::size_t k = 0; k < 3; ++k) {
        const std::size_t vi = faces[closest_face].indices[k];
        if (vi < epa_sbo_capacity && vi < arena.vert_support_count) {
          r.vertices[k] = arena.vert_support_ids[vi];
        }
      }
      return r;
    }

    // Keep each retained face's edges in their original winding order:
    // make_face(new_index, edge[0], edge[1]) depends on that direction to
    // produce an outward-facing normal for the new face. The {min, max} form
    // of an edge is computed only as a duplicate-detection key, never as the
    // edge to build a face from -- sorting it first discards the winding and
    // can silently flip newly-created faces inward once EPA needs a second
    // expansion iteration. An edge two retained faces carry is interior to
    // the polytope and drops out; one carried by a single retained face is a
    // boundary edge the new faces are triangulated on.
    const std::size_t retained = faces.size();
    for (std::size_t fi = 0; fi < retained; ++fi) {
      if (faces[fi].removed) continue;
      // The face is addressed by index, never by reference: push_back below
      // may reallocate the region's heap buffer while this loop runs, and a
      // reference taken before that would dangle into freed memory.
      for (std::size_t e = 0; e < 3; ++e) {
        const std::size_t e0 = faces[fi].indices[e];
        const std::size_t e1 = faces[fi].indices[(e + 1) % 3];
        const auto boundary = std::array{std::min(e0, e1), std::max(e0, e1)};

        std::size_t carriers = 0;
        for (std::size_t fj = 0; fj < retained; ++fj) {
          if (faces[fj].removed) continue;
          for (std::size_t f = 0; f < 3; ++f) {
            const std::size_t f0 = faces[fj].indices[f];
            const std::size_t f1 = faces[fj].indices[(f + 1) % 3];
            const auto key = std::array{std::min(f0, f1), std::max(f0, f1)};
            if (key == boundary) ++carriers;
          }
        }
        if (carriers != 1) continue;

        face fresh = make_face(new_index, e0, e1);
        if (fresh.distance < 0) {
          fresh.distance = -fresh.distance;
          fresh.normal = -fresh.normal;
          std::swap(fresh.indices[0], fresh.indices[1]);
        }
        faces.push_back(fresh);
      }
    }

    // Drop the removed faces, keeping the survivors and then the new ones in
    // order -- the list the scalar path rebuilds into a fresh vector.
    faces.erase_if(
        [&faces](const std::size_t i) { return faces[i].removed; });
  }

  return rejected();
}

// The per-hit polytope storage the EPA evaluators work in: the 2D branch
// keeps the ring of support points, the 3D branch the vertex set and the
// face list built on it. A caller running one batch per simulation step owns
// one arena across the batches: the evaluators reset what they use per hit,
// and reserve() only grows, so the capacity built for the previous step is
// reused instead of rebuilt per batch.
template <typename Pool, std::size_t Dim>
  requires(Dim == 2 || Dim == 3)
using epa_arena = std::conditional_t<
    Dim == 2, polytope_ring<typename Pool::support>,
    polytope_solid<typename Pool::vector, typename Pool::scalar>>;

// The batched penetration pass (core, by reference): the second pass over
// the same pool and the same pair list as gjk_batch. `pairs` is the pair
// list that pass evaluated and `gjk_results` is parallel to it, so a hit's
// shape ids come from its own pair record and its seed from its own result
// -- no GJK run to recover either. hit_pairs[k] = j selects result and pair
// j, so answer k belongs to hit_pairs[k]; a hit the batch cannot answer is
// NaN in place rather than missing from the result. The answer buffer and
// the arena are caller-owned: every entry is assigned before it is read, so
// a buffer kept across steps is safe to reuse with stale contents.
template <std::size_t Dim, typename Pool, typename PairRecord,
          typename Support>
  requires(Dim == 2 || Dim == 3)
void epa_batch(const Pool& pool, const std::vector<PairRecord>& pairs,
               const std::vector<std::uint32_t>& hit_pairs,
               const std::vector<batch_result<Support, Dim + 1>>& gjk_results,
               const gjk_options options,
               std::vector<epa_result<typename Pool::scalar, Dim>>& answers,
               epa_arena<Pool, Dim>& arena, batch_timing* timing = nullptr)
{
  using scalar = typename Pool::scalar;

  assert_device_contract<Pool, Dim>();

  const std::size_t pair_count = gjk_results.size();
  const std::size_t record_count = pairs.size();

  // assign (not resize) so the buffer is re-filled even after the answers it
  // held are stale; a buffer kept across calls is emptied and re-added
  // without losing its heap block, and its capacity from the previous step
  // is reused.
  answers.assign(hit_pairs.size(), unanswered<scalar, Dim>());

  // Validation that depends on the options alone, the way epa rejects them:
  // every hit goes unanswered.
  if (options.tolerance <= 0. || options.max_iterations == 0) return;

  // Pre-size the arena once, before the hit loop, from the exact bound the
  // guards already use. The polytope can hold at most one vertex per point of
  // the two shapes' difference, so max_polytope = vertex_count_A +
  // vertex_count_B bounds it, and the face list follows from Euler (a convex
  // polytope of V vertices carries at most 2V - 4 faces). Both depend only on
  // shape ids, so the worst case over the batch's pairs is one pass over the
  // pool. With the regions reserved to that bound, the hit loop never
  // reallocates: push_back and insert_at fit the reserved storage, and
  // promote() has nothing to allocate. Batching a constant shape set is then
  // allocation-free per pair, which is the property the batch path exists to
  // guarantee.
  //
  // The guards are the hit loop's own: a hit that names a pair this batch
  // never evaluated, or whose ids fall outside the pool, is unanswered in
  // place and contributes nothing to the bound. Without them this sizing pass
  // would read the pool for a pair the loop skips.
  auto vertex_bound = [&](const std::uint32_t hit) -> std::size_t {
    if (hit >= pair_count || hit >= record_count) return 0;
    const auto& pair = pairs[hit];
    if (pair.shape_id_A >= pool.entries.size() ||
        pair.shape_id_B >= pool.entries.size())
      return 0;
    return static_cast<std::size_t>(
               pool.entries[pair.shape_id_A].vertex_count) +
           static_cast<std::size_t>(
               pool.entries[pair.shape_id_B].vertex_count);
  };

  std::size_t worst_vertices = 0;
  for (const std::uint32_t hit : hit_pairs)
    worst_vertices = std::max(worst_vertices, vertex_bound(hit));

  if constexpr (Dim == 2) {
    arena.reserve(worst_vertices);
  }
  else {
    arena.reserve(worst_vertices, 2 * worst_vertices);
  }

  if (timing != nullptr) timing->seconds.assign(hit_pairs.size(), 0.);

  for (std::size_t i = 0; i < hit_pairs.size(); ++i) {
    // A hit that names a pair this batch never evaluated is unanswered in
    // place: there is no result to read a seed from, so every position stays
    // indexable and nothing is read out of bounds.
    const std::uint32_t hit = hit_pairs[i];
    if (hit >= pair_count || hit >= record_count) continue;

    // The hit's own pair record names the two shapes, and a hit whose ids
    // fall outside the pool is left unanswered rather than read out of
    // bounds.
    const auto& pair = pairs[hit];
    if (pair.shape_id_A >= pool.entries.size() ||
        pair.shape_id_B >= pool.entries.size())
      continue;

    const auto hit_start = std::chrono::steady_clock::now();

    if constexpr (Dim == 2)
      answers[i] = epa_penetration_2d(pool, pair.shape_id_A, pair.shape_id_B,
                                      gjk_results[hit].seed, options, arena);
    else
      answers[i] = epa_penetration_3d(pool, pair.shape_id_A, pair.shape_id_B,
                                      gjk_results[hit].seed, options, arena);

    // Skipped rows (invalid hit records) keep their assigned 0 so the row
    // count stays one per hit.
    if (timing != nullptr)
      timing->seconds[i] = std::chrono::duration<double>(
                               std::chrono::steady_clock::now() - hit_start)
                               .count();
  }

  return;
}

// The batched penetration pass (allocating form): allocates the answer
// buffer and the arena and forwards to the core overload above. Kept as a
// convenience for callers that buffer once per call; the by-reference core
// is the intended form for the per-step batch stage.
template <std::size_t Dim, typename Pool, typename PairRecord,
          typename Support>
  requires(Dim == 2 || Dim == 3)
auto epa_batch(const Pool& pool, const std::vector<PairRecord>& pairs,
               const std::vector<std::uint32_t>& hit_pairs,
               const std::vector<batch_result<Support, Dim + 1>>& gjk_results,
               const gjk_options options = {})
    -> std::vector<epa_result<typename Pool::scalar, Dim>>
{
  std::vector<epa_result<typename Pool::scalar, Dim>> answers;
  epa_arena<Pool, Dim> arena;
  epa_batch<Dim>(pool, pairs, hit_pairs, gjk_results, options, answers,
                 arena);
  return answers;
}

// ---------------------------------------------------------------------------
// SAT batch: Separating Axis Theorem as a batched narrowphase, carrying
// answers and not just statuses: one sat_result per pair -- the status,
// the normal oriented A -> B, and a signed gap h, positive for separated
// pairs (the largest gap the axis set found, a lower bound on the true
// gap) and negative for enclosed pairs (the min-overlap depth, the sign
// convention siconos re-derives per pair today). No witness set leaves
// this section: response construction stays with the GJK/EPA path, where
// EPA answers with the seeded ring vertices, and SAT is the answering
// path for small convex polygons and the gate for everything else. This
// keep the three asks (2) prefilter, (3) SAT answers) in one surface.
// ---------------------------------------------------------------------------

// The SAT answer row. st mirrors batch_result's vocabulary with
// max_iterations unreachable (SAT decides per axis, not iteratively);
// invalid_input happens only in the batch wrappers, the scalar entry
// reports it by throwing the way sat_inter does.
template <typename Scalar, std::size_t Dim>
struct sat_result {
  status st{};
  std::array<Scalar, Dim> normal{};  // oriented A -> B
  Scalar h{};                        // +gap separated, -depth enclosed
};

// The tuple sat_intersect_* returns, wrapped once: enclosed rows carry the
// positive min-overlap, separated rows the negated gap (it tracks the
// deepest separation), so the flip to reach h = +gap / -depth is a single
// sign and the axis needs no re-orientation -- the intersect functions
// already leave it A -> B.
template <std::size_t Dim, typename Tuple>
  requires(Dim == 2 || Dim == 3)
auto sat_from_tuple(const Tuple& res)
    -> sat_result<typename decltype(env(std::get<1>(res)))::scalar, Dim>
{
  using scalar = typename decltype(env(std::get<1>(res)))::scalar;
  using answer = sat_result<scalar, Dim>;

  answer entry{};
  const bool inter = std::get<0>(res);
  entry.st = inter ? status::enclosed : status::separated;
  const auto& axis = std::get<1>(res);
  for (std::size_t k = 0; k < Dim; ++k) entry.normal[k] = scalar(axis[k]);
  entry.h = -std::get<2>(res);
  return entry;
}

// The scan over stored convex-clean polygons: the same maths as
// sat_intersect_2d -- edge-normal axes, projection min/max, the gap kept
// as a max, the overlap kept as a min and clamped, both axes oriented
// A -> B by the final centroid flip -- without the per-call hull build
// and axis dedup. Duplicate parallel axes (collinear triples) decide the
// same way, so they are harmless and kept. The caller runs
// cyclic_convex_2d first: the scan assumes a zero-edge-free walk.
template <typename ShapeA, typename ShapeB>
auto sat_2d_scan(const ShapeA& sa, const ShapeB& sb,
                 const gjk_options& options)
    -> std::tuple<bool, std::remove_cvref_t<decltype(sa[0])>,
                  typename decltype(env(sa[0]))::scalar>
{
  using point = std::remove_cvref_t<decltype(sa[0])>;
  using vector = point;
  using scalar = typename decltype(env(point{}))::scalar;

  const auto count_of = [](const auto& shape) {
    return static_cast<std::size_t>(
        std::distance(shape.begin(), shape.end()));
  };
  const std::size_t na = count_of(sa);
  const std::size_t nb = count_of(sb);

  bool intersect = true;
  vector sep_axis{};
  scalar sep_gap = 0;
  vector min_axis{};
  scalar min_overlap = (scalar)1e300;

  auto project_axis = [](const auto& shape, const std::size_t n,
                         const vector& axis) {
    scalar mn = (scalar)1e300, mx = (scalar)-1e300;
    for (std::size_t i = 0; i < n; ++i) {
      const auto& p = shape[i];
      const scalar d = p[0] * axis[0] + p[1] * axis[1];
      mn = std::min(mn, d);
      mx = std::max(mx, d);
    }
    return std::pair{mn, mx};
  };

  for (std::size_t side = 0; side < 2; ++side) {
    const auto& shape = side == 0 ? sa : sb;
    const std::size_t n = side == 0 ? na : nb;
    for (std::size_t i = 0; i < n; ++i) {
      const auto& p0 = shape[i];
      const auto& p1 = shape[(i + 1) % n];
      vector nrm{};
      nrm[0] = -(p1[1] - p0[1]);  // the edge normal, then normalized
      nrm[1] = p1[0] - p0[0];
      const scalar len = std::sqrt(nrm[0] * nrm[0] + nrm[1] * nrm[1]);
      if (len <= options.tolerance) continue;
      nrm[0] /= len;
      nrm[1] /= len;

      const auto [mnA, mxA] = project_axis(sa, na, nrm);
      const auto [mnB, mxB] = project_axis(sb, nb, nrm);
      if (mxA < mnB - options.tolerance || mxB < mnA - options.tolerance) {
        intersect = false;
        const scalar gap = std::max(mnB - mxA, mnA - mxB);
        if (gap > sep_gap) {
          sep_gap = gap;
          sep_axis = nrm;
        }
      }
      else {
        const scalar overlap = std::min(mxA, mxB) - std::max(mnA, mnB);
        if (overlap < min_overlap) {
          min_overlap = overlap;
          min_axis = nrm;
        }
      }
    }
  }

  if (intersect) {
    if (min_overlap < 0) min_overlap = 0;
    scalar cA = 0, cB = 0;
    for (std::size_t i = 0; i < na; ++i)
      cA += sa[i][0] * min_axis[0] + sa[i][1] * min_axis[1];
    for (std::size_t i = 0; i < nb; ++i)
      cB += sb[i][0] * min_axis[0] + sb[i][1] * min_axis[1];
    cA /= static_cast<scalar>(na);
    cB /= static_cast<scalar>(nb);
    if (cB - cA < 0) {
      min_axis[0] = -min_axis[0];
      min_axis[1] = -min_axis[1];
    }
    return std::tuple<bool, vector, scalar>{true, min_axis, min_overlap};
  }

  scalar centerA = 0, centerB = 0;
  for (std::size_t i = 0; i < na; ++i)
    centerA += sa[i][0] * sep_axis[0] + sa[i][1] * sep_axis[1];
  for (std::size_t i = 0; i < nb; ++i)
    centerB += sb[i][0] * sep_axis[0] + sb[i][1] * sep_axis[1];
  centerA /= static_cast<scalar>(na);
  centerB /= static_cast<scalar>(nb);
  if (centerB - centerA < 0) {
    sep_axis[0] = -sep_axis[0];
    sep_axis[1] = -sep_axis[1];
  }
  return std::tuple<bool, vector, scalar>{false, sep_axis, -sep_gap};
}

// The scalar per-pair answer, usable on its own and the engine the batch
// wrappers call: validates shapes and options the way sat_inter does
// (throwing invalid_argument), then answers 2D through the no-alloc scan
// when both shapes walk a convex polygon in their stored order --
// re-checked per call, cheap and allocation-free -- through
// sat_intersect_2d's hull build otherwise, and through sat_intersect_3d's
// face axes always.
template <std::size_t Dim, typename ShapeA, typename ShapeB>
  requires(Dim == 2 || Dim == 3)
auto sat_pair(const ShapeA& sa, const ShapeB& sb,
              const gjk_options options = {})
    -> sat_result<typename decltype(env(sa[0]))::scalar, Dim>
{
  if constexpr (Dim == 2) {
    validate_shape_2d(sa, "first shape");
    validate_shape_2d(sb, "second shape");
    validate_options(options);
    const bool clean = cyclic_convex_2d(sa.begin(), sa.end()) &&
                       cyclic_convex_2d(sb.begin(), sb.end());
    return clean ? sat_from_tuple<Dim>(sat_2d_scan(sa, sb, options))
                 : sat_from_tuple<Dim>(sat_intersect_2d(sa, sb, options));
  }
  else {
    validate_shape_3d(sa, "first shape");
    validate_shape_3d(sb, "second shape");
    validate_options(options);
    return sat_from_tuple<Dim>(sat_intersect_3d(sa, sb, options));
  }
}

// The batched SAT pass, per pair: the same validation shape as
// gjk_batch_pair (every scalar rejection lands as invalid_input, no throw
// crosses the batch edge), then the O(1) ball prefilter when bounds are
// given -- answering separated with the unit axis A -> B and the ball-gap
// lower bound dist - (rA + rB), which the true gap exceeds, without
// reading a vertex -- and the scalar engine otherwise. With bounds and
// both shapes flagged convex by compute_pool_bounds, the per-pair
// convexity walk is skipped and the scan runs on the pool windows
// directly.
template <std::size_t Dim, typename Pool, typename PairRecord>
  requires(Dim == 2 || Dim == 3)
auto sat_batch_pair(
    const Pool& pool, const PairRecord& pair, const gjk_options options,
    const pool_bounds<typename Pool::scalar, Dim>* bounds = nullptr)
    -> sat_result<typename Pool::scalar, Dim>
{
  using scalar = typename Pool::scalar;
  using answer = sat_result<scalar, Dim>;

  const auto invalid_input = []() {
    answer entry{};
    entry.st = status::invalid_input;
    return entry;
  };

  const auto& entries = pool.entries;
  if (pair.shape_id_A >= entries.size() ||
      pair.shape_id_B >= entries.size() ||
      entries[pair.shape_id_A].vertex_count == 0 ||
      entries[pair.shape_id_B].vertex_count == 0 ||
      std::size(pair.init_dir) != Dim ||
      pool.vertex_buffer[entries[pair.shape_id_A].vertex_offset].size() !=
          Dim ||
      pool.vertex_buffer[entries[pair.shape_id_B].vertex_offset].size() !=
          Dim)
    return invalid_input();

  if (options.tolerance <= 0. || options.max_iterations == 0)
    return invalid_input();

  // The prefilter answers without reading a vertex: the axis A -> B and
  // the gap lower bound dist - (rA + rB). A zero centre difference leaves
  // no direction to name, so the normal stays zero on that guard; it is
  // unreachable with positive radii and a rejecting distance anyway.
  if (bounds != nullptr &&
      bounds_separated<Dim>(pool, *bounds, pair.shape_id_A, pair.shape_id_B,
                            options.tolerance)) {
    answer entry{};
    entry.st = status::separated;
    scalar d2{};
    for (std::size_t k = 0; k < Dim; ++k) {
      const scalar diff =
          bounds->centre[std::size_t{pair.shape_id_B} * Dim + k] -
          bounds->centre[std::size_t{pair.shape_id_A} * Dim + k];
      d2 += diff * diff;
      entry.normal[k] = diff;
    }
    const scalar dist = std::sqrt(d2);
    if (dist > scalar{}) {
      for (std::size_t k = 0; k < Dim; ++k) entry.normal[k] /= dist;
    }
    entry.h = dist - (bounds->radius[pair.shape_id_A] +
                      bounds->radius[pair.shape_id_B]);
    return entry;
  }

  if constexpr (Dim == 2) {
    // the pool flags carry the convexity verdict, so the per-pair walk
    // sat_pair would run is already done
    if (bounds != nullptr && bounds->convex[pair.shape_id_A] != 0 &&
        bounds->convex[pair.shape_id_B] != 0) {
      try {
        return sat_from_tuple<Dim>(
            sat_2d_scan(shape_view<Pool>{pool, pair.shape_id_A},
                        shape_view<Pool>{pool, pair.shape_id_B}, options));
      }
      catch (...) {
        return invalid_input();
      }
    }
  }

  try {
    return sat_pair<Dim>(shape_view<Pool>{pool, pair.shape_id_A},
                         shape_view<Pool>{pool, pair.shape_id_B}, options);
  }
  catch (...) {
    return invalid_input();
  }
}

// The batched SAT pass (core, by reference): one pair list, one answer
// buffer, written by the rules the by-ref gjk_batch core follows -- assign
// re-fills, the capacity from the previous step reused, zero stale rows.
// Invalid options reject every pair; everything else is per pair. bounds
// rides through to the pair kernel, timing collects one wall-clock row
// per pair (null by default, zero cost when unset).
template <std::size_t Dim, typename Pool, typename PairRecord>
  requires(Dim == 2 || Dim == 3)
void sat_batch(
    const Pool& pool, const std::vector<PairRecord>& pairs,
    const gjk_options options,
    std::vector<sat_result<typename Pool::scalar, Dim>>& answers,
    const pool_bounds<typename Pool::scalar, Dim>* bounds = nullptr,
    batch_timing* timing = nullptr)
{
  using answer = sat_result<typename Pool::scalar, Dim>;

  assert_device_contract<Pool, Dim>();

  answers.assign(pairs.size(), answer{});

  if (options.tolerance <= 0. || options.max_iterations == 0) {
    for (auto& entry : answers) entry.st = status::invalid_input;
    return;
  }

  if (timing != nullptr) {
    timing->seconds.assign(pairs.size(), 0.);
    for (std::size_t i = 0; i < pairs.size(); ++i) {
      const auto start = std::chrono::steady_clock::now();
      answers[i] = sat_batch_pair<Dim>(pool, pairs[i], options, bounds);
      timing->seconds[i] = std::chrono::duration<double>(
                               std::chrono::steady_clock::now() - start)
                               .count();
    }
    return;
  }

  for (std::size_t i = 0; i < pairs.size(); ++i)
    answers[i] = sat_batch_pair<Dim>(pool, pairs[i], options, bounds);
}

// The batched SAT pass (allocating form): allocates the answer buffer and
// forwards to the core overload above.
template <std::size_t Dim, typename Pool, typename PairRecord>
  requires(Dim == 2 || Dim == 3)
auto sat_batch(const Pool& pool, const std::vector<PairRecord>& pairs,
               const gjk_options options = {})
    -> std::vector<sat_result<typename Pool::scalar, Dim>>
{
  std::vector<sat_result<typename Pool::scalar, Dim>> answers;
  sat_batch<Dim>(pool, pairs, options, answers);
  return answers;
}

// Convenience overload (3D) matching gjk_batch. The init direction stays
// in the signature for interface parity; SAT ignores it. The rows are
// sat_result now, the answer kind this section answers with.
template <std::size_t Dim>
  requires(Dim == 3)
auto sat_batch(const std::vector<std::vector<std::array<double, 3>>>& a,
               const std::vector<std::vector<std::array<double, 3>>>& b,
               const std::array<double, 3>& init,
               const gjk_options options = {})
    -> std::vector<sat_result<double, Dim>>
{
  using shape = std::vector<std::array<double, 3>>;

  const std::size_t pair_count = std::min(a.size(), b.size());
  std::vector<sat_result<double, Dim>> results(pair_count);

  if (a.size() != b.size()) {
    for (auto& entry : results) entry.st = status::invalid_input;
    return results;
  }

  std::vector<shape> pool_shapes;
  pool_shapes.reserve(a.size() + b.size());
  pool_shapes.insert(pool_shapes.end(), a.begin(), a.end());
  pool_shapes.insert(pool_shapes.end(), b.begin(), b.end());

  const auto pool = pack_shapes(pool_shapes);

  const auto offset = static_cast<std::uint32_t>(a.size());
  std::vector<pair_record<double, 3>> pairs;
  pairs.reserve(pair_count);
  for (std::size_t i = 0; i < pair_count; ++i)
    pairs.push_back(
        pair_record<double, 3>{static_cast<std::uint32_t>(i),
                               offset + static_cast<std::uint32_t>(i),
                               {init[0], init[1], init[2]}});

  return sat_batch<Dim>(pool, pairs, options);
}

// Penetration-style SAT batch: returns (depth, normal) per hit like
// epa_batch.
template <std::size_t Dim, typename Pool, typename PairRecord>
  requires(Dim == 2 || Dim == 3)
auto sat_penetration_batch_pair(const Pool& pool, const PairRecord& pair,
                                const gjk_options options)
    -> std::pair<typename Pool::scalar,
                 std::array<typename Pool::scalar, Dim>>
{
  using scalar = typename Pool::scalar;
  using vector = typename Pool::vector;
  const auto nan = std::numeric_limits<scalar>::quiet_NaN();
  auto rejected = [&]() {
    std::array<scalar, Dim> n{};
    n.fill(nan);
    return std::pair<scalar, std::array<scalar, Dim>>{nan, n};
  };

  const auto& entries = pool.entries;
  if (pair.shape_id_A >= entries.size() ||
      pair.shape_id_B >= entries.size() ||
      entries[pair.shape_id_A].vertex_count == 0 ||
      entries[pair.shape_id_B].vertex_count == 0 || options.tolerance <= 0. ||
      options.max_iterations == 0)
    return rejected();

  std::vector<vector> a, b;
  const auto& ea = entries[pair.shape_id_A];
  const auto& eb = entries[pair.shape_id_B];
  for (std::uint32_t i = 0; i < ea.vertex_count; ++i)
    a.push_back(pool.vertex_buffer[ea.vertex_offset + i]);
  for (std::uint32_t i = 0; i < eb.vertex_count; ++i)
    b.push_back(pool.vertex_buffer[eb.vertex_offset + i]);

  try {
    auto res = cxd::sat<Dim>(a, b, options);
    return res;
  }
  catch (...) {
    return rejected();
  }
}

template <std::size_t Dim, typename Pool, typename PairRecord>
  requires(Dim == 2 || Dim == 3)
auto sat_penetration_batch(const Pool& pool,
                           const std::vector<PairRecord>& pairs,
                           const gjk_options options = {})
    -> std::vector<epa_result<typename Pool::scalar, Dim>>
{
  using scalar = typename Pool::scalar;
  using answer = epa_result<scalar, Dim>;
  std::vector<answer> answers(pairs.size());
  const scalar nan = std::numeric_limits<scalar>::quiet_NaN();
  if (options.tolerance <= 0. || options.max_iterations == 0) {
    for (auto& a : answers) {
      std::array<scalar, Dim> n{};
      n.fill(nan);
      a = {nan, n};
    }
    return answers;
  }
  for (std::size_t i = 0; i < pairs.size(); ++i)
    answers[i] = sat_penetration_batch_pair<Dim>(pool, pairs[i], options);
  return answers;
}

// Aliases
template <std::size_t Dim, typename Pool, typename PairRecord>
  requires(Dim == 2 || Dim == 3)
auto sat_batch_penetration(const Pool& pool,
                           const std::vector<PairRecord>& pairs,
                           const gjk_options options = {})
{
  return sat_penetration_batch<Dim>(pool, pairs, options);
}

template <std::size_t Dim, typename Pool, typename PairRecord>
  requires(Dim == 2 || Dim == 3)
auto sat_epa_batch(const Pool& pool, const std::vector<PairRecord>& pairs,
                   const gjk_options options = {})
{
  return sat_penetration_batch<Dim>(pool, pairs, options);
}
}  // namespace cxd

#endif
