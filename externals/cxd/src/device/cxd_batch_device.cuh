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
// The device (CUDA) batch backend, phase 5.
//
// Semantics contract: this file mirrors cxd_batch.hpp's per-pair evaluator
// and both penetration branches line for line -- same status values (the
// cxd::status numbering), same epsilon uses (options.tolerance stays double
// and the comparisons promote the same way C++ promotes them on the host),
// same seed handover from the GJK pass to the EPA pass. The parity bar is
// the committed golden corpus: statuses exact, {depth, normal} within the
// phase-4 bands (double 1e-6/1e-4, float 1e-4/1e-3).
//
// Device model (the phase-2 SBO work is the storage model this follows):
// one thread per pair for both passes, no heap, fixed per-thread arenas:
//   - 2D: a ring of up to cxd_dev_ring_cap support positions (corpus worst
//     is 14);
//   - 3D: a vertex set of up to cxd_dev_vert_cap positions and a face list
//     of up to cxd_dev_face_cap faces (corpus worst is 31 vertices / 61
//     faces).
// The host's capacity guards (vertex_count_A + vertex_count_B, twice that
// for the face list) are kept exactly; the fixed caps sit on top, so a pair
// whose polytope would exceed the caps is refused (NaN in place) rather
// than reading out of bounds. The corpus never trips the caps, and a cap
// trip that moved an answer would surface as a failed golden grade.
//
// Deliberately C++17 with no core-header include: the C++20 constrained
// core headers (cxd.hpp / cxd_batch.hpp) cannot be eaten by this nvcc
// front end, and the seam phase 3 defined is buffers plus fixed strides,
// not types. The host behaviours mirrored here live in src/cxd_batch.hpp;
// the corpus lives in src/test/batch_parity_corpus.hpp (shared with the
// core suite, so the golden pin stays one corpus).
//
// Credits carried over from the throwaway spike (gpu-toolchain-spike
// branch): the thread-per-pair kernel shape, the host-compaction two-pass
// pipeline and the port lesson (carry the validated code; diff by
// execution when a variant misbehaves).
#ifndef CXD_BATCH_DEVICE_HPP
#define CXD_BATCH_DEVICE_HPP

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include <cuda_runtime.h>

// ---------------------------------------------------------------------------
// ABI: the transfer units of the phase-3 contract, restated for the device
// (no core types travel; the layouts are the flat ones the contract names).
// ---------------------------------------------------------------------------

// cxd::status's numbering, mirrored: a device row's status byte grades
// directly against the golden statuses through the name mapping in the
// test.
enum cxd_dev_status {
  cxd_dev_separated = 0,
  cxd_dev_enclosed = 1,
  cxd_dev_max_iterations = 2,
  cxd_dev_invalid_input = 3
};

struct cxd_dev_entry {
  std::uint32_t vertex_offset;
  std::uint32_t vertex_count;
};

template <class Scalar, std::size_t Dim>
struct cxd_dev_pair {
  std::uint32_t shape_id_A;
  std::uint32_t shape_id_B;
  Scalar init_dir[Dim];
};

struct cxd_dev_options {
  double tolerance;
  std::size_t max_iterations;
};

// One row per pair: the status byte plus the witness simplex the enclosed
// path hands to the EPA pass (positions only; the pass reads nothing else).
// Nothing but the status byte crosses back to the host; the seed handover
// stays inside the device, exactly where the phase-3 contract's marshaling
// caveat wants it.
template <class Scalar, std::size_t Dim>
struct cxd_dev_gjk_row {
  unsigned char st;
  Scalar seed[(Dim + 1) * Dim];
};

// One row per hit: the penetration answer. Refusal is the host's NaN in
// place: depth NaN and every normal component NaN.
template <class Scalar, std::size_t Dim>
struct cxd_dev_answer {
  Scalar depth;
  Scalar normal[Dim];
};

// Fixed per-thread arena caps (see the header comment): generous over the
// corpus; a cap trip that moved an answer would fail a golden grade loudly
// instead of being absorbed by a tolerance.
inline constexpr std::size_t cxd_dev_ring_cap = 64;   // 2D ring (corpus 14)
inline constexpr std::size_t cxd_dev_vert_cap = 96;   // 3D verts (corpus 31)
inline constexpr std::size_t cxd_dev_face_cap = 192;  // 3D faces (corpus 61)

// ---------------------------------------------------------------------------
// Device math layer: the small vocabulary cxd_batch's expressions use.
// ---------------------------------------------------------------------------

namespace cxd_dev {

__device__ inline double dev_sqrt(double x) { return ::sqrt(x); }
__device__ inline float dev_sqrt(float x) { return ::sqrtf(x); }

__device__ inline double dev_abs(double x) { return ::fabs(x); }
__device__ inline float dev_abs(float x) { return ::fabsf(x); }

__device__ inline double dev_inf_double()
{
  return __longlong_as_double(0x7ff0000000000000LL);
}
__device__ inline float dev_inf_float() { return __int_as_float(0x7f800000); }
__device__ inline double dev_nan_double()
{
  return __longlong_as_double(0x7ff8000000000000ULL);
}
__device__ inline float dev_nan_float() { return __int_as_float(0x7fc00000); }

template <class Scalar>
__device__ inline Scalar dev_inf()
{
  if constexpr (sizeof(Scalar) == 8) return Scalar(dev_inf_double());
  else return Scalar(dev_inf_float());
}

template <class Scalar>
__device__ inline Scalar dev_nan()
{
  if constexpr (sizeof(Scalar) == 8) return dev_nan_double();
  else return dev_nan_float();
}

template <class Scalar, std::size_t Dim>
using vec = std::array<Scalar, Dim>;

// A pooled shape seen from a kernel: base pointer into the flat buffer plus
// vertex count. Vertices of the shape are base[Dim * i + c].
template <class Scalar, std::size_t Dim>
struct shape_window {
  const Scalar* base;
  std::uint32_t count;
};

template <class Scalar, std::size_t Dim>
struct pool_view {
  const Scalar* vertices;
  const cxd_dev_entry* entries;
  std::uint32_t entry_count;

  __device__ inline shape_window<Scalar, Dim> window(
      const std::uint32_t id) const
  {
    const cxd_dev_entry& entry = entries[id];
    return shape_window<Scalar, Dim>{vertices + std::size_t(entry.vertex_offset) * Dim,
                                     entry.vertex_count};
  }
};

template <class Scalar, std::size_t Dim>
__device__ inline Scalar dot(const vec<Scalar, Dim>& a,
                             const vec<Scalar, Dim>& b)
{
  // Mirror the host's util::sum_ct dot (cxd_std.hpp). The host folds
  // ((f(a[I], b[I])) + ...) and gcc 13 expands that fold right-associatively,
  // E0 + (E1 + (E2 + ...)), which is observable when a cancelling pair of
  // equal-magnitude products meets a third smaller term (the collinear
  // corpus rows rely on such ties). Written as an explicit recursion so the
  // association does not depend on the incremental compiler's fold behaviour.
  Scalar prod[Dim];
  for (std::size_t k = 0; k < Dim; ++k) prod[k] = a[k] * b[k];
  if constexpr (Dim == 1)
    return prod[0];
  else {
    Scalar sum = prod[Dim - 1];
    for (std::size_t k = Dim - 1; k-- > 0;) sum = prod[k] + sum;
    return sum;
  }
}

template <class Scalar, std::size_t Dim>
__device__ inline vec<Scalar, Dim> sub(const vec<Scalar, Dim>& a,
                                       const vec<Scalar, Dim>& b)
{
  vec<Scalar, Dim> out{};
  for (std::size_t k = 0; k < Dim; ++k) out[k] = a[k] - b[k];
  return out;
}

template <class Scalar, std::size_t Dim>
__device__ inline vec<Scalar, Dim> add(const vec<Scalar, Dim>& a,
                                       const vec<Scalar, Dim>& b)
{
  vec<Scalar, Dim> out{};
  for (std::size_t k = 0; k < Dim; ++k) out[k] = a[k] + b[k];
  return out;
}

template <class Scalar, std::size_t Dim>
__device__ inline vec<Scalar, Dim> scale(const vec<Scalar, Dim>& a, Scalar f)
{
  vec<Scalar, Dim> out{};
  for (std::size_t k = 0; k < Dim; ++k) out[k] = a[k] * f;
  return out;
}

template <class Scalar, std::size_t Dim>
__device__ inline vec<Scalar, Dim> neg(const vec<Scalar, Dim>& a)
{
  vec<Scalar, Dim> out{};
  for (std::size_t k = 0; k < Dim; ++k) out[k] = -a[k];
  return out;
}

// The 3D cross product, in the operand order the core uses.
template <class Scalar>
__device__ inline vec<Scalar, 3> cross(const vec<Scalar, 3>& a,
                                       const vec<Scalar, 3>& b)
{
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
          a[0] * b[1] - a[1] * b[0]};
}

// ( a x b ) x c, the expansion cxd.hpp's cross3 carries. The host formula is
// component-free, so it ports directly for both 2D and 3D.
template <class Scalar, std::size_t Dim>
__device__ inline vec<Scalar, Dim> cross3(const vec<Scalar, Dim>& a,
                                          const vec<Scalar, Dim>& b,
                                          const vec<Scalar, Dim>& c)
{
  const Scalar cb = dot(c, b);
  const Scalar ca = dot(c, a);
  vec<Scalar, Dim> out{};
  for (std::size_t k = 0; k < Dim; ++k)
    out[k] = Scalar(0) - cb * a[k] + ca * b[k];
  return out;
}

// The host's normalize: 1/sqrt(dot) times the vector, no tolerance handling
// (the callers pre-check their lengths).
template <class Scalar, std::size_t Dim>
__device__ inline vec<Scalar, Dim> normalize(const vec<Scalar, Dim>& v)
{
  const Scalar invnorm = Scalar(1) / dev_sqrt(dot(v, v));
  return scale(v, invnorm);
}

// The host's support scan semantics: strict <, so ties keep the lowest
// vertex index. The fold starts at minus infinity; starting from the first
// point's own dot is the same result, and finite-only in this corpus.
template <class Scalar, std::size_t Dim>
__device__ inline std::uint32_t support_indice(const shape_window<Scalar, Dim>& w,
                                               const vec<Scalar, Dim>& v)
{
  vec<Scalar, Dim> p0{};
  for (std::size_t c = 0; c < Dim; ++c) p0[c] = w.base[c];
  Scalar best = dot(p0, v);
  std::uint32_t indice = 0;
  for (std::uint32_t i = 1; i < w.count; ++i) {
    vec<Scalar, Dim> p{};
    for (std::size_t c = 0; c < Dim; ++c) p[c] = w.base[std::size_t(i) * Dim + c];
    const Scalar cdot = dot(p, v);
    if (best < cdot) {
      best = cdot;
      indice = i;
    }
  }
  return indice;
}

template <class Scalar, std::size_t Dim>
__device__ inline vec<Scalar, Dim> support_point(
    const shape_window<Scalar, Dim>& w, const std::uint32_t i)
{
  vec<Scalar, Dim> p{};
  for (std::size_t c = 0; c < Dim; ++c) p[c] = w.base[std::size_t(i) * Dim + c];
  return p;
}

// One support query over the pool: the exact pool_diff_support behaviour (a
// scan kept lowest-index on ties) reduced to the position the algorithm
// reads -- every comparison in the core goes through support_vector.
template <class Scalar, std::size_t Dim>
__device__ inline vec<Scalar, Dim> diff_support(
    const pool_view<Scalar, Dim>& pool, const std::uint32_t id_a,
    const std::uint32_t id_b, const vec<Scalar, Dim>& direction)
{
  const shape_window<Scalar, Dim> wa = pool.window(id_a);
  const shape_window<Scalar, Dim> wb = pool.window(id_b);
  const std::uint32_t ia = support_indice(wa, direction);
  const std::uint32_t ib = support_indice(wb, neg(direction));
  return sub(support_point(wa, ia), support_point(wb, ib));
}

__device__ inline double dev_clamp01(double x)
{
  return x < 0. ? 0. : (x > 1. ? 1. : x);
}

__device__ inline float dev_clamp01(float x)
{
  return x < 0.f ? 0.f : (x > 1.f ? 1.f : x);
}

}  // namespace cxd_dev

// ---------------------------------------------------------------------------
// The per-pair GJK evaluator, 2D branch: cxd_batch.hpp's gjk_batch_pair<2>
// carried to the device, control flow for control flow.
// ---------------------------------------------------------------------------

namespace cxd_dev_2d {

using cxd_dev::diff_support;
using cxd_dev::dot;
using cxd_dev::sub;


template <class Scalar>
__device__ inline cxd_dev_status put_seed(cxd_dev_gjk_row<Scalar, 2>& row,
                                          const cxd_dev::vec<Scalar, 2>& s0, const cxd_dev::vec<Scalar, 2>& s1,
                                          const cxd_dev::vec<Scalar, 2>& s2)
{
  for (std::size_t c = 0; c < 2; ++c) {
    row.seed[c] = s0[c];
    row.seed[2 + c] = s1[c];
    row.seed[4 + c] = s2[c];
  }
  return cxd_dev_enclosed;
}

// The enclosed statuses carry the simplex's current slots, repeats
// included, exactly as the host does.
template <class Scalar>
__device__ inline cxd_dev_status dev_gjk_pair(
    const cxd_dev::pool_view<Scalar, 2>& pool,
    const cxd_dev_pair<Scalar, 2>& pair, const cxd_dev_options& options,
    cxd_dev_gjk_row<Scalar, 2>& row)
{
  // Validation mirrors gjk_batch_pair: ids outside the pool, an empty
  // shape, a zero initial direction. (The extent checks are structural on
  // the device: the flat buffer has one stride and pair records carry
  // exactly the record's Dim directions.)
  if (pair.shape_id_A >= pool.entry_count ||
      pair.shape_id_B >= pool.entry_count ||
      pool.entries[pair.shape_id_A].vertex_count == 0 ||
      pool.entries[pair.shape_id_B].vertex_count == 0)
    return cxd_dev_invalid_input;

  const cxd_dev::vec<Scalar, 2> init{pair.init_dir[0], pair.init_dir[1]};

  const double tolerance_squared = options.tolerance * options.tolerance;
  if (dot(init, init) <= tolerance_squared) return cxd_dev_invalid_input;

  const auto is_near_zero = [&tolerance_squared](const cxd_dev::vec<Scalar, 2>& v) {
    return dot(v, v) <= tolerance_squared;
  };

  const cxd_dev::vec<Scalar, 2> simplex0 = diff_support(pool, pair.shape_id_A, pair.shape_id_B, init);
  if (dot(simplex0, init) < -options.tolerance) return cxd_dev_separated;

  cxd_dev::vec<Scalar, 2> simplex[3] = {simplex0, simplex0, simplex0};
  std::size_t simplex_size = 1;
  cxd_dev::vec<Scalar, 2> direction = cxd_dev::neg(simplex0);

  for (std::size_t iteration = 0; iteration < options.max_iterations;
       ++iteration) {
    if (is_near_zero(direction))
      return put_seed(row, simplex[0], simplex[1], simplex[2]);

    const cxd_dev::vec<Scalar, 2> candidate =
        diff_support(pool, pair.shape_id_A, pair.shape_id_B, direction);
    if (dot(candidate, direction) < -options.tolerance)
      return cxd_dev_separated;

    for (std::size_t i = 0; i < simplex_size; ++i) {
      const cxd_dev::vec<Scalar, 2> difference = sub(candidate, simplex[i]);
      if (is_near_zero(difference)) return cxd_dev_separated;
    }

    simplex[2] = simplex[1];
    simplex[1] = simplex[0];
    simplex[0] = candidate;
    ++simplex_size;

    const cxd_dev::vec<Scalar, 2> a = simplex[0];
    const cxd_dev::vec<Scalar, 2> a0 = cxd_dev::neg(a);
    if (simplex_size == 2) {
      const cxd_dev::vec<Scalar, 2> ab = sub(simplex[1], a);
      if (dot(ab, a0) <= options.tolerance) {
        simplex_size = 1;
        direction = a0;
        continue;
      }

      direction = cxd_dev::cross3(ab, a0, ab);
      if (is_near_zero(direction))
        return put_seed(row, simplex[0], simplex[1], simplex[2]);
      continue;
    }

    const cxd_dev::vec<Scalar, 2> ab = sub(simplex[1], a);
    const cxd_dev::vec<Scalar, 2> ac = sub(simplex[2], a);
    const cxd_dev::vec<Scalar, 2> ab_direction = cxd_dev::cross3(ac, ab, ab);
    if (dot(ab_direction, a0) > options.tolerance) {
      simplex_size = 2;
      direction = ab_direction;
      continue;
    }

    const cxd_dev::vec<Scalar, 2> ac_direction = cxd_dev::cross3(ab, ac, ac);
    if (dot(ac_direction, a0) > options.tolerance) {
      simplex[1] = simplex[2];
      simplex_size = 2;
      direction = ac_direction;
      continue;
    }

    return put_seed(row, simplex[0], simplex[1], simplex[2]);
  }

  return cxd_dev_max_iterations;
}

}  // namespace cxd_dev_2d

// ---------------------------------------------------------------------------
// The per-pair GJK evaluator, 3D branch: gjk_batch_pair<3> carried to the
// device, control flow for control flow -- the loop-state stall detector,
// the duplicate-support handling, the face drops and the containment
// decisions all in the host's order.
// ---------------------------------------------------------------------------

namespace cxd_dev_3d {

using cxd_dev::cross;
using cxd_dev::diff_support;
using cxd_dev::dot;
using cxd_dev::sub;
using cxd_dev::add;
using cxd_dev::scale;
using cxd_dev::neg;
using cxd_dev::dev_abs;
using cxd_dev::dev_sqrt;
using cxd_dev::dev_inf;
using cxd_dev::dev_clamp01;

// Swap-remove dedup starting at index 1, for the face-drop reductions:
// exactly the host's inner loops there.
template <class Scalar>
__device__ inline std::size_t dedup_from_1(cxd_dev::vec<Scalar, 3>* sim, std::size_t size,
                                           const double tolerance_squared);

// Swap-remove dedup over the live simplex (positions): exactly the host's
// dedup_simplex lambda, including the trailing slots it leaves behind. The
// caller updates its size from the returned count, the way the captured
// --simplex_size does on the host.
template <class Scalar>
__device__ inline std::size_t dedup_simplex(cxd_dev::vec<Scalar, 3>* sim, std::size_t size,
                                            const double tolerance_squared)
{
  for (std::size_t i = 0; i < size; ++i) {
    for (std::size_t j = i + 1; j < size;) {
      const cxd_dev::vec<Scalar, 3> difference = sub(sim[i], sim[j]);
      if (dot(difference, difference) <= tolerance_squared) {
        sim[j] = sim[size - 1];
        --size;
      }
      else {
        ++j;
      }
    }
  }
  return size;
}

// The host's compact_simplex: dedups the LIVE simplex first (the captured
// --simplex_size side effect), then copies into a padded fresh array. The
// returned count is the deduped size; the caller updates its own variable.
template <class Scalar>
__device__ inline std::size_t compact_simplex(cxd_dev::vec<Scalar, 3>* sim, std::size_t size,
                                              const double tolerance_squared,
                                              cxd_dev::vec<Scalar, 3>* out)
{
  const std::size_t sz = dedup_simplex(sim, size, tolerance_squared);
  for (std::size_t i = 0; i < sz; ++i) out[i] = sim[i];
  for (std::size_t i = sz; i < 4; ++i) out[i] = out[sz - 1];
  return sz;
}

template <class Scalar>
__device__ inline cxd_dev_status put_seed(cxd_dev_gjk_row<Scalar, 3>& row,
                                          const cxd_dev::vec<Scalar, 3>& p0, const cxd_dev::vec<Scalar, 3>& p1,
                                          const cxd_dev::vec<Scalar, 3>& p2, const cxd_dev::vec<Scalar, 3>& p3)
{
  for (std::size_t c = 0; c < 3; ++c) {
    row.seed[c] = p0[c];
    row.seed[3 + c] = p1[c];
    row.seed[6 + c] = p2[c];
    row.seed[9 + c] = p3[c];
  }
  return cxd_dev_enclosed;
}

// The per-pair evaluator: gjk_batch_pair<3> carried to the device, control
// flow for control flow -- the loop-state stall detector, the
// duplicate-support handling, the face drops and the containment decisions
// in the host's order.
template <class Scalar>
__device__ inline cxd_dev_status dev_gjk_pair(
    const cxd_dev::pool_view<Scalar, 3>& pool,
    const cxd_dev_pair<Scalar, 3>& pair, const cxd_dev_options& options,
    cxd_dev_gjk_row<Scalar, 3>& row)
{
  // Validation mirrors gjk_batch_pair: ids outside the pool, an empty
  // shape, a zero initial direction.
  if (pair.shape_id_A >= pool.entry_count ||
      pair.shape_id_B >= pool.entry_count ||
      pool.entries[pair.shape_id_A].vertex_count == 0 ||
      pool.entries[pair.shape_id_B].vertex_count == 0)
    return cxd_dev_invalid_input;

  const cxd_dev::vec<Scalar, 3> init{pair.init_dir[0], pair.init_dir[1], pair.init_dir[2]};

  const double tolerance_squared = options.tolerance * options.tolerance;
  if (dot(init, init) <= tolerance_squared) return cxd_dev_invalid_input;

  const auto is_near_zero = [&tolerance_squared](const cxd_dev::vec<Scalar, 3>& v) {
    return dot(v, v) <= tolerance_squared;
  };

  const cxd_dev::vec<Scalar, 3> simplex0 = diff_support(pool, pair.shape_id_A, pair.shape_id_B, init);
  if (dot(simplex0, init) < -options.tolerance) return cxd_dev_separated;

  cxd_dev::vec<Scalar, 3> simplex[4] = {simplex0, simplex0, simplex0, simplex0};
  std::size_t simplex_size = 1;
  cxd_dev::vec<Scalar, 3> direction = cxd_dev::neg(simplex0);

  cxd_dev::vec<Scalar, 3> last_loop_dir{};
  std::size_t last_loop_size = 0;
  bool has_loop_state = false;
  std::size_t same_loop_state = 0;

  for (std::size_t iteration = 0; iteration < options.max_iterations;
       ++iteration) {
    // Loop-state stall detector, exactly as the host: when the search
    // revisits an identical (simplex size, direction) pair, answer from the
    // converged simplex.
    const Scalar dir_dot_last = dot(direction, last_loop_dir);
    const Scalar dir_norm2 = dot(direction, direction);
    const Scalar last_norm2 = dot(last_loop_dir, last_loop_dir);
    const bool same_axis =
        has_loop_state &&
        dir_dot_last * dir_dot_last >= 0.9 * dir_norm2 * last_norm2;
    if (same_axis && simplex_size == last_loop_size)
      ++same_loop_state;
    else
      same_loop_state = 0;
    has_loop_state = true;
    last_loop_dir = direction;
    last_loop_size = simplex_size;

    if (same_loop_state >= 1 && simplex_size == 3) {
      // Containment judge over the converged triangle, then either the
      // face-normal support separation certificate or the axis probes.
      const cxd_dev::vec<Scalar, 3> t0 = simplex[0];
      const cxd_dev::vec<Scalar, 3> t1 = simplex[1];
      const cxd_dev::vec<Scalar, 3> t2 = simplex[2];
      const cxd_dev::vec<Scalar, 3> t01 = sub(t1, t0);
      const cxd_dev::vec<Scalar, 3> t02 = sub(t2, t0);
      const cxd_dev::vec<Scalar, 3> t12 = sub(t2, t1);
      const cxd_dev::vec<Scalar, 3> tn = cross(t01, t02);
      const Scalar scale = dot(t01, t01) + dot(t02, t02) + dot(tn, tn);
      const cxd_dev::vec<Scalar, 3> t0n = cxd_dev::neg(t0);
      const bool origin_in_triangle =
          cxd_dev::dev_abs(dot(tn, t0n)) <= options.tolerance * scale &&
          dot(cross(t01, tn), t0n) <= options.tolerance * scale &&
          dot(cross(tn, t02), t0n) <= options.tolerance * scale &&
          dot(cross(t12, tn), cxd_dev::neg(t1)) <= options.tolerance * scale;
      if (origin_in_triangle) {
        cxd_dev::vec<Scalar, 3> compacted[4];
        compact_simplex(simplex, simplex_size, tolerance_squared, compacted);
        return put_seed(row, compacted[0], compacted[1], compacted[2],
                        compacted[3]);
      }
      cxd_dev::vec<Scalar, 3> toward_origin = tn;
      if (dot(tn, t0) > Scalar(0)) toward_origin = cxd_dev::neg(tn);
      {
        const cxd_dev::vec<Scalar, 3> sep_witness =
            diff_support(pool, pair.shape_id_A, pair.shape_id_B, toward_origin);
        const Scalar sep_proj = dot(sep_witness, toward_origin);
        cxd_dev::vec<Scalar, 3> probe = toward_origin;
        bool separated = sep_proj < -options.tolerance;
        if (!separated) {
          const Scalar probe_candidates[6][3] = {{1, 0, 0}, {-1, 0, 0},
                                                 {0, 1, 0}, {0, -1, 0},
                                                 {0, 0, 1}, {0, 0, -1}};
          for (std::size_t pi = 0; pi < 6 && !separated; ++pi) {
            probe = cxd_dev::vec<Scalar, 3>{Scalar(probe_candidates[pi][0]),
                        Scalar(probe_candidates[pi][1]),
                        Scalar(probe_candidates[pi][2])};
            const cxd_dev::vec<Scalar, 3> probe_w = diff_support(pool, pair.shape_id_A,
                                             pair.shape_id_B, probe);
            separated = dot(probe_w, probe) < -options.tolerance;
          }
        }
        cxd_dev::vec<Scalar, 3> compacted[4];
        compact_simplex(simplex, simplex_size, tolerance_squared, compacted);
        if (separated) return cxd_dev_separated;
        return put_seed(row, compacted[0], compacted[1], compacted[2],
                        sep_witness);
      }
    }
    if (is_near_zero(direction)) {
      simplex_size = dedup_simplex(simplex, simplex_size, tolerance_squared);
      if (simplex_size >= 3) {
        const cxd_dev::vec<Scalar, 3> a = simplex[0];
        const cxd_dev::vec<Scalar, 3> ab = sub(simplex[1], a);
        const cxd_dev::vec<Scalar, 3> ac = sub(simplex[2], a);
        const cxd_dev::vec<Scalar, 3> abc = cross(ab, ac);
        if (!is_near_zero(abc)) {
          direction = dot(abc, a) > Scalar(0) ? cxd_dev::neg(abc) : abc;
          continue;
        }
      }
      if (simplex_size >= 2) {
        const cxd_dev::vec<Scalar, 3> a = simplex[0];
        const cxd_dev::vec<Scalar, 3> ab = sub(simplex[1], a);
        direction = cross(ab, cxd_dev::vec<Scalar, 3>{Scalar(1), Scalar(0), Scalar(0)});
        if (is_near_zero(direction))
          direction = cross(ab, cxd_dev::vec<Scalar, 3>{Scalar(0), Scalar(1), Scalar(0)});
        if (dot(direction, a) < Scalar(0)) direction = cxd_dev::neg(direction);
        continue;
      }
      direction = cxd_dev::neg(simplex[0]);
      continue;
    }

    const cxd_dev::vec<Scalar, 3> candidate =
        diff_support(pool, pair.shape_id_A, pair.shape_id_B, direction);
    if (dot(candidate, direction) < -options.tolerance)
      return cxd_dev_separated;

    bool duplicate = false;
    for (std::size_t i = 0; i < simplex_size; ++i) {
      const cxd_dev::vec<Scalar, 3> difference = sub(candidate, simplex[i]);
      if (is_near_zero(difference)) {
        duplicate = true;
        break;
      }
    }
    if (duplicate) {
      const std::size_t final_size =
          dedup_simplex(simplex, simplex_size, tolerance_squared);
      simplex_size = final_size;
      if (final_size < 3) {
        simplex_size = 1;
        direction = cxd_dev::neg(simplex[0]);
        continue;
      }
      if (simplex_size == 4) {
        // Containment-aware stall handling: only an origin-enclosing
        // tetrahedron witnesses an intersection.
        const cxd_dev::vec<Scalar, 3> a = simplex[0];
        const cxd_dev::vec<Scalar, 3> b = simplex[1];
        const cxd_dev::vec<Scalar, 3> c = simplex[2];
        const cxd_dev::vec<Scalar, 3> d = simplex[3];
        const cxd_dev::vec<Scalar, 3> ab = sub(b, a);
        const cxd_dev::vec<Scalar, 3> ac = sub(c, a);
        const cxd_dev::vec<Scalar, 3> ad = sub(d, a);
        const cxd_dev::vec<Scalar, 3> bc = sub(c, b);
        const cxd_dev::vec<Scalar, 3> bd = sub(d, b);
        const bool encloses =
            dot(cross(ab, ac), ad) * dot(cross(ab, ac), cxd_dev::neg(a)) >
                Scalar(0) &&
            dot(cross(ab, ad), ac) * dot(cross(ab, ad), cxd_dev::neg(a)) >
                Scalar(0) &&
            dot(cross(ac, ad), ab) * dot(cross(ac, ad), cxd_dev::neg(a)) >
                Scalar(0) &&
            dot(cross(bd, bc), sub(a, b)) * dot(cross(bd, bc), cxd_dev::neg(b)) >
                Scalar(0);
        if (encloses) {
          cxd_dev::vec<Scalar, 3> compacted[4];
          compact_simplex(simplex, simplex_size, tolerance_squared, compacted);
          return put_seed(row, compacted[0], compacted[1], compacted[2],
                          diff_support(pool, pair.shape_id_A, pair.shape_id_B,
                                       direction));
        }
        simplex_size = 3;
      }
      // Stalled three-point simplex: prism containment or the stall third
      // slot. (The host computes a fresh fourth support along the converged
      // face normal into a variable it then never reads -- the returned
      // seed is the padded compaction. This port keeps the returned seed
      // identical, so the query stays out; a comment carries the note.)
      {
        const cxd_dev::vec<Scalar, 3> a = simplex[0];
        const cxd_dev::vec<Scalar, 3> a0 = cxd_dev::neg(a);
        const cxd_dev::vec<Scalar, 3> ab = sub(simplex[1], a);
        const cxd_dev::vec<Scalar, 3> ac = sub(simplex[2], a);
        const cxd_dev::vec<Scalar, 3> abc = cross(ab, ac);
        const bool origin_in_prism =
            dot(cross(ab, abc), a0) <= options.tolerance &&
            dot(cross(abc, ac), a0) <= options.tolerance;
        cxd_dev::vec<Scalar, 3> stall_compacted[4];
        compact_simplex(simplex, simplex_size, tolerance_squared,
                        stall_compacted);
        return put_seed(row, stall_compacted[0], stall_compacted[1],
                        stall_compacted[2], stall_compacted[3]);
      }
    }

    if (simplex_size < 4) {
      for (std::size_t i = simplex_size; i > 0; --i) {
        simplex[i] = simplex[i - 1];
      }
      simplex[0] = candidate;
      ++simplex_size;
    }

    simplex_size = dedup_simplex(simplex, simplex_size, tolerance_squared);

    const cxd_dev::vec<Scalar, 3> a = simplex[0];
    const cxd_dev::vec<Scalar, 3> a0 = cxd_dev::neg(a);

    if (simplex_size == 2) {
      const cxd_dev::vec<Scalar, 3> ab = sub(simplex[1], a);
      const Scalar ab2 = dot(ab, ab);
      const Scalar ab_ao = dot(ab, a0);
      if (ab_ao <= options.tolerance) {
        simplex_size = 1;
        direction = a0;
        continue;
      }
      if (ab_ao >= ab2 - options.tolerance) {
        simplex[0] = simplex[1];
        simplex_size = 1;
        direction = cxd_dev::neg(simplex[0]);
        continue;
      }
      direction = sub(a0, cxd_dev::scale(ab, ab_ao / ab2));
      if (is_near_zero(direction)) {
        cxd_dev::vec<Scalar, 3> perp = cross(ab, cxd_dev::vec<Scalar, 3>{Scalar(1), Scalar(0), Scalar(0)});
        if (is_near_zero(perp)) perp = cross(ab, cxd_dev::vec<Scalar, 3>{Scalar(0), Scalar(1), Scalar(0)});
        direction = perp;
      }
      continue;
    }

    if (simplex_size == 3) {
      const cxd_dev::vec<Scalar, 3> ab = sub(simplex[1], a);
      const cxd_dev::vec<Scalar, 3> ac = sub(simplex[2], a);
      const cxd_dev::vec<Scalar, 3> abc = cross(ab, ac);
      const Scalar abc_ao = dot(abc, a0);

      const cxd_dev::vec<Scalar, 3> ab_cross_abc = cross(ab, abc);
      const cxd_dev::vec<Scalar, 3> ac_cross_abc = cross(abc, ac);

      if (dot(ab_cross_abc, a0) > options.tolerance) {
        simplex_size = 2;
        const cxd_dev::vec<Scalar, 3> ab2 = sub(simplex[1], simplex[0]);
        const Scalar ab2_len = dot(ab2, ab2);
        const Scalar ab_ao2 = dot(ab2, cxd_dev::neg(simplex[0]));
        if (ab_ao2 <= options.tolerance) {
          simplex_size = 1;
          direction = cxd_dev::neg(simplex[0]);
          continue;
        }
        direction = sub(cxd_dev::neg(simplex[0]),
                        cxd_dev::scale(ab2, ab_ao2 / ab2_len));
        continue;
      }

      if (dot(ac_cross_abc, a0) > options.tolerance) {
        simplex[1] = simplex[2];
        simplex_size = 2;
        const cxd_dev::vec<Scalar, 3> ac2 = sub(simplex[1], simplex[0]);
        const Scalar ac2_len = dot(ac2, ac2);
        const Scalar ac_ao2 = dot(ac2, cxd_dev::neg(simplex[0]));
        if (ac_ao2 <= options.tolerance) {
          simplex_size = 1;
          direction = cxd_dev::neg(simplex[0]);
          continue;
        }
        direction = sub(cxd_dev::neg(simplex[0]),
                        cxd_dev::scale(ac2, ac_ao2 / ac2_len));
        continue;
      }

      if (abc_ao > options.tolerance) {
        simplex_size = 3;
        direction = cxd_dev::neg(abc);
      }
      else
        direction = abc;

      continue;
    }

    if (simplex_size == 4) {
      const cxd_dev::vec<Scalar, 3> ab = sub(simplex[1], a);
      const cxd_dev::vec<Scalar, 3> ac = sub(simplex[2], a);
      const cxd_dev::vec<Scalar, 3> ad = sub(simplex[3], a);
      const cxd_dev::vec<Scalar, 3> abc = cross(ab, ac);
      const cxd_dev::vec<Scalar, 3> abd = cross(ab, ad);
      const cxd_dev::vec<Scalar, 3> acd = cross(ad, ac);

      const Scalar abc_ao = dot(abc, a0) * dot(abc, ad);
      const Scalar abd_ao = dot(abd, a0) * dot(abd, ac);
      const Scalar acd_ao = dot(acd, a0) * dot(acd, ab);
      if (abc_ao < -options.tolerance * cxd_dev::dev_abs(dot(abc, ad))) {
        simplex_size = 3;
        direction = cxd_dev::neg(abc);
        simplex_size = dedup_from_1(simplex, simplex_size, tolerance_squared);
        if (simplex_size < 3) {
          simplex_size = 1;
          direction = cxd_dev::neg(simplex[0]);
        }
        continue;
      }

      if (abd_ao < -options.tolerance * cxd_dev::dev_abs(dot(abd, ac))) {
        simplex[2] = simplex[3];
        simplex_size = 3;
        direction = cxd_dev::neg(abd);
        simplex_size = dedup_from_1(simplex, simplex_size, tolerance_squared);
        if (simplex_size < 3) {
          simplex_size = 1;
          direction = cxd_dev::neg(simplex[0]);
        }
        continue;
      }

      if (acd_ao < -options.tolerance * cxd_dev::dev_abs(dot(acd, ab))) {
        simplex[1] = simplex[2];
        simplex[2] = simplex[3];
        simplex_size = 3;
        direction = cxd_dev::neg(acd);
        simplex_size = dedup_from_1(simplex, simplex_size, tolerance_squared);
        if (simplex_size < 3) {
          simplex_size = 1;
          direction = cxd_dev::neg(simplex[0]);
        }
        continue;
      }

      // The fourth face, BCD, was never tested by the three above.
      {
        const cxd_dev::vec<Scalar, 3> bc4 = sub(simplex[2], simplex[1]);
        const cxd_dev::vec<Scalar, 3> bd4 = sub(simplex[3], simplex[1]);
        const cxd_dev::vec<Scalar, 3> bcd4 = cross(bd4, bc4);
        const Scalar bcd_bo =
            dot(bcd4, cxd_dev::neg(simplex[1])) *
            dot(bcd4, sub(simplex[0], simplex[1]));
        if (bcd_bo < -options.tolerance *
                         cxd_dev::dev_abs(dot(bcd4, sub(simplex[0], simplex[1])))) {
          simplex[0] = simplex[1];
          simplex[1] = simplex[2];
          simplex[2] = simplex[3];
          simplex_size = 3;
          direction = cxd_dev::neg(bcd4);
          simplex_size = dedup_from_1(simplex, simplex_size, tolerance_squared);
          if (simplex_size < 3) {
            simplex_size = 1;
            direction = cxd_dev::neg(simplex[0]);
          }
          continue;
        }
      }

      cxd_dev::vec<Scalar, 3> compacted2[4];
      compact_simplex(simplex, simplex_size, tolerance_squared, compacted2);
      const cxd_dev::vec<Scalar, 3> a2 = compacted2[0];
      const cxd_dev::vec<Scalar, 3> b2 = compacted2[1];
      const cxd_dev::vec<Scalar, 3> c2 = compacted2[2];
      const cxd_dev::vec<Scalar, 3> d2 = compacted2[3];
      const cxd_dev::vec<Scalar, 3> ab2 = sub(b2, a2);
      const cxd_dev::vec<Scalar, 3> ac2 = sub(c2, a2);
      const cxd_dev::vec<Scalar, 3> ad2 = sub(d2, a2);
      const Scalar volume_sq = dot(cross(ab2, ac2), ad2);
      if (volume_sq <= options.tolerance * options.tolerance) {
        simplex_size = 3;
        direction = cxd_dev::neg(cross(ab2, ac2));
        continue;
      }

      return put_seed(row, a2, b2, c2, d2);
    }

    cxd_dev::vec<Scalar, 3> compacted3[4];
    compact_simplex(simplex, simplex_size, tolerance_squared, compacted3);
    return put_seed(row, compacted3[0], compacted3[1], compacted3[2],
                    diff_support(pool, pair.shape_id_A, pair.shape_id_B,
                                 direction));
  }

  return cxd_dev_max_iterations;
}

// Swap-remove dedup starting at index 1 (indices 0..3), for the
// face-drop reductions: exactly the host's inner loops there.
template <class Scalar>
__device__ inline std::size_t dedup_from_1(cxd_dev::vec<Scalar, 3>* sim, std::size_t size,
                                           const double tolerance_squared)
{
  for (std::size_t i = 1; i < size; ++i) {
    for (std::size_t j = i + 1; j < size;) {
      const cxd_dev::vec<Scalar, 3> difference = sub(sim[i], sim[j]);
      if (dot(difference, difference) <= tolerance_squared) {
        sim[j] = sim[size - 1];
        --size;
      }
      else {
        ++j;
      }
    }
  }
  return size;
}

}  // namespace cxd_dev_3d

// ---------------------------------------------------------------------------
// The penetration pass, 2D branch: epa_penetration_2d on the device.
// ---------------------------------------------------------------------------

namespace cxd_dev {

// The batch's unanswered rows, in a shared spot: depth NaN, every normal
// component NaN, exactly what the host puts in place of an answer.
template <class Scalar, std::size_t Dim>
__device__ inline cxd_dev_answer<Scalar, Dim> refused()
{
  cxd_dev_answer<Scalar, Dim> out{};
  out.depth = dev_nan<Scalar>();
  for (std::size_t k = 0; k < Dim; ++k) out.normal[k] = dev_nan<Scalar>();
  return out;
}

}  // namespace cxd_dev

namespace cxd_dev_2d {

using cxd_dev::add;
using cxd_dev::scale;

template <class Scalar>
__device__ inline cxd_dev_answer<Scalar, 2> dev_epa_pair(
    const cxd_dev::pool_view<Scalar, 2>& pool,
    const cxd_dev_pair<Scalar, 2>& pair, const cxd_dev::vec<Scalar, 2> seeds[3],
    const cxd_dev_options& options)
{
  // Validation mirrors epa_penetration_2d: an option pair epa would reject,
  // or a shape outside the pool.
  if (options.tolerance <= 0. || options.max_iterations == 0)
    return cxd_dev::refused<Scalar, 2>();
  if (pair.shape_id_A >= pool.entry_count ||
      pair.shape_id_B >= pool.entry_count ||
      pool.entries[pair.shape_id_A].vertex_count == 0 ||
      pool.entries[pair.shape_id_B].vertex_count == 0)
    return cxd_dev::refused<Scalar, 2>();

  // Non-const so the degenerate bootstrap and canonicalization can reassign.
  cxd_dev::vec<Scalar, 2> a = seeds[0];
  cxd_dev::vec<Scalar, 2> b = seeds[1];
  cxd_dev::vec<Scalar, 2> c = seeds[2];
  const Scalar area =
      (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
  if (cxd_dev::dev_abs(area) <= options.tolerance) {
    // Degenerate witness bootstrap, mirroring epa_penetration_2d's recovery.
    // The 2D GJK search direction collapses when the origin lies on the
    // segment between two supports, so the seed triangle has zero area.
    cxd_dev::vec<Scalar, 2> distinct[3];
    std::size_t distinct_count = 0;
    const cxd_dev::vec<Scalar, 2> seed_arr[3] = {a, b, c};
    for (std::size_t i = 0; i < 3; ++i) {
      bool found = false;
      for (std::size_t j = 0; j < distinct_count; ++j) {
        const cxd_dev::vec<Scalar, 2> d = sub(seed_arr[i], distinct[j]);
        if (dot(d, d) <= options.tolerance * options.tolerance) {
          found = true;
          break;
        }
      }
      if (!found) distinct[distinct_count++] = seed_arr[i];
    }

    if (distinct_count == 1) {
      // One support: the difference is a point at the origin (touching
      // contact with no direction).
      const cxd_dev::vec<Scalar, 2> p = distinct[0];
      const Scalar distance = cxd_dev::dev_sqrt(dot(p, p));
      cxd_dev_answer<Scalar, 2> out{};
      out.depth = distance;
      if (distance > options.tolerance) {
        for (std::size_t k = 0; k < 2; ++k) out.normal[k] = p[k] / distance;
      } else {
        out.normal[0] = Scalar(0);
        out.normal[1] = Scalar(0);
      }
      return out;
    }

    // distinct_count == 2: span a line, probe off it for a third support.
    cxd_dev::vec<Scalar, 2> p0 = distinct[0];
    cxd_dev::vec<Scalar, 2> p1 = distinct[1];
    const cxd_dev::vec<Scalar, 2> line = sub(p1, p0);
    const Scalar line_length = cxd_dev::dev_sqrt(dot(line, line));
    bool bootstrapped = false;
    if (line_length > options.tolerance) {
      const cxd_dev::vec<Scalar, 2> perp{-line[1], line[0]};
      const cxd_dev::vec<Scalar, 2> oblique_a{line[0] + perp[0],
                                             line[1] + perp[1]};
      const cxd_dev::vec<Scalar, 2> oblique_b{line[0] - perp[0],
                                             line[1] - perp[1]};
      const cxd_dev::vec<Scalar, 2> neg_perp{-perp[0], -perp[1]};
      const cxd_dev::vec<Scalar, 2> probes[4] = {
          perp, neg_perp, oblique_a, oblique_b};
      for (std::size_t pi = 0; pi < 4 && !bootstrapped; ++pi) {
        const cxd_dev::vec<Scalar, 2>& probe = probes[pi];
        if (dot(probe, probe) <=
            options.tolerance * options.tolerance)
          continue;
        const cxd_dev::vec<Scalar, 2> candidate =
            cxd_dev::diff_support(pool, pair.shape_id_A,
                                  pair.shape_id_B, probe);
        const cxd_dev::vec<Scalar, 2> off = sub(candidate, p0);
        const Scalar perpendicular_distance =
            cxd_dev::dev_abs(line[0] * off[1] - line[1] * off[0]) /
            line_length;
        if (perpendicular_distance > options.tolerance) {
          a = p0;
          b = p1;
          c = candidate;
          bootstrapped = true;
        }
      }
    }

    if (!bootstrapped) {
      // The difference polytope really is a segment: closest point on
      // it to the origin.
      const cxd_dev::vec<Scalar, 2> d = line;
      const Scalar dd = dot(d, d);
      Scalar t = Scalar(0);
      if (dd > Scalar(0)) {
        t = -dot(p0, d) / dd;
        if (t < Scalar(0)) t = Scalar(0);
        if (t > Scalar(1)) t = Scalar(1);
      }
      const cxd_dev::vec<Scalar, 2> closest =
          add(p0, scale(d, t));
      const Scalar distance = cxd_dev::dev_sqrt(dot(closest, closest));
      cxd_dev::vec<Scalar, 2> normal{};
      if (distance > options.tolerance) {
        for (std::size_t k = 0; k < 2; ++k)
          normal[k] = closest[k] / distance;
      } else if (dd > Scalar(0)) {
        const Scalar dlen = cxd_dev::dev_sqrt(dd);
        for (std::size_t k = 0; k < 2; ++k)
          normal[k] = d[k] / dlen;
      } else {
        normal[0] = Scalar(0);
        normal[1] = Scalar(0);
      }
      cxd_dev_answer<Scalar, 2> out{};
      out.depth = distance;
      out.normal[0] = normal[0];
      out.normal[1] = normal[1];
      return out;
    }
  }

  // Orientation canonicalization: canonicalize to counterclockwise so the
  // exterior normal points away from the origin regardless of the witness
  // order GJK handed over. After this, cw is always false.
  const Scalar cross =
      b[0] * c[1] + a[0] * b[1] + a[1] * c[0] -
      (a[1] * b[0] + b[1] * c[0] + a[0] * c[1]);
  if (cross < Scalar(0)) {
    const cxd_dev::vec<Scalar, 2> tmp = b;
    b = c;
    c = tmp;
  }
  const bool cw = false;

  // The polytope can hold one vertex per point of the two shapes'
  // difference; the fixed ring sits on top of that bound.
  const std::size_t max_polytope =
      pool.entries[pair.shape_id_A].vertex_count +
      pool.entries[pair.shape_id_B].vertex_count;

  cxd_dev::vec<Scalar, 2> ring[cxd_dev_ring_cap];
  std::size_t ring_size = 0;
  ring[ring_size++] = a;
  ring[ring_size++] = b;
  ring[ring_size++] = c;

  for (std::size_t iteration = 0; iteration < options.max_iterations;
       ++iteration) {
    Scalar closest_distance = cxd_dev::dev_inf<Scalar>();
    cxd_dev::vec<Scalar, 2> closest_normal{};
    std::size_t closest_edge = 0;

    for (std::size_t i = 0; i < ring_size; ++i) {
      const cxd_dev::vec<Scalar, 2> p1 = ring[i];
      const cxd_dev::vec<Scalar, 2> p2 = ring[(i + 1) % ring_size];
      const cxd_dev::vec<Scalar, 2> edge = sub(p2, p1);
      if (dot(edge, edge) <= options.tolerance * options.tolerance)
        return cxd_dev::refused<Scalar, 2>();

      // The host's exterior_normal: !cw ? {edge[1], -edge[0]} : the
      // opposite, normalized.
      const cxd_dev::vec<Scalar, 2> normal = cxd_dev::normalize(cw
          ? cxd_dev::vec<Scalar, 2>{-edge[1], edge[0]}
          : cxd_dev::vec<Scalar, 2>{edge[1], -edge[0]});
      const Scalar distance = cxd_dev::dev_abs(dot(normal, p1));
      if (distance < closest_distance) {
        closest_distance = distance;
        closest_normal = normal;
        closest_edge = i;
      }
    }

    const cxd_dev::vec<Scalar, 2> candidate =
        cxd_dev::diff_support(pool, pair.shape_id_A, pair.shape_id_B,
                              closest_normal);
    const Scalar candidate_distance = dot(candidate, closest_normal);
    if (candidate_distance - closest_distance <= options.tolerance) {
      cxd_dev_answer<Scalar, 2> out{};
      out.depth = candidate_distance;
      out.normal[0] = closest_normal[0];
      out.normal[1] = closest_normal[1];
      return out;
    }

    for (std::size_t i = 0; i < ring_size; ++i) {
      const cxd_dev::vec<Scalar, 2> difference = sub(candidate, ring[i]);
      if (dot(difference, difference) <=
          options.tolerance * options.tolerance) {
        cxd_dev_answer<Scalar, 2> out{};
        out.depth = closest_distance;
        out.normal[0] = closest_normal[0];
        out.normal[1] = closest_normal[1];
        return out;
      }
    }

    if (ring_size == max_polytope || ring_size >= cxd_dev_ring_cap)
      return cxd_dev::refused<Scalar, 2>();

    // The host's insert_at: shift the tail up by one.
    for (std::size_t i = ring_size; i > closest_edge + 1; --i)
      ring[i] = ring[i - 1];
    ring[closest_edge + 1] = candidate;
    ++ring_size;
  }

  return cxd_dev::refused<Scalar, 2>();
}

}  // namespace cxd_dev_2d

// ---------------------------------------------------------------------------
// The penetration pass, 3D branch: epa_penetration_3d on the device -- the
// degenerate-witness bootstrap, the portal refinement, the flat-shape
// projection, the witness repair and the face walk.
// ---------------------------------------------------------------------------

namespace cxd_dev {

// The analyze() result of the degenerate bootstrap.
template <class Scalar>
struct degeneracy_info {
  bool have_line = false;
  vec<Scalar, 3> line_dir{};
  vec<Scalar, 3> line_point{};
  bool have_plane = false;
  vec<Scalar, 3> normal{};
  vec<Scalar, 3> plane_point{};
  bool have_volume = false;
  vec<Scalar, 3> volume_point{};
};

// The project_onto_points() result.
template <class Scalar>
struct epa_projection {
  Scalar depth;
  vec<Scalar, 3> normal;
};

}  // namespace cxd_dev

namespace cxd_dev_3d {

using cxd_dev::cross;
using cxd_dev::diff_support;
using cxd_dev::dot;
using cxd_dev::sub;
using cxd_dev::add;
using cxd_dev::scale;
using cxd_dev::neg;
using cxd_dev::dev_abs;
using cxd_dev::dev_sqrt;
using cxd_dev::dev_inf;
using cxd_dev::dev_clamp01;

template <class Scalar>
struct dev_face {
  std::size_t idx[3];
  cxd_dev::vec<Scalar, 3> normal{};
  Scalar distance{};
  bool removed = false;
};

template <class Scalar>
__device__ inline cxd_dev::degeneracy_info<Scalar> dev_analyze(
    const cxd_dev::vec<Scalar, 3>* pts, const std::size_t n, const double tolerance)
{
  cxd_dev::degeneracy_info<Scalar> info;
  if (n < 2) return info;
  info.have_line = true;
  info.line_point = pts[1];
  const cxd_dev::vec<Scalar, 3> e1 = sub(pts[1], pts[0]);
  const Scalar e1_len2 = dot(e1, e1);
  // The host divides every component; not a multiply.
  for (std::size_t c = 0; c < 3; ++c) info.line_dir[c] = e1[c] / dev_sqrt(e1_len2);
  for (std::size_t k = 2; k < n; ++k) {
    const cxd_dev::vec<Scalar, 3> ek = sub(pts[k], pts[0]);
    if (!info.have_plane) {
      const cxd_dev::vec<Scalar, 3> n_raw = cross(info.line_dir, ek);
      const Scalar n_len2 = dot(n_raw, n_raw);
      if (n_len2 > tolerance * tolerance) {
        const Scalar n_len = dev_sqrt(n_len2);
        for (std::size_t c = 0; c < 3; ++c) info.normal[c] = n_raw[c] / n_len;
        info.have_plane = true;
        info.plane_point = pts[k];
      }
    }
    else if (!info.have_volume) {
      const Scalar off_plane = dot(info.normal, ek);
      if (dev_abs(off_plane) > tolerance) {
        info.have_volume = true;
        info.volume_point = pts[k];
      }
    }
  }
  return info;
}

// The flat-shape projection: the origin clipped onto the triangle fan built
// from dv[0]'s chain, then the survivor edges; the host's
// project_onto_points carried over.
template <class Scalar>
__device__ inline cxd_dev::epa_projection<Scalar> dev_project_onto_points(
    const cxd_dev::vec<Scalar, 3>* dv, const std::size_t n, const double tolerance)
{
  constexpr std::size_t k_tri_cap = 5;  // one per interior index, n <= 7
  std::array<std::array<std::size_t, 3>, k_tri_cap> triangles{};
  std::size_t triangle_count = 0;
  for (std::size_t i = 1; i + 1 < n; ++i) {
    const cxd_dev::vec<Scalar, 3> e10 = sub(dv[i], dv[0]);
    const cxd_dev::vec<Scalar, 3> e20 = sub(dv[i + 1], dv[0]);
    const cxd_dev::vec<Scalar, 3> normal_v = cross(e20, e10);
    if (dot(normal_v, normal_v) > tolerance * tolerance)
      triangles[triangle_count++] = {std::size_t(0), i, i + 1};
  }

  const cxd_dev::vec<Scalar, 3> v0 = dv[triangles[0][0]];
  const cxd_dev::vec<Scalar, 3> v1 = dv[triangles[0][1]];
  const cxd_dev::vec<Scalar, 3> v2 = dv[triangles[0][2]];
  const cxd_dev::vec<Scalar, 3> e10 = sub(v1, v0);
  const cxd_dev::vec<Scalar, 3> e20 = sub(v2, v0);
  cxd_dev::vec<Scalar, 3> normal = cross(e20, e10);
  {
    const Scalar len = dev_sqrt(dot(normal, normal));
    for (std::size_t c = 0; c < 3; ++c) normal[c] = normal[c] / len;
  }

  const cxd_dev::vec<Scalar, 3> p = sub(cxd_dev::vec<Scalar, 3>{}, scale(normal, -dot(normal, v0)));
  Scalar best_dist2 = dev_inf<Scalar>();
  cxd_dev::vec<Scalar, 3> best_point = p;
  for (std::size_t tri = 0; tri < triangle_count; ++tri) {
    const cxd_dev::vec<Scalar, 3>& ta = dv[triangles[tri][0]];
    const cxd_dev::vec<Scalar, 3>& tb = dv[triangles[tri][1]];
    const cxd_dev::vec<Scalar, 3>& tc = dv[triangles[tri][2]];
    const cxd_dev::vec<Scalar, 3> tab = sub(tb, ta);
    const cxd_dev::vec<Scalar, 3> tac = sub(tc, ta);
    const cxd_dev::vec<Scalar, 3> tap = sub(p, ta);
    const Scalar d1 = dot(tab, tap);
    const Scalar d2 = dot(tac, tap);
    const Scalar d3 = dot(tab, tab);
    const Scalar d4 = dot(tac, tab);
    const Scalar d5 = dot(tab, tac);
    const Scalar d6 = dot(tac, tac);
    const Scalar denom = d3 * d6 - d4 * d5;
    if (dev_abs(denom) > tolerance * tolerance) {
      const Scalar u = (d1 * d6 - d2 * d4) / denom;
      const Scalar v = (d3 * d2 - d5 * d1) / denom;
      if (u >= Scalar(0) && v >= Scalar(0) && u + v <= Scalar(1)) {
        const cxd_dev::vec<Scalar, 3> cp = add(add(ta, scale(tab, u)), scale(tac, v));
        const Scalar dist2 = dot(cp, cp);
        if (dist2 < best_dist2) {
          best_dist2 = dist2;
          best_point = cp;
        }
      }
    }
    const std::size_t edges[3][2] = {
        {triangles[tri][0], triangles[tri][1]},
        {triangles[tri][1], triangles[tri][2]},
        {triangles[tri][2], triangles[tri][0]}};
    for (int escape = 0; escape < 3; ++escape) {
      const cxd_dev::vec<Scalar, 3> ea = dv[edges[escape][0]];
      const cxd_dev::vec<Scalar, 3> eb = dv[edges[escape][1]];
      const cxd_dev::vec<Scalar, 3> e = sub(eb, ea);
      const Scalar t2 =
          Scalar(dev_clamp01(-dot(e, sub(p, ea)) / dot(e, e)));
      const cxd_dev::vec<Scalar, 3> cp = add(ea, scale(e, t2));
      const Scalar dist2 = dot(cp, cp);
      if (dist2 < best_dist2) {
        best_dist2 = dist2;
        best_point = cp;
      }
    }
  }
  cxd_dev::epa_projection<Scalar> out{};
  out.depth = dev_sqrt(best_dist2);
  cxd_dev::vec<Scalar, 3> out_normal = best_point;
  const Scalar out_len = dev_sqrt(dot(out_normal, out_normal));
  if (out_len > Scalar(tolerance)) {
    for (std::size_t c = 0; c < 3; ++c) out_normal[c] = out_normal[c] / out_len;
  }
  else {
    out_normal = normal;
  }
  out.normal = out_normal;
  return out;
}

template <class Scalar>
__device__ inline cxd_dev::epa_projection<Scalar> dev_closest_on_segments(
    const cxd_dev::vec<Scalar, 3>* dv, const std::size_t n, const double tolerance)
{
  // Still collinear/coincident after bootstrapping: closest point over the
  // union of the segments the distinct set spans, then the point-line
  // answer the host falls back to.
  const cxd_dev::vec<Scalar, 3> zero{};
  Scalar best_dist2 = dev_inf<Scalar>();
  cxd_dev::vec<Scalar, 3> best_point = dv[0];
  if (n == 1) {
    best_point = dv[0];
    best_dist2 = dot(best_point, best_point);
  }
  else {
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = i + 1; j < n; ++j) {
        const cxd_dev::vec<Scalar, 3> ea = dv[i];
        const cxd_dev::vec<Scalar, 3> eb = dv[j];
        const cxd_dev::vec<Scalar, 3> e = sub(eb, ea);
        const Scalar t2 = Scalar(
            dev_clamp01(-dot(e, sub(ea, zero)) / dot(e, e)));
        const cxd_dev::vec<Scalar, 3> cp = add(ea, scale(e, t2));
        const Scalar dist2 = dot(cp, cp);
        if (dist2 < best_dist2) {
          best_dist2 = dist2;
          best_point = cp;
        }
      }
    }
  }
  cxd_dev::epa_projection<Scalar> out{};
  out.depth = dev_sqrt(best_dist2);
  cxd_dev::vec<Scalar, 3> out_normal = best_point;
  const Scalar out_len = dev_sqrt(dot(out_normal, out_normal));
  for (std::size_t c = 0; c < 3; ++c) {
    if (out_len > Scalar(tolerance))
      out_normal[c] = out_normal[c] / out_len;
    else
      out_normal[c] = (c == 0) ? Scalar(1) : Scalar(0);
  }
  out.normal = out_normal;
  return out;
}

template <class Scalar>
__device__ inline cxd_dev_answer<Scalar, 3> dev_epa_pair(
    const cxd_dev::pool_view<Scalar, 3>& pool,
    const cxd_dev_pair<Scalar, 3>& pair, const cxd_dev::vec<Scalar, 3> seeds[4],
    const cxd_dev_options& options)
{
  if (options.tolerance <= 0. || options.max_iterations == 0)
    return cxd_dev::refused<Scalar, 3>();
  if (pair.shape_id_A >= pool.entry_count ||
      pair.shape_id_B >= pool.entry_count ||
      pool.entries[pair.shape_id_A].vertex_count == 0 ||
      pool.entries[pair.shape_id_B].vertex_count == 0)
    return cxd_dev::refused<Scalar, 3>();

  cxd_dev::vec<Scalar, 3> a = seeds[0];
  cxd_dev::vec<Scalar, 3> b = seeds[1];
  cxd_dev::vec<Scalar, 3> c = seeds[2];
  cxd_dev::vec<Scalar, 3> d = seeds[3];

  const cxd_dev::vec<Scalar, 3> ab0 = sub(b, a);
  const cxd_dev::vec<Scalar, 3> ac0 = sub(c, a);
  const cxd_dev::vec<Scalar, 3> ad0 = sub(d, a);
  if (cxd_dev::dev_abs(dot(cross(ab0, ac0), ad0)) <=
      options.tolerance * options.tolerance * options.tolerance) {
    // Degenerate witness bootstrap, portal refinement and the flat/proj
    // fallbacks, in the host's order and bounds.
    constexpr std::size_t k_distinct_cap = 4 + 3;
    cxd_dev::vec<Scalar, 3> distinct[k_distinct_cap];
    std::size_t distinct_size = 0;
    const cxd_dev::vec<Scalar, 3> seed_pts[4] = {a, b, c, d};
    for (std::size_t vi = 0; vi < 4; ++vi) {
      const cxd_dev::vec<Scalar, 3>& v = seed_pts[vi];
      bool found = false;
      for (std::size_t u = 0; u < distinct_size; ++u) {
        const cxd_dev::vec<Scalar, 3> diff = sub(distinct[u], v);
        if (dot(diff, diff) <= options.tolerance * options.tolerance) {
          found = true;
          break;
        }
      }
      if (!found) distinct[distinct_size++] = v;
    }

    constexpr std::size_t max_bootstrap_attempts = 3;
    for (std::size_t attempt = 0; attempt < max_bootstrap_attempts;
         ++attempt) {
      const cxd_dev::degeneracy_info<Scalar> info =
          dev_analyze<Scalar>(distinct, distinct_size, options.tolerance);
      if (info.have_volume) break;

      constexpr std::size_t k_candidate_cap = 2;
      cxd_dev::vec<Scalar, 3> candidates[k_candidate_cap];
      std::size_t candidate_count = 0;
      if (!info.have_line) {
        candidates[candidate_count++] = cxd_dev::vec<Scalar, 3>{Scalar(1), Scalar(0), Scalar(0)};
      }
      else if (!info.have_plane) {
        // Escape collinearity: the host's oblique sum of every non-zero
        // axis-perpendicular, then the first perpendicular.
        const cxd_dev::vec<Scalar, 3> axes[3] = {cxd_dev::vec<Scalar, 3>{Scalar(1), Scalar(0), Scalar(0)},
                             cxd_dev::vec<Scalar, 3>{Scalar(0), Scalar(1), Scalar(0)},
                             cxd_dev::vec<Scalar, 3>{Scalar(0), Scalar(0), Scalar(1)}};
        cxd_dev::vec<Scalar, 3> perp_first{};
        cxd_dev::vec<Scalar, 3> perp_sum{};
        bool have_perp_first = false;
        for (std::size_t ai = 0; ai < 3; ++ai) {
          const cxd_dev::vec<Scalar, 3> cand = cross(info.line_dir, axes[ai]);
          const Scalar cand_len2 = dot(cand, cand);
          if (cand_len2 <= options.tolerance * options.tolerance) continue;
          const Scalar cand_len = dev_sqrt(cand_len2);
          cxd_dev::vec<Scalar, 3> unit_cand{};
          for (std::size_t ck = 0; ck < 3; ++ck)
            unit_cand[ck] = cand[ck] / cand_len;
          if (!have_perp_first) {
            perp_first = unit_cand;
            have_perp_first = true;
          }
          perp_sum = add(perp_sum, unit_cand);
        }
        const Scalar perp_sum_len2 = dot(perp_sum, perp_sum);
        if (perp_sum_len2 > options.tolerance * options.tolerance) {
          const Scalar perp_sum_len = dev_sqrt(perp_sum_len2);
          cxd_dev::vec<Scalar, 3> unit_sum{};
          for (std::size_t ck = 0; ck < 3; ++ck)
            unit_sum[ck] = perp_sum[ck] / perp_sum_len;
          candidates[candidate_count++] = unit_sum;
        }
        if (have_perp_first) candidates[candidate_count++] = perp_first;
      }
      else {
        // Escape coplanarity: the plane normal, then its opposite.
        candidates[candidate_count++] = info.normal;
        candidates[candidate_count++] = cxd_dev::neg(info.normal);
      }

      bool added = false;
      for (std::size_t ci = 0; ci < candidate_count; ++ci) {
        const cxd_dev::vec<Scalar, 3> dir = candidates[ci];
        const cxd_dev::vec<Scalar, 3> cv = cxd_dev::diff_support(pool, pair.shape_id_A,
                                            pair.shape_id_B, dir);
        bool is_new = true;
        for (std::size_t u = 0; u < distinct_size; ++u) {
          const cxd_dev::vec<Scalar, 3> diff = sub(distinct[u], cv);
          if (dot(diff, diff) <= options.tolerance * options.tolerance) {
            is_new = false;
            break;
          }
        }
        if (!is_new) continue;

        distinct[distinct_size++] = cv;
        added = true;
        break;
      }
      if (!added) break;
    }

    const cxd_dev::degeneracy_info<Scalar> final_info =
        dev_analyze<Scalar>(distinct, distinct_size, options.tolerance);

    if (final_info.have_volume) {
      // Portal refinement: replace the opposite vertex of the face the
      // origin is least enclosed by, until strict enclosure or budget.
      cxd_dev::vec<Scalar, 3> tet[4] = {distinct[0], final_info.line_point,
                    final_info.plane_point, final_info.volume_point};
      constexpr std::size_t face_defs[4][4] = {
          {0, 1, 2, 3}, {0, 1, 3, 2}, {0, 2, 3, 1}, {1, 2, 3, 0}};

      auto worst_face = [&]() {
        std::size_t worst = 0;
        Scalar worst_val = -dev_inf<Scalar>();
        cxd_dev::vec<Scalar, 3> worst_normal{};
        for (std::size_t fi = 0; fi < 4; ++fi) {
          const cxd_dev::vec<Scalar, 3> fa = tet[face_defs[fi][0]];
          const cxd_dev::vec<Scalar, 3> fb = tet[face_defs[fi][1]];
          const cxd_dev::vec<Scalar, 3> fc = tet[face_defs[fi][2]];
          const cxd_dev::vec<Scalar, 3> opp = tet[face_defs[fi][3]];
          cxd_dev::vec<Scalar, 3> n = cross(sub(fb, fa), sub(fc, fa));
          const Scalar n_len2 = dot(n, n);
          if (n_len2 <= options.tolerance * options.tolerance) continue;
          const Scalar n_len = dev_sqrt(n_len2);
          for (std::size_t ck = 0; ck < 3; ++ck) n[ck] = n[ck] / n_len;
          if (dot(n, sub(opp, fa)) > Scalar(0)) n = cxd_dev::neg(n);
          const Scalar val = dot(n, cxd_dev::neg(fa));
          if (val > worst_val) {
            worst_val = val;
            worst = fi;
            worst_normal = n;
          }
        }
        struct worst_t { std::size_t worst; Scalar val; cxd_dev::vec<Scalar, 3> normal; };
        return worst_t{worst, worst_val, worst_normal};
      };

      bool enclosing = false;
      constexpr std::size_t max_refine_attempts = 4;
      for (std::size_t refine = 0; refine < max_refine_attempts; ++refine) {
        const auto result = worst_face();
        if (result.val < -options.tolerance) {
          enclosing = true;
          break;
        }
        const cxd_dev::vec<Scalar, 3> cv = cxd_dev::diff_support(
            pool, pair.shape_id_A, pair.shape_id_B, result.normal);
        const cxd_dev::vec<Scalar, 3> diff = sub(cv, tet[face_defs[result.worst][3]]);
        if (dot(diff, diff) <= options.tolerance * options.tolerance) break;
        tet[face_defs[result.worst][3]] = cv;
      }

      if (enclosing) {
        a = tet[0];
        b = tet[1];
        c = tet[2];
        d = tet[3];
      }
      else {
        // The fan that existed before the volume-completing point came in.
        const auto proj = dev_project_onto_points<Scalar>(
            distinct, distinct_size - 1, options.tolerance);
        cxd_dev_answer<Scalar, 3> out{};
        out.depth = proj.depth;
        for (std::size_t ck = 0; ck < 3; ++ck) out.normal[ck] = proj.normal[ck];
        return out;
      }
    }
    else if (final_info.have_plane) {
      const auto proj = dev_project_onto_points<Scalar>(
          distinct, distinct_size, options.tolerance);
      cxd_dev_answer<Scalar, 3> out{};
      out.depth = proj.depth;
      for (std::size_t ck = 0; ck < 3; ++ck) out.normal[ck] = proj.normal[ck];
      return out;
    }
    else {
      const auto proj =
          dev_closest_on_segments<Scalar>(distinct, distinct_size,
                                          options.tolerance);
      cxd_dev_answer<Scalar, 3> out{};
      out.depth = proj.depth;
      for (std::size_t ck = 0; ck < 3; ++ck) out.normal[ck] = proj.normal[ck];
      return out;
    }
  }

  // Witness repair: a simplex that only touches the origin gets replaced by
  // the hull of four supports along the six axes that does enclose it.
  {
    auto encloses = [&](const cxd_dev::vec<Scalar, 3>* t) {
      const std::size_t face_list[4][4] = {
          {0, 1, 2, 3}, {0, 1, 3, 2}, {0, 2, 3, 1}, {1, 2, 3, 0}};
      for (std::size_t fi = 0; fi < 4; ++fi) {
        const cxd_dev::vec<Scalar, 3> p = t[face_list[fi][0]];
        const cxd_dev::vec<Scalar, 3> q = t[face_list[fi][1]];
        const cxd_dev::vec<Scalar, 3> r = t[face_list[fi][2]];
        const cxd_dev::vec<Scalar, 3> opp = t[face_list[fi][3]];
        cxd_dev::vec<Scalar, 3> n = cross(sub(q, p), sub(r, p));
        if (dot(n, n) <= options.tolerance * options.tolerance) return false;
        if (dot(n, sub(opp, p)) > Scalar(0)) n = cxd_dev::neg(n);
        if (dot(n, cxd_dev::neg(p)) >= -options.tolerance) return false;
      }
      return true;
    };

    cxd_dev::vec<Scalar, 3> witness[4] = {a, b, c, d};
    if (!encloses(witness)) {
      constexpr std::size_t k_points_cap = 10;  // 4 witnesses + 6 axes
      cxd_dev::vec<Scalar, 3> points[k_points_cap];
      std::size_t points_size = 0;
      for (std::size_t vi = 0; vi < 4; ++vi) points[points_size++] = witness[vi];
      const cxd_dev::vec<Scalar, 3> axes[6] = {cxd_dev::vec<Scalar, 3>{Scalar(1), Scalar(0), Scalar(0)},
                           cxd_dev::vec<Scalar, 3>{Scalar(-1), Scalar(0), Scalar(0)},
                           cxd_dev::vec<Scalar, 3>{Scalar(0), Scalar(1), Scalar(0)},
                           cxd_dev::vec<Scalar, 3>{Scalar(0), Scalar(-1), Scalar(0)},
                           cxd_dev::vec<Scalar, 3>{Scalar(0), Scalar(0), Scalar(1)},
                           cxd_dev::vec<Scalar, 3>{Scalar(0), Scalar(0), Scalar(-1)}};
      for (std::size_t ai = 0; ai < 6; ++ai) {
        const cxd_dev::vec<Scalar, 3> v = cxd_dev::diff_support(pool, pair.shape_id_A,
                                           pair.shape_id_B, axes[ai]);
        bool known = false;
        for (std::size_t u = 0; u < points_size; ++u) {
          const cxd_dev::vec<Scalar, 3> diff = sub(v, points[u]);
          if (dot(diff, diff) <= options.tolerance * options.tolerance) {
            known = true;
            break;
          }
        }
        if (!known) points[points_size++] = v;
      }
      for (std::size_t i = 0; i + 3 < points_size; ++i)
        for (std::size_t j = i + 1; j + 2 < points_size; ++j)
          for (std::size_t k = j + 1; k + 1 < points_size; ++k)
            for (std::size_t l = k + 1; l < points_size; ++l) {
              if (encloses(witness)) break;
              const cxd_dev::vec<Scalar, 3> trial[4] = {points[i], points[j], points[k],
                                    points[l]};
              if (encloses(trial)) {
                for (std::size_t wi = 0; wi < 4; ++wi) witness[wi] = trial[wi];
              }
            }
      a = witness[0];
      b = witness[1];
      c = witness[2];
      d = witness[3];
    }
  }

  // The face walk on the polytope arena.
  cxd_dev::vec<Scalar, 3> vertices[cxd_dev_vert_cap] = {};
  dev_face<Scalar> faces[cxd_dev_face_cap] = {};
  std::size_t vertex_count = 0;
  std::size_t face_count = 0;

  auto make_face = [&](const std::size_t i0, const std::size_t i1,
                       const std::size_t i2) {
    const cxd_dev::vec<Scalar, 3> v0 = vertices[i0];
    const cxd_dev::vec<Scalar, 3> v1 = vertices[i1];
    const cxd_dev::vec<Scalar, 3> v2 = vertices[i2];
    const cxd_dev::vec<Scalar, 3> e10 = sub(v1, v0);
    const cxd_dev::vec<Scalar, 3> e20 = sub(v2, v0);
    cxd_dev::vec<Scalar, 3> normal = cross(e20, e10);
    const Scalar len2 = dot(normal, normal);
    const Scalar len = dev_sqrt(len2);
    if (len > Scalar(options.tolerance)) {
      for (std::size_t ck = 0; ck < 3; ++ck) normal[ck] = normal[ck] / len;
    }
    dev_face<Scalar> f{};
    f.idx[0] = i0;
    f.idx[1] = i1;
    f.idx[2] = i2;
    f.normal = normal;
    f.distance = dot(normal, v0);
    return f;
  };

  vertices[vertex_count++] = a;
  vertices[vertex_count++] = b;
  vertices[vertex_count++] = c;
  vertices[vertex_count++] = d;

  faces[face_count++] = make_face(0, 2, 1);  // ACB
  faces[face_count++] = make_face(0, 1, 3);  // ABD
  faces[face_count++] = make_face(0, 3, 2);  // ADC
  faces[face_count++] = make_face(1, 2, 3);  // BCD

  for (std::size_t i = 0; i < face_count; ++i) {
    if (faces[i].distance < Scalar(0)) {
      faces[i].distance = -faces[i].distance;
      faces[i].normal = cxd_dev::neg(faces[i].normal);
      const std::size_t swap = faces[i].idx[0];
      faces[i].idx[0] = faces[i].idx[1];
      faces[i].idx[1] = swap;
    }
  }

  const std::size_t max_polytope =
      pool.entries[pair.shape_id_A].vertex_count +
      pool.entries[pair.shape_id_B].vertex_count;

  const auto answer_from = [&](const Scalar depth,
                               const cxd_dev::vec<Scalar, 3>& normal) {
    cxd_dev_answer<Scalar, 3> out{};
    out.depth = depth;
    for (std::size_t ck = 0; ck < 3; ++ck) out.normal[ck] = normal[ck];
    return out;
  };

  for (std::size_t iteration = 0; iteration < options.max_iterations;
       ++iteration) {
    Scalar closest_distance = dev_inf<Scalar>();
    std::size_t closest_face = 0;

    for (std::size_t i = 0; i < face_count; ++i) {
      if (faces[i].distance < closest_distance) {
        closest_distance = faces[i].distance;
        closest_face = i;
      }
    }

    const cxd_dev::vec<Scalar, 3>& closest_normal = faces[closest_face].normal;
    const cxd_dev::vec<Scalar, 3> candidate = cxd_dev::diff_support(
        pool, pair.shape_id_A, pair.shape_id_B, closest_normal);
    const Scalar candidate_distance =
        dot(candidate, closest_normal);
    if (candidate_distance - closest_distance <= options.tolerance)
      return answer_from(closest_distance, closest_normal);

    bool duplicate = false;
    for (std::size_t i = 0; i < vertex_count; ++i) {
      const cxd_dev::vec<Scalar, 3> diff = sub(vertices[i], candidate);
      if (dot(diff, diff) <= options.tolerance * options.tolerance) {
        duplicate = true;
        break;
      }
    }
    if (duplicate)
      return answer_from(closest_distance, closest_normal);

    if (vertex_count >= max_polytope || vertex_count >= cxd_dev_vert_cap)
      return cxd_dev::refused<Scalar, 3>();

    const std::size_t new_index = vertex_count;
    vertices[vertex_count++] = candidate;

    for (std::size_t i = 0; i < face_count; ++i) {
      const cxd_dev::vec<Scalar, 3> v0 = vertices[faces[i].idx[0]];
      const Scalar dot_v = dot(faces[i].normal, sub(candidate, v0));
      faces[i].removed = dot_v > options.tolerance;
    }

    bool all_faces_removed = true;
    for (std::size_t i = 0; i < face_count; ++i) {
      if (!faces[i].removed) {
        all_faces_removed = false;
        break;
      }
    }
    if (all_faces_removed)
      return answer_from(closest_distance, closest_normal);

    // Boundary edges carried by exactly one retained face; the new faces
    // are triangulated onto them, addressed by index (the region's storage
    // never moves in this fixed arena, matching the host's index rule).
    const std::size_t retained = face_count;
    for (std::size_t fi = 0; fi < retained; ++fi) {
      if (faces[fi].removed) continue;
      for (std::size_t e = 0; e < 3; ++e) {
        const std::size_t e0 = faces[fi].idx[e];
        const std::size_t e1 = faces[fi].idx[(e + 1) % 3];
        const std::size_t boundary0 = e0 < e1 ? e0 : e1;
        const std::size_t boundary1 = e0 < e1 ? e1 : e0;

        std::size_t carriers = 0;
        for (std::size_t fj = 0; fj < retained; ++fj) {
          if (faces[fj].removed) continue;
          for (std::size_t f = 0; f < 3; ++f) {
            const std::size_t f0 = faces[fj].idx[f];
            const std::size_t f1 = faces[fj].idx[(f + 1) % 3];
            if ((f0 < f1 ? f0 : f1) == boundary0 &&
                (f0 < f1 ? f1 : f0) == boundary1)
              ++carriers;
          }
        }
        if (carriers != 1) continue;

        if (face_count + 1 >= cxd_dev_face_cap)
          return cxd_dev::refused<Scalar, 3>();
        dev_face<Scalar> fresh = make_face(new_index, e0, e1);
        if (fresh.distance < Scalar(0)) {
          fresh.distance = -fresh.distance;
          fresh.normal = cxd_dev::neg(fresh.normal);
          const std::size_t swap = fresh.idx[0];
          fresh.idx[0] = fresh.idx[1];
          fresh.idx[1] = swap;
        }
        faces[face_count++] = fresh;
      }
    }

    // Drop the removed faces, keeping the survivors and then the new ones
    // in order.
    {
      std::size_t kept = 0;
      for (std::size_t i = 0; i < face_count; ++i) {
        if (!faces[i].removed) {
          if (kept != i) faces[kept] = faces[i];
          ++kept;
        }
      }
      face_count = kept;
    }
  }

  return cxd_dev::refused<Scalar, 3>();
}

}  // namespace cxd_dev_3d

// ---------------------------------------------------------------------------
// Kernels and the host launchers: one launch per pass, host-side compaction
// between them (the makeup the validated spike used).
// ---------------------------------------------------------------------------

template <class Scalar, std::size_t Dim>
__global__ void cxd_dev_gjk_kernel(
    cxd_dev::pool_view<Scalar, Dim> pool, const cxd_dev_pair<Scalar, Dim>* pairs,
    std::size_t pair_count, cxd_dev_options options,
    cxd_dev_gjk_row<Scalar, Dim>* rows)
{
  const std::size_t i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
  if (i >= pair_count) return;
  cxd_dev_gjk_row<Scalar, Dim> row{};
  // The batch-level option rejection: every pair invalid, exactly as
  // gjk_batch rejects the whole list.
  if (options.tolerance <= 0. || options.max_iterations == 0) {
    row.st = cxd_dev_invalid_input;
  }
  else if constexpr (Dim == 2) {
    row.st = cxd_dev_2d::dev_gjk_pair<Scalar>(pool, pairs[i], options, row);
  }
  else {
    row.st = cxd_dev_3d::dev_gjk_pair<Scalar>(pool, pairs[i], options, row);
  }
  rows[i] = row;
}

template <class Scalar, std::size_t Dim>
__global__ void cxd_dev_epa_kernel(
    cxd_dev::pool_view<Scalar, Dim> pool, const cxd_dev_pair<Scalar, Dim>* pairs,
    std::size_t pair_count, const std::uint32_t* hit_ids,
    std::size_t hit_count, const cxd_dev_gjk_row<Scalar, Dim>* rows,
    cxd_dev_options options, cxd_dev_answer<Scalar, Dim>* answers)
{
  const std::size_t i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
  if (i >= hit_count) return;
  // Unanswered in place: a hit naming a pair this batch never evaluated
  // keeps the NaN row, the way epa_batch pre-fills its answers.
  cxd_dev_answer<Scalar, Dim> out = cxd_dev::refused<Scalar, Dim>();
  const std::uint32_t hit = hit_ids[i];
  if (hit < pair_count) {
    cxd_dev::vec<Scalar, Dim> seed[Dim + 1];
    for (std::size_t s = 0; s <= Dim; ++s)
      for (std::size_t c = 0; c < Dim; ++c)
        seed[s][c] = rows[hit].seed[s * Dim + c];
    if constexpr (Dim == 2) {
      out = cxd_dev_2d::dev_epa_pair<Scalar>(pool, pairs[hit], seed, options);
    }
    else {
      out = cxd_dev_3d::dev_epa_pair<Scalar>(pool, pairs[hit], seed, options);
    }
  }
  answers[i] = out;
}

// The host side of the backend: pool upload owned once per corpus, one
// launch per pass. Same seam vocabulary as the phase-3 contract: buffers in
// (flat vertices, entries, pairs), rows out (gjk statuses with seeds
// device-side, answers).
namespace cxd_device {

inline void check(const cudaError_t err, const char* what)
{
  if (err != cudaSuccess) {
    throw std::runtime_error(std::string(what) + ": " +
                             cudaGetErrorString(err));
  }
}

template <class Scalar, std::size_t Dim>
class device_pool {
 public:
  void upload(const Scalar* host_vertices, const std::size_t n_vertices,
              const cxd_dev_entry* host_entries, const std::size_t n_entries)
  {
    free_all();
    vertex_count_ = n_vertices;
    entry_count_ = n_entries;
    if (n_vertices != 0) {
      check(cudaMalloc(&dev_vertices_, sizeof(Scalar) * n_vertices),
            "cudaMalloc pool vertices");
      check(cudaMemcpy(dev_vertices_, host_vertices, sizeof(Scalar) * n_vertices,
                       cudaMemcpyHostToDevice),
            "upload pool vertices");
    }
    if (n_entries != 0) {
      check(cudaMalloc(&dev_entries_, sizeof(cxd_dev_entry) * n_entries),
            "cudaMalloc pool entries");
      check(cudaMemcpy(dev_entries_, host_entries,
                       sizeof(cxd_dev_entry) * n_entries,
                       cudaMemcpyHostToDevice),
            "upload pool entries");
    }
  }

  device_pool() = default;
  device_pool(const device_pool&) = delete;
  device_pool& operator=(const device_pool&) = delete;
  ~device_pool() { free_all(); }

  cxd_dev::pool_view<Scalar, Dim> view() const
  {
    return cxd_dev::pool_view<Scalar, Dim>{dev_vertices_, dev_entries_,
                                           static_cast<std::uint32_t>(
                                               entry_count_)};
  }

 private:
  void free_all()
  {
    if (dev_vertices_ != nullptr) {
      cudaFree(dev_vertices_);
      dev_vertices_ = nullptr;
    }
    if (dev_entries_ != nullptr) {
      cudaFree(dev_entries_);
      dev_entries_ = nullptr;
    }
    vertex_count_ = 0;
    entry_count_ = 0;
  }

  Scalar* dev_vertices_ = nullptr;
  cxd_dev_entry* dev_entries_ = nullptr;
  std::size_t vertex_count_ = 0;
  std::size_t entry_count_ = 0;
};

// The first pass: one status row per pair, seeds retained device-side.
template <class Scalar, std::size_t Dim>
std::vector<cxd_dev_gjk_row<Scalar, Dim>> gjk_batch(
    const device_pool<Scalar, Dim>& pool,
    const std::vector<cxd_dev_pair<Scalar, Dim>>& pairs,
    const cxd_dev_options& options)
{
  std::vector<cxd_dev_gjk_row<Scalar, Dim>> results(pairs.size());
  const std::size_t n = pairs.size();
  if (n == 0) return results;

  cxd_dev_pair<Scalar, Dim>* dev_pairs = nullptr;
  cxd_dev_gjk_row<Scalar, Dim>* dev_rows = nullptr;
  check(cudaMalloc(&dev_pairs, sizeof(cxd_dev_pair<Scalar, Dim>) * n),
        "cudaMalloc pairs");
  check(cudaMemcpy(dev_pairs, pairs.data(),
                   sizeof(cxd_dev_pair<Scalar, Dim>) * n,
                   cudaMemcpyHostToDevice),
        "upload pairs");
  check(cudaMalloc(&dev_rows, sizeof(cxd_dev_gjk_row<Scalar, Dim>) * n),
        "cudaMalloc rows");

  const unsigned block = 256;
  const unsigned grid =
      static_cast<unsigned>((n + block - 1) / block);
  cxd_dev_gjk_kernel<Scalar, Dim><<<grid, block>>>(
      pool.view(), dev_pairs, n, options, dev_rows);
  check(cudaGetLastError(), "gjk kernel launch");
  check(cudaDeviceSynchronize(), "gjk kernel sync");
  check(cudaMemcpy(results.data(), dev_rows,
                   sizeof(cxd_dev_gjk_row<Scalar, Dim>) * n,
                   cudaMemcpyDeviceToHost),
        "download rows");

  cudaFree(dev_pairs);
  cudaFree(dev_rows);
  return results;
}

// The second pass: one answer row per hit, over the GJK rows the first pass
// wrote (kept device-side through this signature).
template <class Scalar, std::size_t Dim>
std::vector<cxd_dev_answer<Scalar, Dim>> epa_batch(
    const device_pool<Scalar, Dim>& pool,
    const std::vector<cxd_dev_pair<Scalar, Dim>>& pairs,
    const std::vector<std::uint32_t>& hit_pairs,
    const std::vector<cxd_dev_gjk_row<Scalar, Dim>>& results,
    const cxd_dev_options& options)
{
  std::vector<cxd_dev_answer<Scalar, Dim>> answers(hit_pairs.size());
  const std::size_t nh = hit_pairs.size();
  if (nh == 0) return answers;

  cxd_dev_pair<Scalar, Dim>* dev_pairs = nullptr;
  std::uint32_t* dev_hits = nullptr;
  cxd_dev_gjk_row<Scalar, Dim>* dev_rows = nullptr;
  cxd_dev_answer<Scalar, Dim>* dev_answers = nullptr;
  check(cudaMalloc(&dev_pairs, sizeof(cxd_dev_pair<Scalar, Dim>) * pairs.size()),
        "cudaMalloc pairs");
  check(cudaMemcpy(dev_pairs, pairs.data(),
                   sizeof(cxd_dev_pair<Scalar, Dim>) * pairs.size(),
                   cudaMemcpyHostToDevice),
        "upload pairs");
  check(cudaMalloc(&dev_hits, sizeof(std::uint32_t) * nh), "cudaMalloc hits");
  check(cudaMemcpy(dev_hits, hit_pairs.data(), sizeof(std::uint32_t) * nh,
                   cudaMemcpyHostToDevice),
        "upload hits");
  check(cudaMalloc(&dev_rows,
                   sizeof(cxd_dev_gjk_row<Scalar, Dim>) * results.size()),
        "cudaMalloc gjk rows");
  check(cudaMemcpy(dev_rows, results.data(),
                   sizeof(cxd_dev_gjk_row<Scalar, Dim>) * results.size(),
                   cudaMemcpyHostToDevice),
        "upload gjk rows");
  check(cudaMalloc(&dev_answers, sizeof(cxd_dev_answer<Scalar, Dim>) * nh),
        "cudaMalloc answers");

  const unsigned block = 256;
  const unsigned grid = static_cast<unsigned>((nh + block - 1) / block);
  cxd_dev_epa_kernel<Scalar, Dim><<<grid, block>>>(
      pool.view(), dev_pairs, pairs.size(), dev_hits, nh, dev_rows, options,
      dev_answers);
  check(cudaGetLastError(), "epa kernel launch");
  check(cudaDeviceSynchronize(), "epa kernel sync");
  check(cudaMemcpy(answers.data(), dev_answers,
                   sizeof(cxd_dev_answer<Scalar, Dim>) * nh,
                   cudaMemcpyDeviceToHost),
        "download answers");

  cudaFree(dev_pairs);
  cudaFree(dev_hits);
  cudaFree(dev_rows);
  cudaFree(dev_answers);
  return answers;
}

}  // namespace cxd_device

#endif  // CXD_BATCH_DEVICE_HPP
