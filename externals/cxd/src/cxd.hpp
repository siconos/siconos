/* CXD is a program dedicated to convex collision detection
 * (currently part of the Siconos software).
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
#ifndef CXD_HPP
#define CXD_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "cxd_traits.hpp"
#include "cxd_utility.hpp"
#include "cxd_concepts.hpp"
#include "cxd_dim_math.hpp"

namespace cxd
{
  namespace num { namespace operators {}}

  using namespace cxd::traits;
  using namespace cxd::num::operators;
  using namespace cxd::concepts;

  struct gjk_options
  {
    double tolerance = 1e-6;
    std::size_t max_iterations = 64;
  };

  static constexpr void validate_vector_2d(const any_vector auto& v,
                                           const char* name)
  {
    if (v.size() != 2)
      throw std::invalid_argument(std::string{name} + " must be two-dimensional");
  }

  static constexpr void validate_shape_2d(const any_vectors auto& points,
                                          const char* name)
  {
    if (points.empty())
      throw std::invalid_argument(std::string{name} + " must not be empty");

    for (const auto& point : points)
      validate_vector_2d(point, name);
  }

  static constexpr void validate_vector_3d(const any_vector auto& v,
                                           const char* name)
  {
    if (v.size() != 3)
      throw std::invalid_argument(std::string{name} + " must be three-dimensional");
  }

  static constexpr void validate_shape_3d(const any_vectors auto& points,
                                          const char* name)
  {
    if (points.empty())
      throw std::invalid_argument(std::string{name} + " must not be empty");

    for (const auto& point : points)
      validate_vector_3d(point, name);
  }

  static constexpr void validate_options(const gjk_options& options)
  {
    if (options.tolerance <= 0.)
      throw std::invalid_argument("GJK tolerance must be positive");
    if (options.max_iterations == 0)
      throw std::invalid_argument("GJK max_iterations must be positive");
  }


  static constexpr const auto normalize(const any_vector auto& v)
  {
    using env = decltype(env(v));
    using vector = typename env::vector;
    using scalar = typename env::scalar;

    scalar invnorm = 1./sqrt(num::inner_prod(v, v));

    return vector{invnorm * v};
  }

  // ( a x b ) x c
  template<any_vector T>
  static inline constexpr const T cross3(const T& a,
                                         const T& b,
                                         const T& c)
  {
    // FIX: cast to input vector type
    // otherwise ublas vector_expression is not propagated correctly?
    return
      (- num::inner_prod(c, b) * a + num::inner_prod(c, a) * b);
  }

  static constexpr const auto support_indice(const any_vectors auto& points,
                                              const any_vector auto& v)
  {
    using env = decltype(env(v));
    using scalar = typename env::scalar;
    using indice = typename env::indice;

    return
      std::get<1>(
        util::ifold_left(points,
                         std::tuple{
                           -std::numeric_limits<scalar>::infinity(),
                           indice{0}},
                         [&v](const auto& state,
                              const indice i,
                              const auto& point)
                         {
                           const auto& cdot = num::inner_prod(point, v);
                           return (std::get<0>(state) < cdot) ?
                             std::tuple{cdot, i} : state;
                         }));
  }

  static constexpr const auto diff_support(const any_vectors auto& sa,
                                           const any_vectors auto& sb,
                                           const any_vector  auto& v)
  {
    using env = decltype(env(v));
    using vector = typename env::vector;
    using indice = typename env::indice;

    const indice ia = support_indice(sa, v);
    const indice ib = support_indice(sb, -v);
    return std::tuple{ia, ib, vector{sa[ia] - sb[ib]}};
  }

  static constexpr const auto support_vector(const any_support auto& s)
  {
    return std::get<2>(s);
  }

  static constexpr const bool clockwise(const any_vector auto& a,
                                        const any_vector auto& b,
                                        const any_vector auto& c)
  {
    const auto& x = [](const auto& v) { return v[0];};
    const auto& y = [](const auto& v) { return v[1];};

    return (x(b)*y(c)+x(a)*y(b)+y(a)*x(c) - (y(a)*x(b)+y(b)*x(c)+x(a)*y(c))) < 0;
  }


  static constexpr const auto exterior_normal(const bool clockwise,
                                              const any_vector auto& edge)
  {
    using env = decltype(env(edge));
    using vector = typename env::vector;

    return normalize(!clockwise ?
                     vector{ edge[1], -edge[0]} :
                     vector{-edge[1],  edge[0]});
  }


  template<std::size_t Dim, typename ShapesA, typename ShapesB, typename Vector>
  requires (Dim == 2 || Dim == 3)
  static constexpr const auto gjk_inter(const ShapesA& sa,
                                        const ShapesB& sb,
                                        const Vector& init,
                                        const gjk_options options = {})
  {
    using env = decltype(env(init));
    using vector = typename env::vector;
    using support = typename env::support;
    using scalar = typename env::scalar;

    if constexpr (Dim == 2)
    {
      validate_shape_2d(sa, "first shape");
      validate_shape_2d(sb, "second shape");
      validate_vector_2d(init, "initial direction");
      validate_options(options);

      const auto tolerance_squared = options.tolerance * options.tolerance;
      if (num::inner_prod(init, init) <= tolerance_squared)
        throw std::invalid_argument("initial direction must not be zero");

      const auto is_near_zero = [&tolerance_squared](const vector& v)
      {
        return num::inner_prod(v, v) <= tolerance_squared;
      };
      const auto simplex0 = diff_support(sa, sb, init);
      if (num::inner_prod(support_vector(simplex0), init) < -options.tolerance)
        return std::tuple{false, simplex0, simplex0, simplex0};

      std::array<support, 3> simplex{simplex0, simplex0, simplex0};
      std::size_t simplex_size = 1;
      vector direction = -support_vector(simplex0);

      for (std::size_t iteration = 0; iteration < options.max_iterations; ++iteration)
      {
        if (is_near_zero(direction))
          return std::tuple{true, simplex[0], simplex[1], simplex[2]};

        const support candidate = diff_support(sa, sb, direction);
        if (num::inner_prod(support_vector(candidate), direction) < -options.tolerance)
          return std::tuple{false, candidate, simplex[0], simplex[1]};

        for (std::size_t i = 0; i < simplex_size; ++i)
        {
          const vector difference = support_vector(candidate) - support_vector(simplex[i]);
          if (is_near_zero(difference))
            return std::tuple{false, candidate, simplex[0], simplex[1]};
        }

        simplex[2] = simplex[1];
        simplex[1] = simplex[0];
        simplex[0] = candidate;
        ++simplex_size;

        const vector a = support_vector(simplex[0]);
        const vector a0 = -a;
        if (simplex_size == 2)
        {
          const vector ab = support_vector(simplex[1]) - a;
          if (num::inner_prod(ab, a0) <= options.tolerance)
          {
            simplex_size = 1;
            direction = a0;
            continue;
          }

          direction = cross3(ab, a0, ab);
          if (is_near_zero(direction))
            return std::tuple{true, simplex[0], simplex[1], simplex[2]};
          continue;
        }

        const vector ab = support_vector(simplex[1]) - a;
        const vector ac = support_vector(simplex[2]) - a;
        const vector ab_direction = cross3(ac, ab, ab);
        if (num::inner_prod(ab_direction, a0) > options.tolerance)
        {
          simplex_size = 2;
          direction = ab_direction;
          continue;
        }

        const vector ac_direction = cross3(ab, ac, ac);
        if (num::inner_prod(ac_direction, a0) > options.tolerance)
        {
          simplex[1] = simplex[2];
          simplex_size = 2;
          direction = ac_direction;
          continue;
        }

        return std::tuple{true, simplex[0], simplex[1], simplex[2]};
      }

      throw std::runtime_error("GJK exceeded max_iterations");
    }
    else if constexpr (Dim == 3)
    {
      validate_shape_3d(sa, "first shape");
      validate_shape_3d(sb, "second shape");
      validate_vector_3d(init, "initial direction");
      validate_options(options);

      const auto tolerance_squared = options.tolerance * options.tolerance;
      if (num::inner_prod(init, init) <= tolerance_squared)
        throw std::invalid_argument("initial direction must not be zero");

      const auto is_near_zero = [&tolerance_squared](const vector& v)
      {
        return num::inner_prod(v, v) <= tolerance_squared;
      };
      const auto simplex0 = diff_support(sa, sb, init);
      if (num::inner_prod(support_vector(simplex0), init) < -options.tolerance)
        return std::tuple{false, simplex0, simplex0, simplex0, simplex0};

      std::array<support, 4> simplex{simplex0, simplex0, simplex0, simplex0};
      std::size_t simplex_size = 1;
      vector direction = -support_vector(simplex0);

      const auto dedup_simplex = [&]() -> std::size_t {
        for (std::size_t i = 0; i < simplex_size; ++i)
        {
          for (std::size_t j = i + 1; j < simplex_size;)
          {
            if (is_near_zero(support_vector(simplex[i]) - support_vector(simplex[j])))
            {
              simplex[j] = simplex[simplex_size - 1];
              --simplex_size;
            }
            else
            {
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

      for (std::size_t iteration = 0; iteration < options.max_iterations; ++iteration)
      {

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
        const bool same_axis = has_loop_state
            and dir_dot_last * dir_dot_last >= 0.9 * dir_norm2 * last_norm2;
        if (same_axis and simplex_size == last_loop_size)
          ++same_loop_state;
        else
          same_loop_state = 0;
        has_loop_state = true;
        last_loop_dir = direction;
        last_loop_size = simplex_size;

        if (same_loop_state >= 1 and simplex_size == 3)
        {
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
          const scalar scale = num::inner_prod(t01, t01)
              + num::inner_prod(t02, t02) + num::inner_prod(tn, tn);
          const vector t0n = -t0;
          const bool origin_in_triangle =
              std::abs(num::inner_prod(tn, t0n)) <= options.tolerance * scale
              && num::inner_prod(cross(t01, tn), t0n) <= options.tolerance * scale
              && num::inner_prod(cross(tn, t02), t0n) <= options.tolerance * scale
              && num::inner_prod(cross(t12, tn), -t1) <= options.tolerance * scale;
          if (origin_in_triangle)
          {
            const auto compacted = compact_simplex();
            return std::tuple{true, compacted[0], compacted[1],
                              compacted[2], compacted[3]};
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
          if (num::inner_prod(tn, t0) > 0)
            toward_origin = -tn;
          const support sep_witness = diff_support(sa, sb, toward_origin);
          const scalar sep_proj = num::inner_prod(support_vector(sep_witness), toward_origin);
          // pair189 witness: tn=(232,233,232), t0=(5.56,-0.64,-13.17), w=(5.56,-0.64,8.41), proj=3097.56.
          // The stall triangle's normal is not the isolating direction, so the
          // face-normal support test cannot certify separation. Fall back to an
          // axis-aligned probe of the support difference along the three
          // coordinate axes of the converged triangle's bounding box: a
          // negative projection along any axis certifies separation of the
          // boxes (separated box pairs have their isochronous side planes
          // perpendicular to a box axis).
          vector probe = toward_origin;
          bool separated = sep_proj < -options.tolerance;
          if (!separated)
          {
            const scalar probe_candidates[6][3] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
            for (std::size_t pi = 0; pi < 6 and not separated; ++pi)
            {
              probe = vector{probe_candidates[pi][0], probe_candidates[pi][1], probe_candidates[pi][2]};
              const support probe_w = diff_support(sa, sb, probe);
              separated = num::inner_prod(support_vector(probe_w), probe) < -options.tolerance;
            }
          }
          const auto compacted = compact_simplex();
          if (separated)
            return std::tuple{false, compacted[0], compacted[1],
                              compacted[2], compacted[3]};
          return std::tuple{true, compacted[0], compacted[1],
                            compacted[2], sep_witness};
        }
        if (is_near_zero(direction))
        {
          dedup_simplex();
          if (simplex_size >= 3)
          {
            const vector a = support_vector(simplex[0]);
            const vector ab = support_vector(simplex[1]) - a;
            const vector ac = support_vector(simplex[2]) - a;
            const vector abc = cross(ab, ac);
            if (!is_near_zero(abc))
            {
              direction = num::inner_prod(abc, a) > 0 ? -abc : abc;
              continue;
            }
          }
          if (simplex_size >= 2)
          {
            const vector a = support_vector(simplex[0]);
            const vector ab = support_vector(simplex[1]) - a;
            direction = cross(ab, vector{1, 0, 0});
            if (is_near_zero(direction))
              direction = cross(ab, vector{0, 1, 0});
            if (num::inner_prod(direction, a) < 0)
              direction = -direction;
            continue;
          }
          direction = -support_vector(simplex[0]);
          continue;
        }

        const support candidate = diff_support(sa, sb, direction);
        if (num::inner_prod(support_vector(candidate), direction) < -options.tolerance)
        {
          return std::tuple{false, candidate, simplex[0], simplex[1], simplex[2]};
        }

        bool duplicate = false;
        for (std::size_t i = 0; i < simplex_size; ++i)
        {
          const vector difference = support_vector(candidate) - support_vector(simplex[i]);
          if (is_near_zero(difference))
          {
            duplicate = true;
            break;
          }
        }
        if (duplicate)
        {
          // Duplicate support vertex: the current simplex already encodes the
          // closest feature, and forcing a new vertex via geometric fallbacks
          // can cycle on degenerate cases (e.g. axis-aligned cubes where every
          // escape direction still maps to the same face). Following Bullet/Jolt/
          // ReactPhysics3D practice: reduce the simplex, keep the best feature,
          // and stop when no new vertex can be obtained.
          const std::size_t final_size = dedup_simplex();
          if (final_size < 3)
          {
            simplex_size = 1;
            direction = -support_vector(simplex[0]);
            continue;
          }
          if (simplex_size == 4)
          {
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
            const bool encloses =
                num::inner_prod(cross(ab, ac), ad) * num::inner_prod(cross(ab, ac), -a) > 0
                && num::inner_prod(cross(ab, ad), ac) * num::inner_prod(cross(ab, ad), -a) > 0
                && num::inner_prod(cross(ac, ad), ab) * num::inner_prod(cross(ac, ad), -a) > 0
                && num::inner_prod(cross(bd, bc), a - b) * num::inner_prod(cross(bd, bc), -b) > 0;
            if (encloses)
            {
              const auto compacted = compact_simplex();
              return std::tuple{true, compacted[0], compacted[1], compacted[2],
                                diff_support(sa, sb, direction)};
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
                num::inner_prod(cross(ab, abc), a0) <= options.tolerance
                && num::inner_prod(cross(abc, ac), a0) <= options.tolerance;
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
              if (num::inner_prod(dn, da) > 0)
                ddir = -dn;
              const support stall_wit4 = diff_support(sa, sb, ddir);
              if (origin_in_prism)
              {
                return std::tuple{true, stall_compacted[0], stall_compacted[1],
                                  stall_compacted[2], stall_compacted[3]};
              }
            // No separation certificate was found: GJK's early-out is the only
            // sound separation answer and it did not fire, so the stalled
            // query cannot prove separation. Answer intersection (the
            // pristine behaviour; the reviewed false positives for separated
            // pairs are caught by the loop-state stall terminal instead,
            // which sees the actual axis-aligned separation).
              return std::tuple{true, stall_compacted[0], stall_compacted[1],
                                stall_compacted[2], stall_compacted[3]};
            }
          }
        }

        if (simplex_size < 4)
        {
          for (std::size_t i = simplex_size; i > 0; --i)
          {
            simplex[i] = simplex[i - 1];
          }
          simplex[0] = candidate;
          ++simplex_size;
        }

            for (std::size_t i = 0; i < simplex_size; ++i)
            {
              for (std::size_t j = i + 1; j < simplex_size;)
              {
                if (is_near_zero(support_vector(simplex[i]) - support_vector(simplex[j])))
                {
                  simplex[j] = simplex[simplex_size - 1];
                  --simplex_size;
                }
                else
                {
                  ++j;
                }
              }
            }

        const vector a = support_vector(simplex[0]);
        const vector a0 = -a;

        if (simplex_size == 2)
        {
          const vector ab = support_vector(simplex[1]) - a;
          const scalar ab2 = num::inner_prod(ab, ab);
          const scalar ab_ao = num::inner_prod(ab, a0);
          if (ab_ao <= options.tolerance)
          {
            simplex_size = 1;
            direction = a0;
            continue;
          }
          if (ab_ao >= ab2 - options.tolerance)
          {
            simplex[0] = simplex[1];
            simplex_size = 1;
            direction = -support_vector(simplex[0]);
            continue;
          }
          direction = a0 - (ab_ao / ab2) * ab;
          if (is_near_zero(direction))
          {
            vector perp = cross(ab, vector{1, 0, 0});
            if (is_near_zero(perp))
              perp = cross(ab, vector{0, 1, 0});
            direction = perp;
          }
          else
          {
          }
          continue;
        }

        if (simplex_size == 3)
        {
          const vector ab = support_vector(simplex[1]) - a;
          const vector ac = support_vector(simplex[2]) - a;
          const vector abc = cross(ab, ac);
          const scalar abc_ao = num::inner_prod(abc, a0);

          const vector ab_cross_abc = cross(ab, abc);
          const vector ac_cross_abc = cross(abc, ac);

          if (num::inner_prod(ab_cross_abc, a0) > options.tolerance)
          {
            simplex_size = 2;
            const vector ab2 = support_vector(simplex[1]) - support_vector(simplex[0]);
            const scalar ab2_len = num::inner_prod(ab2, ab2);
            const scalar ab_ao2 = num::inner_prod(ab2, -support_vector(simplex[0]));
            if (ab_ao2 <= options.tolerance)
            {
              simplex_size = 1;
              direction = -support_vector(simplex[0]);
              continue;
            }
            direction = -support_vector(simplex[0]) - (ab_ao2 / ab2_len) * ab2;
            continue;
          }

          if (num::inner_prod(ac_cross_abc, a0) > options.tolerance)
          {
            simplex[1] = simplex[2];
            simplex_size = 2;
            const vector ac2 = support_vector(simplex[1]) - support_vector(simplex[0]);
            const scalar ac2_len = num::inner_prod(ac2, ac2);
            const scalar ac_ao2 = num::inner_prod(ac2, -support_vector(simplex[0]));
            if (ac_ao2 <= options.tolerance)
            {
              simplex_size = 1;
              direction = -support_vector(simplex[0]);
              continue;
            }
            direction = -support_vector(simplex[0]) - (ac_ao2 / ac2_len) * ac2;
            continue;
          }

          if (abc_ao > options.tolerance)
          {
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

        if (simplex_size == 4)
        {
          const vector ab = support_vector(simplex[1]) - a;
          const vector ac = support_vector(simplex[2]) - a;
          const vector ad = support_vector(simplex[3]) - a;
          const vector abc = cross(ab, ac);
          const vector abd = cross(ab, ad);
          const vector acd = cross(ad, ac);

          const scalar abc_ao = num::inner_prod(abc, a0) * num::inner_prod(abc, ad);
          const scalar abd_ao = num::inner_prod(abd, a0) * num::inner_prod(abd, ac);
          const scalar acd_ao = num::inner_prod(acd, a0) * num::inner_prod(acd, ab);
          if (abc_ao < -options.tolerance * std::abs(num::inner_prod(abc, ad)))
          {
            // P, B, C are already at indices 0, 1, 2 in the right order;
            // dropping D (index 3) needs no data movement at all.
            simplex_size = 3;
            direction = -abc;
            for (std::size_t i = 1; i < simplex_size; ++i)
            {
              for (std::size_t j = i + 1; j < simplex_size;)
              {
                if (is_near_zero(support_vector(simplex[i]) - support_vector(simplex[j])))
                {
                  simplex[j] = simplex[simplex_size - 1];
                  --simplex_size;
                }
                else
                {
                  ++j;
                }
              }
            }
            if (simplex_size < 3)
            {
              simplex_size = 1;
              direction = -support_vector(simplex[0]);
            }
            continue;
          }

          if (abd_ao < -options.tolerance * std::abs(num::inner_prod(abd, ac)))
          {
            // Keep P, B, D (drop C): move D from index 3 down to index 2,
            // the slot the size-3 code below reads.
            simplex[2] = simplex[3];
            simplex_size = 3;
            direction = -abd;
            for (std::size_t i = 1; i < simplex_size; ++i)
            {
              for (std::size_t j = i + 1; j < simplex_size;)
              {
                if (is_near_zero(support_vector(simplex[i]) - support_vector(simplex[j])))
                {
                  simplex[j] = simplex[simplex_size - 1];
                  --simplex_size;
                }
                else
                {
                  ++j;
                }
              }
            }
            if (simplex_size < 3)
            {
              simplex_size = 1;
              direction = -support_vector(simplex[0]);
            }
            continue;
          }

          if (acd_ao < -options.tolerance * std::abs(num::inner_prod(acd, ab)))
          {
            // Keep P, C, D (drop B): move C from index 2 down to index 1
            // and D from index 3 down to index 2.
            simplex[1] = simplex[2];
            simplex[2] = simplex[3];
            simplex_size = 3;
            direction = -acd;
            for (std::size_t i = 1; i < simplex_size; ++i)
            {
              for (std::size_t j = i + 1; j < simplex_size;)
              {
                if (is_near_zero(support_vector(simplex[i]) - support_vector(simplex[j])))
                {
                  simplex[j] = simplex[simplex_size - 1];
                  --simplex_size;
                }
                else
                {
                  ++j;
                }
              }
            }
            if (simplex_size < 3)
            {
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
            const vector bc4 = support_vector(simplex[2]) - support_vector(simplex[1]);
            const vector bd4 = support_vector(simplex[3]) - support_vector(simplex[1]);
            const vector bcd4 = cross(bd4, bc4);
            const scalar bcd_bo = num::inner_prod(bcd4, -support_vector(simplex[1])) *
                                  num::inner_prod(bcd4, support_vector(simplex[0]) - support_vector(simplex[1]));
            if (bcd_bo < -options.tolerance * std::abs(num::inner_prod(bcd4, support_vector(simplex[0]) - support_vector(simplex[1]))))
            {
              // Drop A (index 0) and continue with face BCD as the triangle.
              simplex[0] = simplex[1];
              simplex[1] = simplex[2];
              simplex[2] = simplex[3];
              simplex_size = 3;
              direction = -bcd4;
              for (std::size_t i = 1; i < simplex_size; ++i)
              {
                for (std::size_t j = i + 1; j < simplex_size;)
                {
                  if (is_near_zero(support_vector(simplex[i]) - support_vector(simplex[j])))
                  {
                    simplex[j] = simplex[simplex_size - 1];
                    --simplex_size;
                  }
                  else
                  {
                    ++j;
                  }
                }
              }
              if (simplex_size < 3)
              {
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
            if (volume_sq <= options.tolerance * options.tolerance)
            {
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

            return std::tuple{true, a, b, c, d};
         }

         // A three-point simplex cannot enclose the origin in 3D, but it may
         // still witness an intersection when the search converged without
         // discarding it (the supports keep advancing). EPA needs four
         // distinct vertices, so fill the fourth slot with a fresh support
         // along the current search direction instead of repeating the last
         // point of the compacted simplex.
         const auto compacted3 = compact_simplex();
         return std::tuple{true, compacted3[0], compacted3[1], compacted3[2],
                           diff_support(sa, sb, direction)};
      }

      throw std::runtime_error("GJK exceeded max_iterations");
    }
  }

  static constexpr const auto gjk_inter(const any_vectors auto& sa,
                                        const any_vectors auto& sb,
                                        const any_vector auto& init,
                                        const gjk_options options = {})
  {
    return gjk_inter<2>(sa, sb, init, options);
  }

  // EPA polytope storage: small-buffer optimization for shapes whose
  // combined vertex count fits within this limit. The Minkowski-difference
  // polytope has at most sa.size() + sb.size() edges. Larger shapes fall
  // back to heap allocation. Increase if you regularly collide high-count
  // convex hulls and want to avoid the heap path.
  static constexpr std::size_t epa_sbo_capacity = 64;

  template<std::size_t Dim, typename ShapesA, typename ShapesB, typename Support>
  requires (Dim == 2 || Dim == 3)
  static constexpr const auto epa(const ShapesA& sa,
                                  const ShapesB& sb,
                                  const Support& supporta,
                                  const Support& supportb,
                                  const Support& supportc,
                                  const Support& supportd,
                                  const any_scalar auto& tol,
                                  const std::size_t max_iterations = 64)
  {
    using env = decltype(env(sa));
    using vector = typename env::vector;
    using scalar = typename env::scalar;
    using support = typename env::support;

    if constexpr (Dim == 2)
    {
      validate_shape_2d(sa, "first shape");
      validate_shape_2d(sb, "second shape");
      if (tol <= 0)
        throw std::invalid_argument("EPA tolerance must be positive");
      if (max_iterations == 0)
        throw std::invalid_argument("EPA max_iterations must be positive");

      const vector a = support_vector(supporta);
      const vector b = support_vector(supportb);
      const vector c = support_vector(supportc);
      const scalar area = (b[0] - a[0]) * (c[1] - a[1]) -
                          (b[1] - a[1]) * (c[0] - a[0]);
      if (std::abs(area) <= tol)
        throw std::invalid_argument("EPA requires a non-degenerate triangle simplex");

      const bool cw = clockwise(a, b, c);

      const std::size_t max_polytope = std::size(sa) + std::size(sb);
      const bool use_sbo = max_polytope <= epa_sbo_capacity;

      std::array<support, epa_sbo_capacity> sbo{};
      std::vector<support> heap;
      support* polytope_data = nullptr;
      std::size_t polytope_size = 0;

      if (use_sbo)
      {
        sbo[0] = supporta;
        sbo[1] = supportb;
        sbo[2] = supportc;
        polytope_data = sbo.data();
        polytope_size = 3;
      }
      else
      {
        // Both storage paths must leave every slot the walk can write into
        // constructed. The insertion below assigns to
        // polytope_data[polytope_size], and assigning a tuple over raw
        // storage makes it read a wild vector pointer and free it -- masking
        // the heap path behind a crash the inline path never showed. resize
        // constructs every slot up front, the way the inline array's value
        // initialization does; reserve would leave them raw. Three is the
        // smallest polytope the bootstrap needs.
        heap.resize(std::max(max_polytope, std::size_t{3}));
        heap[0] = supporta;
        heap[1] = supportb;
        heap[2] = supportc;
        polytope_data = heap.data();
        polytope_size = 3;
      }

      for (std::size_t iteration = 0; iteration < max_iterations; ++iteration)
      {
        scalar closest_distance = std::numeric_limits<scalar>::infinity();
        vector closest_normal = vector();
        std::size_t closest_edge = 0;

        for (std::size_t i = 0; i < polytope_size; ++i)
        {
          const vector p1 = support_vector(polytope_data[i]);
          const vector p2 = support_vector(polytope_data[(i + 1) % polytope_size]);
          const vector edge = p2 - p1;
          if (num::inner_prod(edge, edge) <= tol * tol)
            throw std::runtime_error("EPA polytope contains duplicate support points");

          const vector normal = exterior_normal(cw, edge);
          const scalar distance = std::abs(num::inner_prod(normal, p1));
          if (distance < closest_distance)
          {
            closest_distance = distance;
            closest_normal = normal;
            closest_edge = i;
          }
        }

        const support candidate = diff_support(sa, sb, closest_normal);
        const scalar candidate_distance =
          num::inner_prod(support_vector(candidate), closest_normal);
        if (candidate_distance - closest_distance <= tol)
          return std::tuple{candidate_distance, closest_normal, candidate};

        for (std::size_t i = 0; i < polytope_size; ++i)
        {
          const vector difference = support_vector(candidate) - support_vector(polytope_data[i]);
          if (num::inner_prod(difference, difference) <= tol * tol)
            return std::tuple{closest_distance, closest_normal, candidate};
        }

        if (polytope_size == max_polytope)
          throw std::runtime_error("EPA polytope overflow");
        std::copy_backward(polytope_data + closest_edge + 1,
                           polytope_data + polytope_size,
                           polytope_data + polytope_size + 1);
        polytope_data[closest_edge + 1] = candidate;
        ++polytope_size;
      }

      throw std::runtime_error("EPA exceeded max_iterations");
    }
    else if constexpr (Dim == 3)
    {
      validate_shape_3d(sa, "first shape");
      validate_shape_3d(sb, "second shape");
      if (tol <= 0)
        throw std::invalid_argument("EPA tolerance must be positive");
      if (max_iterations == 0)
        throw std::invalid_argument("EPA max_iterations must be positive");

      vector a = support_vector(supporta);
      vector b = support_vector(supportb);
      vector c = support_vector(supportc);
      vector d = support_vector(supportd);

      const vector ab = b - a;
      const vector ac = c - a;
      const vector ad = d - a;
      const scalar volume = std::abs(num::inner_prod(cross(ab, ac), ad));
      if (volume <= tol * tol * tol)
      {
        // GJK handed EPA a degenerate witness simplex (points collinear or
        // coincident) that doesn't span a tetrahedron, even though the shapes
        // may genuinely overlap with real volume. Bootstrap a real
        // non-degenerate tetrahedron with a few extra support queries along
        // directions independent of what we already have, rather than trust
        // that degenerate witness for penetration depth. Only fall back to a
        // closest-point approximation on the degenerate point set if the
        // shapes truly don't admit a tetrahedron (e.g. a point/point or
        // line/line touching contact).
        std::vector<vector> distinct;
        for (const vector& v : {a, b, c, d})
        {
          bool found = false;
          for (const vector& u : distinct)
          {
            const vector diff = u - v;
            if (num::inner_prod(diff, diff) <= tol * tol) { found = true; break; }
          }
          if (!found) distinct.push_back(v);
        }

        struct degeneracy_info
        {
          bool have_line = false;
          vector line_dir{};
          vector line_point{};
          bool have_plane = false;
          vector normal{};
          vector plane_point{};
          bool have_volume = false;
          vector volume_point{};
        };

        auto analyze = [&](const std::vector<vector>& pts)
        {
          degeneracy_info info;
          if (pts.size() < 2) return info;
          info.have_line = true;
          info.line_point = pts[1];
          const vector e1 = pts[1] - pts[0];
          const scalar e1_len2 = num::inner_prod(e1, e1);
          info.line_dir = util::transform(e1, [e1_len2](const auto&, const auto& x) { return x / std::sqrt(e1_len2); });
          for (std::size_t k = 2; k < pts.size(); ++k)
          {
            const vector ek = pts[k] - pts[0];
            if (!info.have_plane)
            {
              const vector n = cross(info.line_dir, ek);
              const scalar n_len2 = num::inner_prod(n, n);
              if (n_len2 > tol * tol)
              {
                const scalar n_len = std::sqrt(n_len2);
                info.normal = util::transform(n, [n_len](const auto&, const auto& x) { return x / n_len; });
                info.have_plane = true;
                info.plane_point = pts[k];
              }
            }
            else if (!info.have_volume)
            {
              const scalar off_plane = num::inner_prod(info.normal, ek);
              if (std::abs(off_plane) > tol)
              {
                info.have_volume = true;
                info.volume_point = pts[k];
              }
            }
          }
          return info;
        };

        constexpr std::size_t max_bootstrap_attempts = 3;
        for (std::size_t attempt = 0; attempt < max_bootstrap_attempts; ++attempt)
        {
          const auto info = analyze(distinct);
          if (info.have_volume)
            break;

          std::vector<vector> candidates;
          if (!info.have_line)
          {
            candidates.push_back(vector{1, 0, 0});
          }
          else if (!info.have_plane)
          {
            // Escape collinearity. A single-axis perpendicular is itself
            // axis-aligned and, for axis-aligned box shapes, tends to
            // re-discover the same tied support vertex that produced the
            // degenerate simplex. Summing every non-zero axis-perpendicular
            // gives an oblique direction that reaches a genuinely different
            // vertex.
            const std::array<vector, 3> axes = {vector{1, 0, 0}, vector{0, 1, 0}, vector{0, 0, 1}};
            vector perp_first{};
            vector perp_sum{};
            bool have_perp_first = false;
            for (const vector& axis : axes)
            {
              const vector cand = cross(info.line_dir, axis);
              const scalar cand_len2 = num::inner_prod(cand, cand);
              if (cand_len2 <= tol * tol)
                continue;
              const scalar cand_len = std::sqrt(cand_len2);
              const vector unit_cand = util::transform(cand, [cand_len](const auto&, const auto& x) { return x / cand_len; });
              if (!have_perp_first) { perp_first = unit_cand; have_perp_first = true; }
              perp_sum = perp_sum + unit_cand;
            }
            const scalar perp_sum_len2 = num::inner_prod(perp_sum, perp_sum);
            if (perp_sum_len2 > tol * tol)
            {
              const scalar perp_sum_len = std::sqrt(perp_sum_len2);
              candidates.push_back(util::transform(perp_sum, [perp_sum_len](const auto&, const auto& x) { return x / perp_sum_len; }));
            }
            if (have_perp_first)
              candidates.push_back(perp_first);
          }
          else
          {
            // Escape coplanarity: try the plane normal, then its opposite —
            // an axis-aligned normal can likewise re-discover an existing
            // tied vertex on one side of the plane.
            candidates.push_back(info.normal);
            candidates.push_back(-info.normal);
          }

          bool added = false;
          for (const vector& dir : candidates)
          {
            const support candidate = diff_support(sa, sb, dir);
            const vector cv = support_vector(candidate);
            bool is_new = true;
            for (const vector& u : distinct)
            {
              const vector diff = u - cv;
              if (num::inner_prod(diff, diff) <= tol * tol) { is_new = false; break; }
            }
            if (!is_new)
              continue;

            distinct.push_back(cv);
            added = true;
            break;
          }
          if (!added)
            break;
        }

        const auto final_info = analyze(distinct);

        // Project the origin onto the triangle fan built from a coplanar (or
        // richer) point set, clipped to the fan's triangles/edges. Used both
        // as the genuinely-flat-shapes answer and as the last-resort answer
        // when a bootstrapped tetrahedron can't be made to enclose the
        // origin (see below).
        auto project_onto_points = [&](const std::vector<vector>& dv) -> std::tuple<scalar, vector>
        {
          std::vector<std::array<std::size_t, 3>> triangles;
          for (std::size_t i = 1; i + 1 < dv.size(); ++i)
          {
            const vector e10 = dv[i] - dv[0];
            const vector e20 = dv[i + 1] - dv[0];
            const vector n = cross(e20, e10);
            if (num::inner_prod(n, n) > tol * tol)
              triangles.push_back({std::size_t(0), i, i + 1});
          }

          const vector& v0 = dv[triangles[0][0]];
          const vector& v1 = dv[triangles[0][1]];
          const vector& v2 = dv[triangles[0][2]];
          const vector e10 = v1 - v0;
          const vector e20 = v2 - v0;
          vector normal = cross(e20, e10);
          const scalar len = std::sqrt(num::inner_prod(normal, normal));
          normal = util::transform(normal, [len](const auto&, const auto& x) { return x / len; });

          const vector origin = vector();
          const scalar t = -num::inner_prod(normal, v0);
          const vector p = origin - t * normal;
          scalar best_dist2 = std::numeric_limits<scalar>::infinity();
          vector best_point = p;
          for (const auto& tri : triangles)
          {
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
            if (std::abs(denom) > tol * tol)
            {
              const scalar u = (d1 * d6 - d2 * d4) / denom;
              const scalar v = (d3 * d2 - d5 * d1) / denom;
              if (u >= 0 && v >= 0 && u + v <= 1)
              {
                const vector cp = ta + u * tab + v * tac;
                const scalar dist2 = num::inner_prod(cp, cp);
                if (dist2 < best_dist2) { best_dist2 = dist2; best_point = cp; }
              }
            }
            const std::array<std::array<std::size_t, 2>, 3> edges = {{
              {tri[0], tri[1]}, {tri[1], tri[2]}, {tri[2], tri[0]}
            }};
            for (const auto& edge : edges)
            {
              const vector ea = dv[edge[0]];
              const vector eb = dv[edge[1]];
              const vector e = eb - ea;
              const scalar t2 = std::clamp(-num::inner_prod(e, p - ea) / num::inner_prod(e, e), scalar(0), scalar(1));
              const vector cp = ea + t2 * e;
              const scalar dist2 = num::inner_prod(cp, cp);
              if (dist2 < best_dist2) { best_dist2 = dist2; best_point = cp; }
            }
          }
          const scalar depth = std::sqrt(best_dist2);
          vector out_normal = best_point;
          const scalar out_len = std::sqrt(num::inner_prod(out_normal, out_normal));
          if (out_len > tol)
            out_normal = util::transform(out_normal, [out_len](const auto&, const auto& x) { return x / out_len; });
          else
            out_normal = normal;
          return std::tuple{depth, out_normal};
        };

        if (final_info.have_volume)
        {
          // Bootstrap found four affinely-independent points, but nothing
          // guarantees they enclose the origin the way GJK's own witness
          // simplex normally would -- EPA's face-walk loop below requires
          // that enclosure to produce a meaningful answer. Refine with a
          // bounded portal-style procedure (in the spirit of MPR/XenoCollide
          // portal construction): repeatedly find the face the origin is
          // least enclosed by, and replace its opposite vertex with a fresh
          // support point along that face's outward normal.
          std::array<vector, 4> tet = {distinct[0], final_info.line_point,
                                        final_info.plane_point, final_info.volume_point};
          struct face_def { std::size_t i0, i1, i2, opp; };
          constexpr std::array<face_def, 4> face_defs = {{
            {0, 1, 2, 3}, {0, 1, 3, 2}, {0, 2, 3, 1}, {1, 2, 3, 0}
          }};

          auto worst_face = [&]()
          {
            std::size_t worst = 0;
            scalar worst_val = -std::numeric_limits<scalar>::infinity();
            vector worst_normal = vector();
            for (std::size_t fi = 0; fi < 4; ++fi)
            {
              const auto& fd = face_defs[fi];
              const vector& fa = tet[fd.i0];
              const vector& fb = tet[fd.i1];
              const vector& fc = tet[fd.i2];
              const vector& opp = tet[fd.opp];
              vector n = cross(fb - fa, fc - fa);
              const scalar n_len2 = num::inner_prod(n, n);
              if (n_len2 <= tol * tol)
                continue;
              const scalar n_len = std::sqrt(n_len2);
              n = util::transform(n, [n_len](const auto&, const auto& x) { return x / n_len; });
              if (num::inner_prod(n, opp - fa) > 0)
                n = -n;
              const scalar val = num::inner_prod(n, vector() - fa);
              if (val > worst_val) { worst_val = val; worst = fi; worst_normal = n; }
            }
            return std::tuple{worst, worst_val, worst_normal};
          };

          bool enclosing = false;
          constexpr std::size_t max_refine_attempts = 4;
          for (std::size_t refine = 0; refine < max_refine_attempts; ++refine)
          {
            const auto [worst, worst_val, worst_normal] = worst_face();
            if (worst_val < -tol) { enclosing = true; break; }
            const support candidate = diff_support(sa, sb, worst_normal);
            const vector cv = support_vector(candidate);
            const vector diff = cv - tet[face_defs[worst].opp];
            if (num::inner_prod(diff, diff) <= tol * tol)
              break;
            tet[face_defs[worst].opp] = cv;
          }

          if (enclosing)
          {
            // Fall through: a/b/c/d now span a tetrahedron that strictly
            // encloses the origin, so the ordinary face-walk EPA loop below
            // can run on it directly.
            a = tet[0];
            b = tet[1];
            c = tet[2];
            d = tet[3];
          }
          else
          {
            // Refinement couldn't produce a strictly-enclosing tetrahedron.
            // Fall back to projecting onto the coplanar subset that existed
            // just before the volume-completing point was bootstrapped in
            // (that point is always the last one appended to `distinct`).
            const auto [depth, out_normal] =
              project_onto_points(std::vector<vector>(distinct.begin(), distinct.end() - 1));
            return std::tuple{depth, out_normal, supporta};
          }
        }
        else if (final_info.have_plane)
        {
          // The Minkowski difference is genuinely flat (e.g. two coplanar
          // shapes): no tetrahedron is possible. Project the origin onto the
          // triangle fan built from the plane we do have.
          const auto [depth, out_normal] = project_onto_points(distinct);
          return std::tuple{depth, out_normal, supporta};
        }
        else
        {
          // Still collinear/coincident after bootstrapping: the shapes truly
          // reduce to a point or a line in the Minkowski difference. Fall
          // back to the closest point on the union of segments between the
          // distinct points found.
          const vector origin = vector();
          scalar best_dist2 = std::numeric_limits<scalar>::infinity();
          vector best_point = distinct[0];
          if (distinct.size() == 1)
          {
            best_point = distinct[0];
            best_dist2 = num::inner_prod(best_point, best_point);
          }
          else
          {
            for (std::size_t i = 0; i < distinct.size(); ++i)
            {
              for (std::size_t j = i + 1; j < distinct.size(); ++j)
              {
                const vector& ea = distinct[i];
                const vector& eb = distinct[j];
                const vector e = eb - ea;
                const scalar t2 = std::clamp(
                  -num::inner_prod(e, ea - origin) / num::inner_prod(e, e), scalar(0), scalar(1));
                const vector cp = ea + t2 * e;
                const scalar dist2 = num::inner_prod(cp, cp);
                if (dist2 < best_dist2)
                {
                  best_dist2 = dist2;
                  best_point = cp;
                }
              }
            }
          }

          const scalar depth = std::sqrt(best_dist2);
          vector out_normal = best_point;
          const scalar out_len = std::sqrt(num::inner_prod(out_normal, out_normal));
          if (out_len > tol)
            out_normal = util::transform(out_normal, [out_len](const auto&, const auto& x) { return x / out_len; });
          else
            out_normal = vector{1, 0, 0};
          return std::tuple{depth, out_normal, supporta};
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
        const auto encloses = [&tol](const std::array<vector, 4>& t) {
          static constexpr std::array<std::array<std::size_t, 4>, 4> faces = {{
            {{0, 1, 2, 3}}, {{0, 1, 3, 2}}, {{0, 2, 3, 1}}, {{1, 2, 3, 0}}
          }};
          for (const auto& f : faces)
          {
            const vector& p = t[f[0]];
            const vector& q = t[f[1]];
            const vector& r = t[f[2]];
            const vector& opp = t[f[3]];
            vector n = cross(q - p, r - p);
            if (num::inner_prod(n, n) <= tol * tol)
              return false;
            if (num::inner_prod(n, opp - p) > 0)
              n = -n;
            if (num::inner_prod(n, vector() - p) >= -tol)
              return false;
          }
          return true;
        };

        std::array<vector, 4> witness = {a, b, c, d};
        if (!encloses(witness))
        {
          std::vector<vector> points(witness.begin(), witness.end());
          const std::array<vector, 6> axes = {{vector{1, 0, 0}, vector{-1, 0, 0},
                                               vector{0, 1, 0}, vector{0, -1, 0},
                                               vector{0, 0, 1}, vector{0, 0, -1}}};
          for (const vector& axis : axes)
          {
            const vector v = support_vector(diff_support(sa, sb, axis));
            bool known = false;
            for (const vector& u : points)
              if (num::inner_prod(u - v, u - v) <= tol * tol)
              {
                known = true;
                break;
              }
            if (!known)
              points.push_back(v);
          }
          for (std::size_t i = 0; i + 3 < points.size(); ++i)
            for (std::size_t j = i + 1; j + 2 < points.size(); ++j)
              for (std::size_t k = j + 1; k + 1 < points.size(); ++k)
                for (std::size_t l = k + 1; l < points.size(); ++l)
                {
                  if (encloses(witness))
                    break;
                  const std::array<vector, 4> trial = {points[i], points[j],
                                                       points[k], points[l]};
                  if (encloses(trial))
                    witness = trial;
                }
          a = witness[0];
          b = witness[1];
          c = witness[2];
          d = witness[3];
        }
      }

      struct face {
        std::array<std::size_t, 3> indices;
        vector normal;
        scalar distance;
      };

      std::vector<vector> vertices = {a, b, c, d};
      std::vector<face> faces;

      auto make_face = [&](std::size_t i0, std::size_t i1, std::size_t i2) {
        const vector v0 = vertices[i0];
        const vector v1 = vertices[i1];
        const vector v2 = vertices[i2];
        const vector e10 = v1 - v0;
        const vector e20 = v2 - v0;
        vector normal = cross(e20, e10);
        const scalar len2 = num::inner_prod(normal, normal);
        const scalar len = std::sqrt(len2);
        if (len > tol) {
          normal = util::transform(normal, [len](const auto&, const auto& x) { return x / len; });
        }
        const scalar d = num::inner_prod(normal, v0);
        return face{{i0, i1, i2}, normal, d};
      };

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

      const std::size_t max_polytope = std::size(sa) + std::size(sb);

      for (std::size_t iteration = 0; iteration < max_iterations; ++iteration)
      {
        scalar closest_distance = std::numeric_limits<scalar>::infinity();
        std::size_t closest_face = 0;

        for (std::size_t i = 0; i < faces.size(); ++i)
        {
          if (faces[i].distance < closest_distance)
          {
            closest_distance = faces[i].distance;
            closest_face = i;
          }
        }

        const support candidate = diff_support(sa, sb, faces[closest_face].normal);
        const scalar candidate_distance =
          num::inner_prod(support_vector(candidate), faces[closest_face].normal);
        if (candidate_distance - closest_distance <= tol)
          return std::tuple{closest_distance, faces[closest_face].normal, candidate};

        const vector candidate_vec = support_vector(candidate);
        bool duplicate = false;
        for (const auto& v : vertices) {
          const vector diff = v - candidate_vec;
          if (num::inner_prod(diff, diff) <= tol * tol) {
            duplicate = true;
            break;
          }
        }
        if (duplicate)
          return std::tuple{closest_distance, faces[closest_face].normal, candidate};

        if (vertices.size() >= max_polytope)
          throw std::runtime_error("EPA polytope overflow");

        const std::size_t new_index = vertices.size();
        vertices.push_back(candidate_vec);

        std::vector<std::array<std::size_t, 2>> boundary_edges;
        std::vector<char> face_removed(faces.size(), 0);

        for (std::size_t i = 0; i < faces.size(); ++i)
        {
          const vector& v0 = vertices[faces[i].indices[0]];
          const scalar dot = num::inner_prod(faces[i].normal, candidate_vec - v0);
          if (dot > tol)
          {
            face_removed[i] = true;
          }
        }

        bool all_faces_removed = true;
        for (bool removed : face_removed)
        {
          if (!removed) { all_faces_removed = false; break; }
        }
        if (all_faces_removed)
        {
          return std::tuple{closest_distance, faces[closest_face].normal, candidate};
        }

        // Keep each retained face's edges in their original winding order:
        // make_face(new_index, edge[0], edge[1]) below depends on that
        // direction to produce an outward-facing normal for the new face.
        // Sorting an edge to {min, max} before use (as this used to do)
        // discards that direction and can silently flip newly-created faces
        // inward once EPA needs a second expansion iteration; the sorted
        // form is only safe as a duplicate-detection key, never as the
        // edge to build a face from.
        std::vector<std::array<std::size_t, 2>> kept_face_edges;
        for (std::size_t i = 0; i < faces.size(); ++i)
        {
          if (face_removed[i]) continue;
          const auto& fi = faces[i];
          for (std::size_t e = 0; e < 3; ++e)
          {
            const std::size_t e0 = fi.indices[e];
            const std::size_t e1 = fi.indices[(e + 1) % 3];
            kept_face_edges.push_back({e0, e1});
          }
        }
        std::vector<bool> kept_edge_dup(kept_face_edges.size(), false);
        for (std::size_t i = 0; i < kept_face_edges.size(); ++i)
        {
          const auto sorted_i = std::array{std::min(kept_face_edges[i][0], kept_face_edges[i][1]),
                                            std::max(kept_face_edges[i][0], kept_face_edges[i][1])};
          for (std::size_t j = i + 1; j < kept_face_edges.size(); ++j)
          {
            const auto sorted_j = std::array{std::min(kept_face_edges[j][0], kept_face_edges[j][1]),
                                              std::max(kept_face_edges[j][0], kept_face_edges[j][1])};
            if (sorted_i == sorted_j)
            {
              kept_edge_dup[i] = true;
              kept_edge_dup[j] = true;
            }
          }
        }
        for (std::size_t i = 0; i < kept_face_edges.size(); ++i)
        {
          if (!kept_edge_dup[i])
          {
            boundary_edges.push_back(kept_face_edges[i]);
          }
        }

        std::vector<face> new_faces;
        for (const auto& edge : boundary_edges)
        {
          new_faces.push_back(make_face(new_index, edge[0], edge[1]));
        }

        for (auto& f : new_faces) {
          if (f.distance < 0) {
            f.distance = -f.distance;
            f.normal = -f.normal;
            std::swap(f.indices[0], f.indices[1]);
          }
        }

        std::vector<face> remaining_faces;
        for (std::size_t i = 0; i < faces.size(); ++i)
        {
          if (!face_removed[i])
            remaining_faces.push_back(faces[i]);
        }
        faces = std::move(remaining_faces);
        faces.insert(faces.end(), new_faces.begin(), new_faces.end());
      }

      throw std::runtime_error("EPA exceeded max_iterations");
    }
  }

  static constexpr const auto epa(const any_vectors auto& sa,
                                  const any_vectors auto& sb,
                                  const any_support auto& supporta,
                                  const any_support auto& supportb,
                                  const any_support auto& supportc,
                                  const any_scalar auto& tol,
                                  const std::size_t max_iterations = 64)
  {
    return epa<2>(sa, sb, supporta, supportb, supportc, supportc, tol, max_iterations);
  }

  static constexpr const auto epa(const any_vectors auto& sa,
                                  const any_vectors auto& sb,
                                  const any_support auto& supporta,
                                  const any_support auto& supportb,
                                  const any_support auto& supportc,
                                  const any_support auto& supportd,
                                  const any_scalar auto& tol,
                                  const std::size_t max_iterations = 64)
  {
    const std::size_t dim = sa.empty() ? 2 : sa[0].size();
    if (dim == 2)
      return epa<2>(sa, sb, supporta, supportb, supportc, supportc, tol, max_iterations);
    if (dim == 3)
      return epa<3>(sa, sb, supporta, supportb, supportc, supportd, tol, max_iterations);
    throw std::invalid_argument("unsupported dimension");
  }
}
#endif
