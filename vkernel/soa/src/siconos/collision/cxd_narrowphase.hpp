#pragma once

/* cxd-backed narrowphase batch engine (pure functions, no vkernel item).

   Semantics (spec §4, exact-value tests in tests/polytope_cxd.cpp):
   - h follows the siconos gap convention: +distance when separated,
     -depth when penetrating.
   - the contact normal points 1 -> 2; equivalently it is the exterior
     normal of a face of the Minkowski-difference polytope.
     separated (origin outside): distance GJK gives the witness foot
     f = a - b, normal = -f/|f|; the witness pair is the supporting-FACE
     construction along the normal (overlap-window centers for flush
     edge-edge contacts, vertices / perpendicular feet otherwise), so the
     cps and h stay mutually consistent: h = n.(cp2 - cp1).
     enclosed (origin inside): cxd's EPA gives the depth and the normal
     of the deepest face of the difference polytope -- its 2D EPA
     bootstraps a non-degenerate seed when the GJK hands a degenerate one,
     which a flush axis-aligned overlap always is; the witnesses are the
     supporting faces along that normal, and h = -depth.
   - cp1/cp2 for both cases are the supporting faces' points (uniquely
     determined, deterministic); cxd's support tuples stay inside the
     GJK internals.

   Separation sub-distance: computed siconos-side with a small 2D distance
   GJK over the pool, using cxd's own support primitives (shape_view +
   cxd::diff_support). cxd's batch deliberately leaves the seed EMPTY on
   separated statuses ("misses keep every slot untouched"), so the
   seed-projection sub-distance of the spec has to travel with the engine.

   max_iterations / invalid_input: h = 0 placeholder + status flag; the
   CALLER (the stage) keeps the relation's previous cached row (status, h,
   normal, cps). A pair whose separation sub-distance gives up, or whose
   EPA depth is non-finite, is reported with such a failure status -- never
   as a successful separated/enclosed row with h = 0. No output is ever
   non-finite. */

#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <vector>
// clang-format off
/* load-bearing order — mirrors cxd's own test harness (test.cpp); do not
   regroup or sort: backends must be registered before cxd_batch.hpp pulls
   in cxd.hpp, which resolves cxd::num::inner_prod through them. */
#include "cxd_std.hpp"
#include "cxd_std_array_vec.hpp"
#include "cxd_ublas.hpp"
#include "cxd_eigen.hpp"
#include "cxd_batch.hpp"
// clang-format on

namespace siconos::collision {

/* The batch engine is scalar-generic: the stage instantiates it with the
   environment's scalar (some::scalar cannot be written in a non-item
   struct -- it is a storage descriptor, not a type), and cxd's packed pool
   needs the point type to be a trivially copyable std::array. */
template <typename Scalar>
struct narrowphase_point {
  Scalar h;  // gap of this point along the contact normal
  std::array<Scalar, 2> cp1, cp2;
};

template <typename Scalar>
struct narrowphase_contact {
  int status;                    // cxd::status value
  Scalar h;                      // siconos gap: +distance / -depth
  std::array<Scalar, 2> normal;  // contact normal, body 1 -> body 2
  std::array<Scalar, 2> cp1, cp2;
  // contact manifold: up to two points with their own gap, from
  // reference / incident face clipping, sorted along the tangent
  // (-normal[1], normal[0]). A contact whose closest features are not a
  // face (vertex against vertex) has the single point (h, cp1, cp2).
  std::size_t npoints;
  std::array<narrowphase_point<Scalar>, 2> points;
};

template <typename Scalar>
struct narrowphase_input {
  std::vector<std::array<Scalar, 2>> vertices;  // shape vertices concatenated, world frame
  std::vector<cxd::shape_entry> entries;        // one per shape id: offset/count
  std::vector<std::uint32_t> ids_a,
      ids_b;                                     // shape ids per pair (id = entries index)
  std::vector<std::array<Scalar, 2>> init_dirs;  // zero -> {1,0}
};

namespace detail {

template <typename V>
using support_t = std::tuple<std::size_t, std::size_t, V>;

template <typename V>
constexpr auto vdot(V a, V b) {
  using scalar = typename V::value_type;
  return static_cast<scalar>(a[0] * b[0] + a[1] * b[1]);
}
template <typename V>
inline auto vnorm(V a) {
  return std::hypot(a[0], a[1]);
}
template <typename V>
inline V vadd(V a, V b) {
  return {a[0] + b[0], a[1] + b[1]};
}
template <typename V>
inline V vscale(typename V::value_type s, V a) {
  return {s * a[0], s * a[1]};
}

template <typename V>
inline const V& support_pos(const support_t<V>& s) {
  return std::get<2>(s);
}

template <typename V>
inline V vsub2(V a, V b) {
  return {a[0] - b[0], a[1] - b[1]};
}

// closest point of the origin to the segment [p, q]: foot and t in [0,1]
// (t == 0 / 1 picks the vertex feature)
template <typename V>
inline std::pair<V, typename V::value_type> foot_on_segment(V p, V q) {
  using scalar = typename V::value_type;
  const V e = vsub2(q, p);
  const scalar ee = vdot(e, e);
  if (ee <= scalar{0}) return {p, scalar{0}};
  scalar t = -vdot(p, e) / ee;
  t = std::max(scalar{0}, std::min(scalar{1}, t));
  return {vadd(p, vscale(t, e)), t};
}

/* One reduction step: closest point of the origin to the convex hull of
   the (up to 3) supports; returns the foot, per-support weights (same
   count as the given simplex range), and a flag whether the foot is on a
   VERTEX (so the simplex reduces to that vertex) or on an EDGE (reduces
   to that edge's two supports). */
template <typename V>
inline auto closest_to_origin(const std::array<support_t<V>, 3>& w, std::size_t n)
    -> std::tuple<V, std::array<typename V::value_type, 3>, bool /*vertex feature*/,
                  std::size_t /*feature index*/> {
  using scalar = typename V::value_type;
  std::array<scalar, 3> alpha{};
  V best_foot{};
  scalar best_norm = std::numeric_limits<scalar>::infinity();
  bool best_is_vertex = false;
  std::size_t best_index = 0;

  // vertex features
  for (std::size_t i = 0; i < n; ++i) {
    const auto& p = support_pos(w[i]);
    if (const scalar d = vnorm(p); d < best_norm) {
      std::array<scalar, 3> a{};
      a[i] = scalar{1};
      best_norm = d;
      best_foot = p;
      best_is_vertex = true;
      best_index = i;
      alpha = a;
    }
  }
  // edge features
  for (std::size_t i = 0; i < n; ++i) {
    const auto& p = support_pos(w[i]);
    const auto& q = support_pos(w[(i + 1) % n]);
    const auto [foot, t] = foot_on_segment(p, q);
    if (const scalar d = vnorm(foot); d < best_norm) {
      std::array<scalar, 3> a{};
      a[i] = scalar{1} - t;
      a[(i + 1) % n] = t;
      best_norm = d;
      best_foot = foot;
      best_is_vertex = (t == scalar{0} || t == scalar{1});  // degenerate to a vertex
      best_index = i;
      alpha = a;
    }
  }
  return {best_foot, alpha, best_is_vertex, best_index};
}

/* 2D distance GJK over one pooled pair; exact separation witness.
   Returns nullopt when no separation certificate was reached (touching,
   enclosed or budget exhausted): the caller then reports the pair as not
   separated and keeps its previous cached h. */
template <typename Pool>
inline auto subdistance_gjk(const Pool& pool, std::uint32_t id_a, std::uint32_t id_b,
                            const typename Pool::vector& dir0, std::size_t budget)
    -> std::optional<
        std::tuple<typename Pool::scalar /*distance*/, typename Pool::vector /*foot*/,
                   std::array<support_t<typename Pool::vector>, 3>, std::size_t /*n*/>> {
  using V = typename Pool::vector;
  using scalar = typename Pool::scalar;
  const auto ds = [&](V d) -> support_t<V> {
    return cxd::diff_support(cxd::shape_view<Pool>{pool, id_a},
                             cxd::shape_view<Pool>{pool, id_b}, d);
  };

  // float environments need a tolerance above their epsilon; double keeps
  // the spec's 1e-10
  const scalar tol =
      std::max(scalar{1e-10}, scalar{8} * std::numeric_limits<scalar>::epsilon());

  std::array<support_t<V>, 3> w{};
  w[0] = ds(dir0);
  std::size_t n = 1;

  for (std::size_t it = 0; it < budget; ++it) {
    const auto [foot, alpha, is_vertex, fidx] = closest_to_origin(w, n);

    // origin on/inside the hull: no separation certificate
    if (vnorm(foot) <= tol) return std::nullopt;

    const V dir = {-foot[0], -foot[1]};  // toward the origin

    // simplex reduction to the current feature (the edge (i, (i+1)%n))
    if (is_vertex) {
      w[0] = w[fidx];
      n = 1;
    } else {
      const support_t<V> wa = w[fidx];
      const support_t<V> wb = w[(fidx + 1) % n];
      w[0] = wa;
      w[1] = wb;
      n = 2;
    }

    const support_t<V> cand = ds(dir);
    // no support strictly beyond the current foot: exact sub-distance
    // (dir = -foot, so the support bound and the simplex bound coincide)
    if (vdot(support_pos(cand), dir) <= -vdot(foot, foot) + tol) {
      scalar distance = vnorm(foot);
      return std::tuple{distance, foot, w, n};
    }

    // insert the candidate at the front of the reduced simplex
    if (n == 3) {
      w[2] = w[1];
      w[1] = w[0];
      w[0] = cand;
    } else {
      w[n] = w[0];
      w[0] = cand;
      ++n;
    }
  }
  return std::nullopt;  // budget exhausted
}

/* Witness pair for a separated pair, from the SUPPORTING FACES of both
   polygons along the contact normal n (unit, 1 -> 2). cxd's support query
   tie-breaks to the lowest vertex index, which would report a corner of a
   flush edge-edge contact; the face-based witness reports the midpoint of
   the overlap window instead (the pinned exact values, e.g. cp (1,0)/(2,0)
   for the two-squares gap-1 case). Vertex features fall out naturally:
   a unique supporting vertex gives that vertex, a vertex-edge pair gives
   the vertex and the foot of the perpendicular. */
template <typename Pool>
inline auto face_witnesses(const Pool& pool, std::uint32_t id_a, std::uint32_t id_b,
                           const typename Pool::vector& n)
    -> std::pair<typename Pool::vector, typename Pool::vector> {
  using V = typename Pool::vector;
  using scalar = typename Pool::scalar;
  const cxd::shape_view<Pool> va{pool, id_a};
  const cxd::shape_view<Pool> vb{pool, id_b};

  const auto nverts = [](const auto& view) {
    return static_cast<std::size_t>(std::distance(view.begin(), view.end()));
  };
  const scalar amax = [&] {
    scalar m = -std::numeric_limits<scalar>::infinity();
    for (std::size_t i = 0; i < nverts(va); ++i) m = std::max(m, vdot(n, va[i]));
    return m;
  }();
  const scalar bmin = [&] {
    scalar m = std::numeric_limits<scalar>::infinity();
    for (std::size_t j = 0; j < nverts(vb); ++j) m = std::min(m, vdot(n, vb[j]));
    return m;
  }();
  const scalar eps = scalar{1e-9} * (scalar{1} + std::abs(amax) + std::abs(bmin));

  std::vector<std::size_t> tie_a, tie_b;
  for (std::size_t i = 0; i < nverts(va); ++i)
    if (vdot(n, va[i]) >= amax - eps) tie_a.push_back(i);
  for (std::size_t j = 0; j < nverts(vb); ++j)
    if (vdot(n, vb[j]) <= bmin + eps) tie_b.push_back(j);

  if (tie_a.size() == 1 && tie_b.size() == 1) return {va[tie_a[0]], vb[tie_b[0]]};

  // the supporting face of a tie set: its two extreme vertices along the
  // tangent. Tied vertices are not consecutive in index order when the face
  // wraps around index 0, so front()/back() of the tie set is not the face.
  const V tng{-n[1], n[0]};
  const auto face_ends = [&](const auto& view, const std::vector<std::size_t>& tie) {
    std::size_t lo = tie.front(), hi = tie.front();
    for (auto i : tie) {
      if (vdot(tng, view[i]) < vdot(tng, view[lo])) lo = i;
      if (vdot(tng, view[i]) > vdot(tng, view[hi])) hi = i;
    }
    return std::pair{lo, hi};
  };
  // closest point of the face [a, b] to the vertex p
  const auto foot_of = [](const V& p, const V& a, const V& b) {
    const auto [foot, t] = foot_on_segment(vsub2(a, p), vsub2(b, p));
    return vadd(foot, p);
  };

  if (tie_a.size() == 1) {  // vertex on A, face on B
    const V p = va[tie_a[0]];
    const auto [b0, b1] = face_ends(vb, tie_b);
    return {p, foot_of(p, vb[b0], vb[b1])};
  }
  if (tie_b.size() == 1) {  // face on A, vertex on B
    const V p = vb[tie_b[0]];
    const auto [a0, a1] = face_ends(va, tie_a);
    return {foot_of(p, va[a0], va[a1]), p};
  }

  // both faces (parallel edges): midpoint of the overlap window along the
  // tangent
  const auto [ia0, ia1] = face_ends(va, tie_a);
  const auto [ib0, ib1] = face_ends(vb, tie_b);
  const scalar tl = vnorm(tng);
  const V e = tl > scalar{0} ? vscale(scalar{1} / tl, tng) : V{1, 0};
  const scalar pa0 = vdot(e, va[ia0]);
  const scalar pa1 = vdot(e, va[ia1]);
  const scalar pb0 = vdot(e, vb[ib0]);
  const scalar pb1 = vdot(e, vb[ib1]);
  const scalar lo = std::max(std::min(pa0, pa1), std::min(pb0, pb1));
  const scalar hi = std::min(std::max(pa0, pa1), std::max(pb0, pb1));
  const scalar tm = scalar{0.5} * (lo + hi);
  const auto lerp_on = [&](const auto& view, std::size_t i0, std::size_t i1, scalar p0,
                           scalar p1) -> V {
    const scalar denom = p1 - p0;
    if (std::abs(denom) < scalar{1e-30}) return view[i0];
    const scalar t = std::max(scalar{0}, std::min(scalar{1}, (tm - p0) / denom));
    return vadd(view[i0], vscale(t, vsub2(view[i1], view[i0])));
  };
  return {lerp_on(va, ia0, ia1, pa0, pa1), lerp_on(vb, ib0, ib1, pb0, pb1)};
}

/* separating-axis contact of two convex vertex sets (a segment is a
   2-vertex set): the edge normal of either shape with the largest
   separation, oriented from A to B. The separation is exact when the shapes
   touch or penetrate, and a lower bound of the distance when they are apart
   (vertex-vertex closest features). Fallback for the pairs the GJK / EPA
   path cannot answer -- in practice exactly touching shapes -- so that no
   pair is left with a zero normal. */
template <typename Pool>
inline auto sat_contact(const Pool& pool, std::uint32_t id_a, std::uint32_t id_b)
    -> std::optional<std::pair<typename Pool::scalar, typename Pool::vector>> {
  using V = typename Pool::vector;
  using scalar = typename Pool::scalar;
  const cxd::shape_view<Pool> va{pool, id_a};
  const cxd::shape_view<Pool> vb{pool, id_b};
  const auto nverts = [](const auto& view) {
    return static_cast<std::size_t>(std::distance(view.begin(), view.end()));
  };
  const auto range = [&](const auto& view, const V& u) {
    scalar lo = std::numeric_limits<scalar>::infinity();
    scalar hi = -std::numeric_limits<scalar>::infinity();
    for (std::size_t i = 0; i < nverts(view); ++i) {
      lo = std::min(lo, vdot(u, view[i]));
      hi = std::max(hi, vdot(u, view[i]));
    }
    return std::pair{lo, hi};
  };
  bool have = false;
  scalar best = scalar{0};
  V best_n{0, 0};
  const auto test_edges = [&](const auto& view) {
    const std::size_t n = nverts(view);
    for (std::size_t i = 0; i < n; ++i) {
      const V e = vsub2(view[(i + 1) % n], view[i]);
      const scalar len = vnorm(e);
      if (len <= scalar{0}) continue;
      const V u{e[1] / len, -e[0] / len};
      const auto [amin, amax] = range(va, u);
      const auto [bmin, bmax] = range(vb, u);
      // separation along +u (A then B) and along -u (B then A)
      const scalar s_pos = bmin - amax;
      const scalar s_neg = amin - bmax;
      const scalar sep = std::max(s_pos, s_neg);
      if (!have || sep > best) {
        have = true;
        best = sep;
        best_n = s_pos >= s_neg ? u : V{-u[0], -u[1]};
      }
    }
  };
  test_edges(va);
  test_edges(vb);
  if (!have) return std::nullopt;
  return std::pair{best, best_n};
}

/* Reference / incident face clipping (Box2D b2CollidePolygons style): the
   contact manifold of a face contact along the normal n (A -> B).

   On each shape the face whose outward normal is closest to +n (on A) or to
   -n (on B); the more parallel one is the reference face, the other the
   incident face (ties go to A, so the choice does not flip between equal
   faces). The incident edge is clipped to the extent of the reference face:
   at most two points, each with its own gap measured along the reference
   normal, and its witness on the reference face. A segment (2 vertices) is a
   single face with both normals. When n is not a face normal of either shape
   (the closest features are two vertices) there is no face contact and 0 is
   returned. */
template <typename Pool, typename Scalar>
inline std::size_t clip_manifold(const Pool& pool, std::uint32_t id_a, std::uint32_t id_b,
                                 const typename Pool::vector& n,
                                 std::array<narrowphase_point<Scalar>, 2>& out) {
  using V = typename Pool::vector;
  using scalar = typename Pool::scalar;
  const cxd::shape_view<Pool> va{pool, id_a};
  const cxd::shape_view<Pool> vb{pool, id_b};
  const auto nverts = [](const auto& view) {
    return static_cast<std::size_t>(std::distance(view.begin(), view.end()));
  };

  struct face {
    std::size_t i0, i1;
    V normal;
    scalar align;
  };
  const auto best_face = [&](const auto& view, const V& d) {
    const std::size_t nv = nverts(view);
    scalar area2 = scalar{0};
    for (std::size_t i = 0; i < nv; ++i) {
      const V& a = view[i];
      const V& b = view[(i + 1) % nv];
      area2 += a[0] * b[1] - b[0] * a[1];
    }
    const scalar orient = area2 >= scalar{0} ? scalar{1} : scalar{-1};
    face best{0, 0, V{0, 0}, -std::numeric_limits<scalar>::infinity()};
    for (std::size_t i = 0; i < nv; ++i) {
      const V e = vsub2(view[(i + 1) % nv], view[i]);
      const scalar len = vnorm(e);
      if (len <= scalar{0}) continue;
      V o{e[1] / len, -e[0] / len};
      if (nv >= 3) {
        o = vscale(orient, o);
      } else if (vdot(o, d) < scalar{0}) {
        o = vscale(scalar{-1}, o);  // a segment faces both ways
      }
      const scalar al = vdot(o, d);
      if (al > best.align) best = face{i, (i + 1) % nv, o, al};
    }
    return best;
  };

  const face fa = best_face(va, n);
  const face fb = best_face(vb, vscale(scalar{-1}, n));
  constexpr scalar face_tol = scalar{1e-6};
  if (std::max(fa.align, fb.align) < scalar{1} - face_tol) return 0;

  const bool ref_a = !(fb.align > fa.align + scalar{1e-3});
  const auto& rview = ref_a ? va : vb;
  const auto& iview = ref_a ? vb : va;
  const face& rf = ref_a ? fa : fb;
  const face& inc = ref_a ? fb : fa;

  const V r0 = rview[rf.i0];
  const V r1 = rview[rf.i1];
  const V p0 = iview[inc.i0];
  const V p1 = iview[inc.i1];
  const V rt = vsub2(r1, r0);
  const scalar len = vnorm(rt);
  if (len <= scalar{0}) return 0;
  const V t = vscale(scalar{1} / len, rt);
  const scalar s0 = vdot(t, vsub2(p0, r0));
  const scalar s1 = vdot(t, vsub2(p1, r0));

  // the part of [p0, p1] whose tangential coordinate is in [0, len]
  scalar lmin, lmax;
  if (std::abs(s1 - s0) <= scalar{1e-12} * (scalar{1} + len)) {
    if (s0 < scalar{0} || s0 > len) return 0;
    lmin = scalar{0};
    lmax = scalar{1};
  } else {
    const scalar la = (scalar{0} - s0) / (s1 - s0);
    const scalar lb = (len - s0) / (s1 - s0);
    lmin = std::max(scalar{0}, std::min(la, lb));
    lmax = std::min(scalar{1}, std::max(la, lb));
    if (lmin > lmax) return 0;
  }
  const V dp = vsub2(p1, p0);
  std::array<V, 2> q{vadd(p0, vscale(lmin, dp)), vadd(p0, vscale(lmax, dp))};
  const std::size_t nq = (lmax - lmin) * vnorm(dp) > scalar{1e-12} * (scalar{1} + len) ? 2 : 1;

  for (std::size_t k = 0; k < nq; ++k) {
    const scalar gap = vdot(rf.normal, vsub2(q[k], r0));
    const V w = vsub2(q[k], vscale(gap, rf.normal));
    out[k].h = static_cast<Scalar>(gap);
    out[k].cp1 = ref_a ? w : q[k];
    out[k].cp2 = ref_a ? q[k] : w;
  }
  if (nq == 2) {
    const V tng{-n[1], n[0]};
    if (vdot(tng, out[1].cp1) < vdot(tng, out[0].cp1)) std::swap(out[0], out[1]);
  }
  return nq;
}

}  // namespace detail

template <typename Scalar>
inline std::vector<narrowphase_contact<Scalar>> narrowphase_batch(
    narrowphase_input<Scalar>&& in, std::size_t subdistance_budget = 64) {
  using cxd::status;
  using V = std::array<Scalar, 2>;

  const std::size_t npairs = in.ids_a.size();
  std::vector<narrowphase_contact<Scalar>> out(npairs);
  for (auto& c : out) {
    c.status = static_cast<int>(cxd::status::invalid_input);
    c.h = Scalar{0};
    c.normal = V{0, 0};
    c.cp1 = V{0, 0};
    c.cp2 = V{0, 0};
    c.npoints = 0;
  }

  if (in.entries.size() == 0) {
    return out;
  }

  // the pool takes ownership of the vertex buffer
  cxd::packed_pool<V> pool{std::move(in.vertices), std::move(in.entries)};

  std::vector<cxd::pair_record<Scalar, 2>> pairs(npairs);
  for (std::size_t k = 0; k < npairs; ++k) {
    V d0 = in.init_dirs[k];
    const Scalar n0 = std::hypot(d0[0], d0[1]);
    if (n0 <= Scalar{0})
      d0 = V{1, 0};
    else {
      d0 = V{d0[0] / n0, d0[1] / n0};
    }
    pairs[k] = cxd::pair_record<Scalar, 2>{in.ids_a[k], in.ids_b[k], {d0[0], d0[1]}};
  }

  // spec §5 tolerance/budget: tight tolerance, default iteration budget.
  // cxd's gjk_options carries the tolerance as a double; keep the float
  // case above epsilon so the GJK can terminate.
  const double tol = static_cast<double>(
      std::max(Scalar{1e-10}, Scalar{8} * std::numeric_limits<Scalar>::epsilon()));
  const cxd::gjk_options opts{.tolerance = tol, .max_iterations = 64};

  const auto results = cxd::gjk_batch<2>(pool, pairs, opts);

  // The enclosed pairs get their depth and normal from cxd's EPA. The 2D
  // EPA bootstraps a non-degenerate seed when the GJK hands one that is
  // degenerate -- the case a flush, axis-aligned overlap always is, since
  // the GJK's direction then collapses onto the line through the origin.
  std::vector<std::uint32_t> hit_pairs;
  for (std::size_t k = 0; k < npairs; ++k)
    if (results[k].st == status::enclosed) hit_pairs.push_back(static_cast<std::uint32_t>(k));
  const auto epa_answers = cxd::epa_batch<2>(pool, pairs, hit_pairs, results, opts);
  std::size_t hit_i = 0;

  // the separating-axis answer for a pair the GJK / EPA path left
  // unanswered; returns false (failure status kept) only for shapes with no
  // edge at all
  const auto sat_fallback = [&](std::size_t k, narrowphase_contact<Scalar>& c) {
    const auto sat = detail::sat_contact(pool, pairs[k].shape_id_A, pairs[k].shape_id_B);
    if (!sat) return false;
    const auto& [sep, nrm] = *sat;
    const auto [cp1, cp2] =
        detail::face_witnesses(pool, pairs[k].shape_id_A, pairs[k].shape_id_B, nrm);
    c.status =
        static_cast<int>(sep > Scalar{0} ? cxd::status::separated : cxd::status::enclosed);
    c.h = sep;
    c.normal = nrm;
    c.cp1 = cp1;
    c.cp2 = cp2;
    return true;
  };

  for (std::size_t k = 0; k < npairs; ++k) {
    auto& c = out[k];
    switch (results[k].st) {
      case status::enclosed: {
        c.status = static_cast<int>(cxd::status::enclosed);
        const auto& [depth, epa_normal] = epa_answers[hit_i++];
        // a touching pair can come back with a finite depth 0 and a zero
        // normal: that is no answer either
        if (std::isfinite(depth) && std::isfinite(epa_normal[0]) &&
            std::isfinite(epa_normal[1]) &&
            std::hypot(epa_normal[0], epa_normal[1]) > Scalar{0.5}) {
          const V nrm{epa_normal[0], epa_normal[1]};
          // the witnesses are the supporting faces along EPA's normal; h is
          // EPA's depth (the spec's mapping: enclosed = -depth)
          const auto [cp1, cp2] =
              detail::face_witnesses(pool, pairs[k].shape_id_A, pairs[k].shape_id_B, nrm);
          c.normal = nrm;
          c.cp1 = cp1;
          c.cp2 = cp2;
          c.h = -depth;
        } else if (!sat_fallback(k, c)) {
          // EPA refused (non-finite depth/normal) and no separating axis:
          // a failure status, so the stage keeps its previous cached row
          c.status = static_cast<int>(cxd::status::invalid_input);
          c.h = Scalar{0};
          c.normal = V{0, 0};
          c.cp1 = V{0, 0};
          c.cp2 = V{0, 0};
        }
        break;
      }
      case status::separated: {
        c.status = static_cast<int>(cxd::status::separated);
        const V dir0{pairs[k].init_dir[0], pairs[k].init_dir[1]};
        const auto sub = detail::subdistance_gjk(
            pool, pairs[k].shape_id_A, pairs[k].shape_id_B, dir0, subdistance_budget);
        if (sub) {
          const auto& [distance, foot, w, n] = *sub;
          const Scalar fn = detail::vnorm(foot);
          c.normal = fn > Scalar{0} ? detail::vscale(-Scalar{1} / fn, foot) : V{0, 0};
          // supporting-face witnesses along the contact normal; h from the
          // witness pair keeps h and the cps mutually consistent
          const auto [cp1, cp2] =
              detail::face_witnesses(pool, pairs[k].shape_id_A, pairs[k].shape_id_B, c.normal);
          c.cp1 = cp1;
          c.cp2 = cp2;
          const Scalar h_wit = detail::vdot(c.normal, detail::vsub2(cp2, cp1));
          c.h = h_wit > Scalar{0} ? h_wit : distance;
        } else if (!sat_fallback(k, c)) {
          // no separation certificate and no separating axis: a failure
          // status, so the stage keeps its previous cached row
          c.status = static_cast<int>(cxd::status::max_iterations);
        }
        break;
      }
      case status::max_iterations:
        // the GJK exhausted its budget: the separating axis still answers
        if (!sat_fallback(k, c)) c.status = static_cast<int>(cxd::status::max_iterations);
        break;
      case status::invalid_input:
      default:
        break;  // keep-previous contract (already zeroed)
    }
  }

  // the contact manifold of every answered pair
  for (std::size_t k = 0; k < npairs; ++k) {
    auto& c = out[k];
    if (c.status != static_cast<int>(cxd::status::separated) &&
        c.status != static_cast<int>(cxd::status::enclosed))
      continue;
    c.npoints = detail::clip_manifold(pool, pairs[k].shape_id_A, pairs[k].shape_id_B, c.normal,
                                      c.points);
    if (c.npoints == 0) {
      c.npoints = 1;
      c.points[0] = {c.h, c.cp1, c.cp2};
    }
  }
  return out;
}

/// gap given to a manifold point that does not exist (the manifold has one
/// point): far beyond any activation threshold, so its row stays inactive
template <typename Scalar>
constexpr Scalar inactive_gap() {
  return Scalar{1e30};
}

/// point i of the manifold of c, as a single-point contact (h, cp1, cp2);
/// a point the manifold does not have is the first point at inactive_gap
template <typename Scalar>
narrowphase_contact<Scalar> manifold_point(const narrowphase_contact<Scalar>& c,
                                           std::size_t i) {
  auto r = c;
  if (c.npoints == 0) return r;  // unanswered pair: the failure row as is
  const auto& p = c.points[i < c.npoints ? i : 0];
  r.h = i < c.npoints ? p.h : inactive_gap<Scalar>();
  r.cp1 = p.cp1;
  r.cp2 = p.cp2;
  return r;
}

}  // namespace siconos::collision

#include "siconos/collision/collision_head.hpp"
#include "siconos/collision/shape_of.hpp"
#include "siconos/model/lagrangian_ds.hpp"
#include "siconos/simul/topology.hpp"

namespace siconos::collision {

template <typename Topo>
struct polytopepolytope_r;  // polytopepolytope_r.hpp (Task 4)

/* ---- the narrowphase stage item (Task 5) ----
   Owns the world's polytope interactions: rebinds generic ds-ds relations
   to polytopepolytope_r items, packs world-vertex data once per call,
   runs the batch engine, writes statuses back and caches the step tag +
   a world checksum so ensure(step) re-runs only when the geometry moved. */
template <typename Topo>
struct cxd_narrowphase : item {
  using topology = Topo;
  using finteraction = typename topology::finteraction;
  using poly_r_t = polytopepolytope_r<Topo>;

  struct attributes {
    some::item_ref<topology> topology;
    some::item_ref<collision::diskdisk_r> dsds_r;  // generic binder to detect
    // the plan's poly_dsds_r ref (cxd_narrowphase -> polytopepolytope_r) is
    // DISABLED on purpose: the mutual item_refs with poly_r's narrowphase
    // ref make storage::make's all_items walker recurse forever (found by
    // bisection, see ledger); polytope relations are identified by dispatch
    // type in rebind/pack, so the back-ref is not needed
    some::indice cache_step;
    // FNV-1a over the pair set and the world vertices: a hash, not a scalar
    some::indice cache_checksum;
    // diagnostic: number of full-world scans (rebind + checksum) performed
    some::indice cache_scan_count;
    // instrumentation (additive, non-behavioural): seconds spent in each
    // phase of the LAST ensure() call. ensure() resets them, so they always
    // describe one call (a cache-hit step reports ~0). t_checksum covers
    // the ONE checksum() call of that ensure() (run() re-stamps the value
    // ensure() already computed instead of scanning the world again).
    some::scalar t_rebind;
    some::scalar t_checksum;
    some::scalar t_pack;
    some::scalar t_batch;
    some::scalar t_writeback;
    // instrumentation: state of the LAST COMPLETED batch (not reset by
    // ensure()): n_pairs = pairs packed into the batch, n_contacts = pairs
    // whose row was written (no failure status), n_enclosed = pairs the
    // batch reported as enclosed.
    some::indice n_pairs;
    some::indice n_contacts;
    some::indice n_enclosed;
    some::map<some::array<some::indice, some::indice_value<2>>,
              some::item_ref<poly_r_t>>
        rel_map;  // ordered ds-index pair -> polytope relation reuse map
  };

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) topology() {
      return storage::make_ref_handle(self()->data(), storage::attr<"topology">(*self()));
    }
    decltype(auto) rel_map() { return storage::attr<"rel_map">(*self()); }
    decltype(auto) dsds_r() { return storage::attr<"dsds_r">(*self()); }
    decltype(auto) cache_scan_count() { return storage::attr<"cache_scan_count">(*self()); }
    decltype(auto) cache_step() { return storage::attr<"cache_step">(*self()); }

    // instrumentation accessors (see the attribute comments)
    decltype(auto) t_rebind() { return storage::attr<"t_rebind">(*self()); }
    decltype(auto) t_checksum() { return storage::attr<"t_checksum">(*self()); }
    decltype(auto) t_pack() { return storage::attr<"t_pack">(*self()); }
    decltype(auto) t_batch() { return storage::attr<"t_batch">(*self()); }
    decltype(auto) t_writeback() { return storage::attr<"t_writeback">(*self()); }
    decltype(auto) n_pairs() { return storage::attr<"n_pairs">(*self()); }
    decltype(auto) n_contacts() { return storage::attr<"n_contacts">(*self()); }
    decltype(auto) n_enclosed() { return storage::attr<"n_enclosed">(*self()); }

    // seconds elapsed since t0 (instrumentation helper)
    static double seconds_since(std::chrono::steady_clock::time_point t0) {
      return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    }

    /* stream a ds' world vertices, CCW: world = c(q) + R(theta) local,
       c(q) = (q(0), q(1)), theta = q(2); the rotation-direction convention
       is pinned by an exact-value test (R(pi/2)*(1,0) == (0,1)). Streaming
       keeps the checksum allocation-free. */
    template <typename F>
    void for_each_world_vertex(auto& ds, auto step, F&& f) {
      using env_t = decltype(self()->env());
      using scalar = typename env_t::scalar;
      auto& data = self()->data();
      auto sh = collision::shape_of<shape::polytope>(ds);
      auto& q = storage::attr<"q">(ds, step);
      const scalar c = std::cos(q(2)), s = std::sin(q(2));
      const auto n = sh.nvertices();
      for (std::size_t i = 0; i < n; ++i) {
        const auto& v = sh.vertices()[i];
        f(std::array<scalar, 2>{q(0) + c * v[0] - s * v[1], q(1) + s * v[0] + c * v[1]});
      }
    }

    auto world_vertices(auto& ds, auto step) {
      using env_t = decltype(self()->env());
      using scalar = typename env_t::scalar;
      std::vector<std::array<scalar, 2>> out;
      for_each_world_vertex(ds, step, [&](const auto& v) { out.push_back(v); });
      return out;
    }

    /* 1. rebind: interactions between two ds whose shapes are polytopes and
       whose relation is still the generic diskdisk_r get a
       polytopepolytope_r item (created or reused from rel_map; reuse keeps
       init_dir); the relation ref then points at this stage item. */
    void rebind(auto step) {
      auto& data = self()->data();
      const auto t0 = std::chrono::steady_clock::now();
      for (auto inter : storage::handles<finteraction>(data, step)) {
        auto ds1 = storage::make_handle(data, storage::prop<"ds1", Topo>(inter));
        auto ds2 = storage::make_handle(data, storage::prop<"ds2", Topo>(inter));
        auto poly_shape = [](auto& ds) { return collision::has_shape<shape::polytope>(ds); };
        auto& relvar = storage::attr<"relation">(inter);
        // Case A: polytope x polytope (ds1 != ds2), generic binder diskdisk_r
        if (ds1.index().value() != ds2.index().value() && poly_shape(ds1) && poly_shape(ds2) &&
            std::holds_alternative<storage::index<collision::diskdisk_r, std::size_t>>(
                relvar)) {
          // generic binder still attached: give this pair a polytope relation
          const auto i1 = ds1.index().value(), i2 = ds2.index().value();
          const std::array<std::decay_t<decltype(i1)>, 2> key =
              i1 < i2 ? std::array{i1, i2} : std::array{i2, i1};
          auto& rmap = rel_map();
          auto rel = [&] {
            if (auto it = rmap.find(key); it != rmap.end()) {
              return storage::make_handle(data, it->second);
            }
            auto rel = storage::add<poly_r_t>(data);
            // fresh vector attrs are uninitialized storage (Eigen default);
            // give the first batch a well-defined direction: the engine maps
            // a zero init_dir to {1, 0} (spec §5). status/h define the
            // keep-previous baseline: no valid result yet.
            rel.status() = static_cast<std::size_t>(cxd::status::invalid_input);
            rel.h() = 0;
            rel.cache_pair_checksum() = 0;
            rel.normal() = {0, 0};
            rel.cp1() = {0, 0};
            rel.cp2() = {0, 0};
            rel.init_dir() = {0, 0};
            rmap[key] = rel;
            return rel;
          }();
          rel.narrowphase() = *self();
          inter.relation() = rel;
        }
        // Case B: polytope x fixed segment (self interaction, ds1==ds2),
        // generic binder diskfsegment_r. Rebind to polytopefsegment_r,
        // carrying the segment ref.
      }
      storage::attr<"t_rebind">(*self()) += seconds_since(t0);
    }

    /* FNV-1a over the PAIR SET and every involved ds' world vertices. The
       pair keys are load-bearing: a new interaction between two ds that are
       already in the checksum (both have another pair) changes no vertex and
       no ds index, but must refresh the cache. */
    auto checksum(auto step) {
      auto& data = self()->data();
      const auto t0 = std::chrono::steady_clock::now();
      std::uint64_t hsh = kFnvBasis;
      std::unordered_map<std::size_t, bool> seen;
      for (auto inter : storage::handles<finteraction>(data, step)) {
        auto& relvar = storage::attr<"relation">(inter);
        if (std::holds_alternative<storage::index<poly_r_t, std::size_t>>(relvar)) {
          auto ds1 = storage::make_handle(data, storage::prop<"ds1", Topo>(inter));
          auto ds2 = storage::make_handle(data, storage::prop<"ds2", Topo>(inter));
          const auto i1 = ds1.index().value();
          const auto i2 = ds2.index().value();
          fnv_mix(hsh, &i1, sizeof i1);
          fnv_mix(hsh, &i2, sizeof i2);
          for (auto* ds : {&ds1, &ds2}) {
            const auto idx = ds->index().value();
            if (seen[idx]) continue;
            seen[idx] = true;
            for_each_world_vertex(*ds, step, [&](const auto& v) {
              fnv_mix(hsh, &v[0], sizeof v[0]);
              fnv_mix(hsh, &v[1], sizeof v[1]);
            });
          }
        }
      }
      storage::attr<"t_checksum">(*self()) += seconds_since(t0);
      return hsh;
    }

    /* FNV-1a mixer over raw bytes; every hash below builds on it */
    static void fnv_mix(std::uint64_t& hsh, const void* p, std::size_t n) {
      const auto* bytes = static_cast<const unsigned char*>(p);
      for (std::size_t i = 0; i < n; ++i) {
        hsh ^= bytes[i];
        hsh *= 1099511628211ull;
      }
    }
    static constexpr std::uint64_t kFnvBasis = 1469598103934665603ull;

    /* FNV over ONE ds' world vertices. The per-pair signature combines the
       two bodies' hashes, and the batch writeback memoises this per-ds
       hash so a body attached to several pairs is hashed once per run */
    auto body_vertex_hash(auto& ds, auto step) {
      std::uint64_t hsh = kFnvBasis;
      for_each_world_vertex(ds, step, [&](const auto& v) {
        fnv_mix(hsh, &v[0], sizeof v[0]);
        fnv_mix(hsh, &v[1], sizeof v[1]);
      });
      return hsh;
    }

    /* combine the two ds indices and their per-body vertex hashes, order
       normalised by ds index so either handle order hashes the same */
    auto pair_signature(std::size_t i1, std::size_t i2, std::uint64_t h1, std::uint64_t h2) {
      std::uint64_t hsh = kFnvBasis;
      if (i1 <= i2) {
        fnv_mix(hsh, &i1, sizeof i1);
        fnv_mix(hsh, &i2, sizeof i2);
        fnv_mix(hsh, &h1, sizeof h1);
        fnv_mix(hsh, &h2, sizeof h2);
      } else {
        fnv_mix(hsh, &i2, sizeof i2);
        fnv_mix(hsh, &i1, sizeof i1);
        fnv_mix(hsh, &h2, sizeof h2);
        fnv_mix(hsh, &h1, sizeof h1);
      }
      return hsh;
    }

    /* per-pair signature over the pair's two bodies only, recomputed from
       the current world state (the relation gate re-derives it in
       O(vertices of the pair); the stage stamps it with memoised body
       hashes). Order normalised by ds index. */
    auto pair_checksum(auto& ds1, auto& ds2, auto step) {
      const auto i1 = ds1.index().value();
      const auto i2 = ds2.index().value();
      const auto h1 = body_vertex_hash(ds1, step);
      const auto h2 = body_vertex_hash(ds2, step);
      return pair_signature(i1, i2, h1, h2);
    }

    /* 2. pack the batch input: candidates are all interactions whose
       relation is a polytopepolytope_r (active or not: h continuity);
       each involved ds contributes its world vertices once */
    auto pack(auto step) {
      using env_t = decltype(self()->env());
      using scalar = typename env_t::scalar;
      using indice_t = typename env_t::indice;
      auto& data = self()->data();
      const auto t0 = std::chrono::steady_clock::now();
      narrowphase_input<scalar> in;
      std::unordered_map<std::size_t, std::uint32_t> slot;
      std::unordered_map<std::size_t, std::uint32_t> seg_slot;
      for (auto inter : storage::handles<finteraction>(data, step)) {
        auto& relvar = storage::attr<"relation">(inter);
        if (std::holds_alternative<storage::index<poly_r_t, std::size_t>>(relvar)) {
          auto rel = storage::make_handle(
              data, std::get<storage::index<poly_r_t, std::size_t>>(relvar));
          auto ds1 = storage::make_handle(data, storage::prop<"ds1", Topo>(inter));
          auto ds2 = storage::make_handle(data, storage::prop<"ds2", Topo>(inter));
          for (auto* ds : {&ds1, &ds2}) {
            if (slot.count(ds->index().value())) continue;
            const std::uint32_t offset = static_cast<std::uint32_t>(in.vertices.size());
            for (const auto& v : world_vertices(*ds, step)) {
              in.vertices.push_back(v);
            }
            in.entries.push_back(
                {offset, static_cast<std::uint32_t>(in.vertices.size() - offset)});
            // shape ids are ENTRY indices, not vertex offsets
            slot[ds->index().value()] = static_cast<std::uint32_t>(in.entries.size()) - 1;
          }
          in.ids_a.push_back(slot[ds1.index().value()]);
          in.ids_b.push_back(slot[ds2.index().value()]);
          in.init_dirs.push_back({rel.init_dir()[0], rel.init_dir()[1]});
        }
      }
      storage::attr<"n_pairs">(*self()) = static_cast<indice_t>(in.ids_a.size());
      storage::attr<"t_pack">(*self()) += seconds_since(t0);
      return in;
    }

    /* 3. batch + write-back: statuses/h/normal/cps go to the relations;
       init_dir becomes the contact normal when enclosed, else the unit
       vector from cp1 toward cp2; max_iterations/invalid_input keep the
       previous cached h and init_dir (keep-previous contract) */
    void run(auto step) {
      using env_t = decltype(self()->env());
      using indice_t = typename env_t::indice;
      auto& data = self()->data();
      auto in = pack(step);
      if (!in.ids_a.empty()) {
        const auto t_batch0 = std::chrono::steady_clock::now();
        auto out = narrowphase_batch(std::move(in));
        storage::attr<"t_batch">(*self()) += seconds_since(t_batch0);
        const auto t_wb0 = std::chrono::steady_clock::now();
        // each ds' vertex hash is memoised: a body attached to several
        // pairs is hashed once per run instead of once per pair
        std::unordered_map<std::size_t, std::uint64_t> ds_h;
        ds_h.reserve(in.ids_a.size());
        auto hash_ds = [&](auto& ds) -> std::uint64_t {
          const auto key = ds.index().value();
          const auto it = ds_h.find(key);
          if (it != ds_h.end()) return it->second;
          const auto h = body_vertex_hash(ds, step);
          ds_h.emplace(key, h);
          return h;
        };
        std::size_t k = 0;
        std::size_t n_written = 0;
        std::size_t n_enclosed = 0;
        for (auto inter : storage::handles<finteraction>(data, step)) {
          auto& relvar = storage::attr<"relation">(inter);
          if (std::holds_alternative<storage::index<poly_r_t, std::size_t>>(relvar)) {
            auto rel = storage::make_handle(
                data, std::get<storage::index<poly_r_t, std::size_t>>(relvar));
            auto ds1 = storage::make_handle(data, storage::prop<"ds1", Topo>(inter));
            auto ds2 = storage::make_handle(data, storage::prop<"ds2", Topo>(inter));
            // the pair is fresh for this configuration, even when the batch
            // kept the previous row: stamp it so relation reads skip the scan
            rel.cache_pair_checksum() = pair_signature(
                ds1.index().value(), ds2.index().value(), hash_ds(ds1), hash_ds(ds2));
            const auto& c = out[k++];
            if (c.status == static_cast<int>(cxd::status::max_iterations) ||
                c.status == static_cast<int>(cxd::status::invalid_input)) {
              continue;  // keep-previous: the whole cached row stays
            }
            rel.status() = c.status;
            rel.h() = c.h;
            rel.normal() = {c.normal[0], c.normal[1]};
            rel.cp1() = {c.cp1[0], c.cp1[1]};
            rel.cp2() = {c.cp2[0], c.cp2[1]};
            ++n_written;
            if (c.status == static_cast<int>(cxd::status::enclosed)) {
              ++n_enclosed;
              rel.init_dir() = {rel.normal()[0], rel.normal()[1]};
            } else if (c.h > 0) {
              // separated: the warm start points from body 1 to body 2, i.e.
              // -(cp1 - cp2)/h (cp2 - cp1 = h * normal for a separated pair)
              using scalar = typename std::decay_t<decltype(self()->env())>::scalar;
              const scalar dx = rel.cp2()[0] - rel.cp1()[0];
              const scalar dy = rel.cp2()[1] - rel.cp1()[1];
              const scalar d = std::hypot(dx, dy);
              if (d > scalar{0}) rel.init_dir() = {dx / d, dy / d};
            }
          }
        }
        storage::attr<"n_contacts">(*self()) = static_cast<indice_t>(n_written);
        storage::attr<"n_enclosed">(*self()) = static_cast<indice_t>(n_enclosed);
        storage::attr<"t_writeback">(*self()) += seconds_since(t_wb0);
      }
      storage::attr<"cache_step">(*self()) = step;
      // the world checksum is NOT recomputed here: ensure() computed it
      // just before calling run() and nothing since changed the ds data
      // the hash reads, so ensure() stamps the value it already has
    }

    /* ensure(step): rebind is idempotent and cheap, so it runs every call;
       the batch re-runs only on a new step tag or a changed world checksum */
    void ensure(auto step) {
      storage::attr<"cache_scan_count">(*self()) += 1;
      // the phase timers describe THIS call: reset them so they never
      // accumulate across steps (a cache-hit step reports ~0)
      storage::attr<"t_rebind">(*self()) = 0;
      storage::attr<"t_checksum">(*self()) = 0;
      storage::attr<"t_pack">(*self()) = 0;
      storage::attr<"t_batch">(*self()) = 0;
      storage::attr<"t_writeback">(*self()) = 0;
      rebind(step);
      const auto ck = checksum(step);
      if (storage::attr<"cache_step">(*self()) == step &&
          storage::attr<"cache_checksum">(*self()) == ck) {
        return;
      }
      run(step);
      // run() changed relations only: the ds data checksum() read is the
      // same, so the stamp is the value already computed above
      storage::attr<"cache_checksum">(*self()) = ck;
    }
  };
};

}  // namespace siconos::collision
