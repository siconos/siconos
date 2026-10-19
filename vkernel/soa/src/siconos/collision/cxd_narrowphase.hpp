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
     enclosed (origin inside): exact minimum translation depth over the
     edge normals of both polygons (see enclosed_depth for why this is
     siconos-side rather than cxd's 2D epa, which rejects the degenerate
     seeds that flush contacts deterministically produce); witnesses from
     the same face construction along the deepest-face normal; the
     ambiguous-tie case (coincident/symmetric bodies) falls back to h = 0.
   - cp1/cp2 for both cases are the supporting faces' points (uniquely
     determined, deterministic); cxd's support tuples stay inside the
     GJK internals.

   Separation sub-distance: computed siconos-side with a small 2D distance
   GJK over the pool, using cxd's own support primitives (shape_view +
   cxd::diff_support). cxd's batch deliberately leaves the seed EMPTY on
   separated statuses ("misses keep every slot untouched"), so the
   seed-projection sub-distance of the spec has to travel with the engine.

   max_iterations / invalid_input: h = 0 placeholder + status flag; the
   CALLER (the relation item) keeps its previous cached h. No output is
   ever non-finite. */

#include <array>
#include <cmath>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <tuple>
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

struct narrowphase_contact {
  int status;                    // cxd::status value
  double h;                      // siconos gap: +distance / -depth
  std::array<double, 2> normal;  // contact normal, body 1 -> body 2
  std::array<double, 2> cp1, cp2;
};

struct narrowphase_input {
  std::vector<std::array<double, 2>>
      vertices;                           // shape vertices concatenated, world frame, double
  std::vector<cxd::shape_entry> entries;  // one per shape id: offset/count
  std::vector<std::uint32_t> ids_a,
      ids_b;                                     // shape ids per pair (id = entries index)
  std::vector<std::array<double, 2>> init_dirs;  // zero -> {1,0}
};

namespace detail {

using V = std::array<double, 2>;
using support_t = std::tuple<std::size_t, std::size_t, V>;

constexpr double vdot(V a, V b) { return a[0] * b[0] + a[1] * b[1]; }
inline double vnorm(V a) { return std::hypot(a[0], a[1]); }
inline V vadd(V a, V b) { return {a[0] + b[0], a[1] + b[1]}; }
inline V vscale(double s, V a) { return {s * a[0], s * a[1]}; }

inline const V& support_pos(const support_t& s) { return std::get<2>(s); }

inline V vsub2(V a, V b) { return {a[0] - b[0], a[1] - b[1]}; }

// closest point of the origin to the segment [p, q]: foot and t in [0,1]
// (t == 0 / 1 picks the vertex feature)
inline std::pair<V, double> foot_on_segment(V p, V q) {
  const V e = vsub2(q, p);
  const double ee = vdot(e, e);
  if (ee <= 0) return {p, 0.0};
  double t = -vdot(p, e) / ee;
  t = std::max(0.0, std::min(1.0, t));
  return {vadd(p, vscale(t, e)), t};
}

/* One reduction step: closest point of the origin to the convex hull of
   the (up to 3) supports; returns the foot, per-support weights (same
   count as the given simplex range), and a flag whether the foot is on a
   VERTEX (so the simplex reduces to that vertex) or on an EDGE (reduces
   to that edge's two supports). */
inline auto closest_to_origin(const std::array<support_t, 3>& w, std::size_t n)
    -> std::tuple<V, std::array<double, 3>, bool /*vertex feature*/,
                  std::size_t /*feature index*/> {
  std::array<double, 3> alpha{};
  V best_foot{0., 0.};
  double best_norm = std::numeric_limits<double>::infinity();
  bool best_is_vertex = false;
  std::size_t best_index = 0;

  // vertex features
  for (std::size_t i = 0; i < n; ++i) {
    const auto& p = support_pos(w[i]);
    if (const double d = vnorm(p); d < best_norm) {
      std::array<double, 3> a{};
      a[i] = 1.0;
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
    if (const double d = vnorm(foot); d < best_norm) {
      std::array<double, 3> a{};
      a[i] = 1.0 - t;
      a[(i + 1) % n] = t;
      best_norm = d;
      best_foot = foot;
      best_is_vertex = (t == 0.0 || t == 1.0);  // degenerate to a vertex
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
inline auto subdistance_gjk(const Pool& pool, std::uint32_t id_a, std::uint32_t id_b, V dir0,
                            std::size_t budget)
    -> std::optional<std::tuple<double /*distance*/, V /*foot*/, std::array<support_t, 3>,
                                std::size_t /*n*/>> {
  const auto ds = [&](V d) -> support_t {
    return cxd::diff_support(cxd::shape_view<Pool>{pool, id_a},
                             cxd::shape_view<Pool>{pool, id_b}, d);
  };

  const double tol = 1e-10;

  std::array<support_t, 3> w{};
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
      const support_t wa = w[fidx];
      const support_t wb = w[(fidx + 1) % n];
      w[0] = wa;
      w[1] = wb;
      n = 2;
    }

    const support_t cand = ds(dir);
    // no support strictly beyond the current foot: exact sub-distance
    // (dir = -foot, so the support bound and the simplex bound coincide)
    if (vdot(support_pos(cand), dir) <= -vdot(foot, foot) + tol) {
      double distance = vnorm(foot);
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

// barycentric cps from the weights + witness tuples
template <typename Pool>
inline auto witness_points(const Pool& pool, std::uint32_t id_a, std::uint32_t id_b,
                           const std::array<support_t, 3>& w,
                           const std::array<double, 3>& alpha, std::size_t n)
    -> std::pair<V, V> {
  const cxd::shape_view<Pool> va{pool, id_a};
  const cxd::shape_view<Pool> vb{pool, id_b};
  V cp1{0., 0.}, cp2{0., 0.};
  for (std::size_t i = 0; i < n; ++i) {
    const auto& [ia, ib, pos] = w[i];
    cp1 = vadd(cp1, vscale(alpha[i], va[ia]));
    cp2 = vadd(cp2, vscale(alpha[i], vb[ib]));
  }
  return {cp1, cp2};
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
inline auto face_witnesses(const Pool& pool, std::uint32_t id_a, std::uint32_t id_b, V n)
    -> std::pair<V, V> {
  const cxd::shape_view<Pool> va{pool, id_a};
  const cxd::shape_view<Pool> vb{pool, id_b};

  const auto nverts = [](const auto& view) {
    return static_cast<std::size_t>(std::distance(view.begin(), view.end()));
  };
  const double amax = [&] {
    double m = -std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < nverts(va); ++i) m = std::max(m, vdot(n, va[i]));
    return m;
  }();
  const double bmin = [&] {
    double m = std::numeric_limits<double>::infinity();
    for (std::size_t j = 0; j < nverts(vb); ++j) m = std::min(m, vdot(n, vb[j]));
    return m;
  }();
  const double eps = 1e-9 * (1.0 + std::abs(amax) + std::abs(bmin));

  std::vector<std::size_t> tie_a, tie_b;
  for (std::size_t i = 0; i < nverts(va); ++i)
    if (vdot(n, va[i]) >= amax - eps) tie_a.push_back(i);
  for (std::size_t j = 0; j < nverts(vb); ++j)
    if (vdot(n, vb[j]) <= bmin + eps) tie_b.push_back(j);

  if (tie_a.size() == 1 && tie_b.size() == 1) return {va[tie_a[0]], vb[tie_b[0]]};

  if (tie_a.size() == 1) {  // vertex on A, face on B: perpendicular foot
    const V p = va[tie_a[0]];
    const auto [foot, t] = foot_on_segment(vb[tie_b.front()], vb[tie_b.back()]);
    return {p, foot};
  }
  if (tie_b.size() == 1) {  // face on A, vertex on B
    const V p = vb[tie_b[0]];
    const auto [foot, t] = foot_on_segment(va[tie_a.front()], va[tie_a.back()]);
    return {foot, p};
  }

  // both faces (parallel edges): midpoint of the overlap window along the
  // common edge direction
  const V e0 = vsub2(vb[tie_b.back()], vb[tie_b.front()]);
  const double el = vnorm(e0);
  const V e = el > 0 ? vscale(1.0 / el, e0) : V{1.0, 0.0};
  const double pa0 = vdot(e, va[tie_a.front()]);
  const double pa1 = vdot(e, va[tie_a.back()]);
  const double pb0 = vdot(e, vb[tie_b.front()]);
  const double pb1 = vdot(e, vb[tie_b.back()]);
  const double lo = std::max(std::min(pa0, pa1), std::min(pb0, pb1));
  const double hi = std::min(std::max(pa0, pa1), std::max(pb0, pb1));
  const double tm = 0.5 * (lo + hi);
  const auto lerp_on = [&](const auto& view, std::size_t i0, std::size_t i1, double p0,
                           double p1) -> V {
    const double denom = p1 - p0;
    if (std::abs(denom) < 1e-30) return view[i0];
    const double t = std::max(0.0, std::min(1.0, (tm - p0) / denom));
    return vadd(view[i0], vscale(t, vsub2(view[i1], view[i0])));
  };
  return {lerp_on(va, tie_a.front(), tie_a.back(), pa0, pa1),
          lerp_on(vb, tie_b.front(), tie_b.back(), pb0, pb1)};
}

/* Enclosed-pair measurement, exact for convex polygons.

   The Minkowski difference's faces all have the normal of some edge of
   A or of B, so the minimum translation depth is the minimum overlap
        depth(n) = max_A (n.a) - min_B (n.b)
   over the unit edge normals of both polygons. The argmin direction is
   the exterior normal of the deepest diff face: unique -> face_witnesses
   along it yields the deep contact pair, h = n.(cp2-cp1) = -depth.
   A TIED argmin means the deepest face is ambiguous (concentric /
   symmetric configurations: for coincident squares four directions tie
   at depth 2): no canonical contact exists, the caller falls back.

   cxd's 2D epa rejects degenerate seeds outright (its repair paths are
   3D only), and flush contacts -- the steady state of the warm-start
   writeback, init_dir = normal -- deterministically produce those
   degenerate seeds, so the depth is computed here instead of epa_batch. */
template <typename Pool>
inline auto enclosed_depth(const Pool& pool, std::uint32_t id_a, std::uint32_t id_b)
    -> std::optional<std::tuple<V, V, V>> {  // normal, cp1, cp2
  const cxd::shape_view<Pool> va{pool, id_a};
  const cxd::shape_view<Pool> vb{pool, id_b};
  const auto nverts = [](const auto& view) {
    return static_cast<std::size_t>(std::distance(view.begin(), view.end()));
  };
  const std::size_t na = nverts(va), nb = nverts(vb);
  if (na < 3 || nb < 3 || (na + nb) > 32) return std::nullopt;

  // overlap along a unit direction n
  const auto overlap = [&va, &vb, na, nb](const V& n) {
    double amax = -std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < na; ++i) amax = std::max(amax, vdot(n, va[i]));
    double bmin = std::numeric_limits<double>::infinity();
    for (std::size_t j = 0; j < nb; ++j) bmin = std::min(bmin, vdot(n, vb[j]));
    return amax - bmin;
  };

  // candidate unit normals: perpendiculars of every edge, both signs.
  // Parallel edges (dup vertices included) give identical candidates, so
  // dedupe before scoring -- ties must be counted over DISTINCT directions.
  std::vector<V> cands;
  const auto cand_tol = 1e-12;
  const auto push_cand = [&](const V& n) {
    for (const auto& m : cands)
      if (n[0] * m[0] + n[1] * m[1] > 1.0 - cand_tol) return;
    cands.push_back(n);
  };
  for (const cxd::shape_view<Pool>* view : {&va, &vb}) {
    for (std::size_t i = 0; i < nverts(*view); ++i) {
      const V p = (*view)[i];
      const V q = (*view)[(i + 1) % nverts(*view)];
      const V e = vsub2(q, p);
      const double el = vnorm(e);
      if (el <= 1e-30) continue;
      const V n = V{-e[1] / el, e[0] / el};
      push_cand(n);
      push_cand(V{-n[0], -n[1]});
    }
  }
  if (cands.empty()) return std::nullopt;

  double best = std::numeric_limits<double>::infinity();
  std::size_t ties = 0;
  V best_n{};
  for (const auto& n : cands) {
    const double d = overlap(n);
    const double tol = 1e-9 * (1.0 + std::abs(best) + std::abs(d));
    if (d < best - tol) {
      best = d;
      best_n = n;
      ties = 1;
    } else if (d <= best + tol) {
      if (d < best) best = d;  // tighten within the tie band
      best_n = ties == 0 ? n : best_n;
      ++ties;
    }
  }

  if (!std::isfinite(best) || best <= 0.0 || ties != 1) return std::nullopt;
  const auto [cp1, cp2] = face_witnesses(pool, id_a, id_b, best_n);
  return std::tuple{best_n, cp1, cp2};
}

}  // namespace detail

inline std::vector<narrowphase_contact> narrowphase_batch(narrowphase_input&& in) {
  using cxd::status;
  using V = std::array<double, 2>;

  const std::size_t npairs = in.ids_a.size();
  std::vector<narrowphase_contact> out(npairs);
  for (auto& c : out) {
    c.status = static_cast<int>(status::invalid_input);
    c.h = 0.0;
    c.normal = {0.0, 0.0};
    c.cp1 = {0.0, 0.0};
    c.cp2 = {0.0, 0.0};
  }

  if (in.entries.size() == 0) return out;

  // the pool takes ownership of the vertex buffer
  cxd::packed_pool<V> pool{std::move(in.vertices), std::move(in.entries)};

  std::vector<cxd::pair_record<double, 2>> pairs(npairs);
  for (std::size_t k = 0; k < npairs; ++k) {
    V d0 = in.init_dirs[k];
    const double n0 = std::hypot(d0[0], d0[1]);
    if (n0 <= 0)
      d0 = {1.0, 0.0};
    else {
      d0 = {d0[0] / n0, d0[1] / n0};
    }
    pairs[k] = cxd::pair_record<double, 2>{in.ids_a[k], in.ids_b[k], {d0[0], d0[1]}};
  }

  // spec §5 tolerance/budget: tight tolerance, default iteration budget
  const cxd::gjk_options opts{.tolerance = 1e-6, .max_iterations = 64};

  const auto results = cxd::gjk_batch<2>(pool, pairs, opts);

  for (std::size_t k = 0; k < npairs; ++k) {
    auto& c = out[k];
    switch (results[k].st) {
      case status::enclosed: {
        c.status = static_cast<int>(status::enclosed);
        // deepest-face measurement, siconos-side (see enclosed_depth:
        // cxd's 2D epa rejects the degenerate seeds that flush contacts
        // deterministically produce)
        const auto m = detail::enclosed_depth(pool, pairs[k].shape_id_A, pairs[k].shape_id_B);
        if (m) {
          const auto& [nrm, cp1, cp2] = *m;
          c.normal = nrm;
          c.cp1 = cp1;
          c.cp2 = cp2;
          c.h = detail::vdot(nrm, detail::vsub2(cp2, cp1));
        } else {  // ambiguous deepest face (e.g. coincident bodies):
                  // keep-previous contract, finite outputs
          c.h = 0.0;
          c.normal = {0.0, 0.0};
          c.cp1 = {0.0, 0.0};
          c.cp2 = {0.0, 0.0};
        }
        break;
      }
      case status::separated: {
        c.status = static_cast<int>(status::separated);
        const V dir0{pairs[k].init_dir[0], pairs[k].init_dir[1]};
        const auto sub =
            detail::subdistance_gjk(pool, pairs[k].shape_id_A, pairs[k].shape_id_B, dir0, 64);
        if (sub) {
          const auto& [distance, foot, w, n] = *sub;
          const double fn = detail::vnorm(foot);
          c.normal = fn > 0 ? detail::vscale(-1.0 / fn, foot) : V{0., 0.};
          // supporting-face witnesses along the contact normal; h from the
          // witness pair keeps h and the cps mutually consistent
          const auto [cp1, cp2] =
              detail::face_witnesses(pool, pairs[k].shape_id_A, pairs[k].shape_id_B, c.normal);
          c.cp1 = cp1;
          c.cp2 = cp2;
          const double h_wit = detail::vdot(c.normal, detail::vsub2(cp2, cp1));
          c.h = h_wit > 0 ? h_wit : distance;
        }
        // else: keep-previous contract (status separated, h untouched)
        break;
      }
      case status::max_iterations:
      case status::invalid_input:
      default:
        break;  // keep-previous contract (already zeroed)
    }
  }
  return out;
}

}  // namespace siconos::collision
