#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace siconos::collision::shape {

struct polytope : item {
  struct attributes {
    // capacity-16 local vertices (v1); variable count via nvertices
    some::array<some::vector<some::scalar, some::indice_value<2>>, some::indice_value<16>>
        vertices;
    some::indice nvertices;
    // bounding radius cached for the generic ds-ds binder: collision::diskdisk_r
    // reads shape.radius() (an lvalue) on every relation compute, so a
    // polytope attached to a ds must expose one; radius() refreshes it from
    // rmax() (the true circumradius of the local polygon)
    some::scalar radius;
  };

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) vertices() { return attr<"vertices">(*self()); }

    decltype(auto) nvertices() { return attr<"nvertices">(*self()); }

    // max circumradius (used by broadphase radius sizing): max distance
    // from the local origin (= ds position q) to a vertex
    decltype(auto) rmax() {
      using env_t = decltype(self()->env());
      using scalar = typename env_t::scalar;
      scalar m = 0;
      for (auto i = 0u; i < nvertices(); ++i) {
        const auto& v = vertices()[i];
        m = std::max(m, std::hypot(v[0], v[1]));
      }
      return m;
    }

    // radius lvalue for the generic ds-ds binder (collision::diskdisk_r binds
    // `auto& r = shape.radius()`): refreshed from rmax() on every access
    decltype(auto) radius() {
      attr<"radius">(*self()) = rmax();
      return attr<"radius">(*self());
    }

    /// set the local vertices, one per row (x, y; a third column is
    /// ignored), relative to the body's centre of mass. Checks 3 to 16
    /// vertices, a convex polygon of non-zero area without repeated vertices,
    /// and stores it counter-clockwise (the orientation the contact code
    /// expects); a clockwise list is reversed.
    template <typename M>
    void set_vertices(M& v) {
      using env_t = decltype(self()->env());
      using scalar = typename env_t::scalar;
      const std::size_t n = algebra::rows(v);
      if (n < 3 || n > std::size(vertices())) {
        throw std::runtime_error("shape::polytope::set_vertices: 3 to 16 vertices expected");
      }
      scalar area2 = 0;
      for (std::size_t i = 0; i < n; ++i) {
        const std::size_t k = (i + 1) % n;
        area2 += v(i, 0) * v(k, 1) - v(k, 0) * v(i, 1);
      }
      if (!(std::abs(area2) > scalar{0})) {
        throw std::runtime_error("shape::polytope::set_vertices: polygon of zero area");
      }
      for (std::size_t i = 0; i < n; ++i) {
        const std::size_t src = area2 > scalar{0} ? i : n - 1 - i;
        vertices()[i] = {static_cast<scalar>(v(src, 0)), static_cast<scalar>(v(src, 1))};
      }
      nvertices() = n;
      check_convex();
    }

    /// convex, counter-clockwise, no repeated consecutive vertex: every
    /// turn is to the left (collinear vertices are accepted)
    void check_convex() {
      using env_t = decltype(self()->env());
      using scalar = typename env_t::scalar;
      const std::size_t n = nvertices();
      if (n > std::size(vertices())) {
        throw std::runtime_error("shape::polytope: nvertices exceeds the 16-vertex capacity");
      }
      for (std::size_t i = 0; i < n; ++i) {
        const auto& a = vertices()[i];
        const auto& b = vertices()[(i + 1) % n];
        const auto& c = vertices()[(i + 2) % n];
        if (a[0] == b[0] && a[1] == b[1]) {
          throw std::runtime_error("shape::polytope: repeated vertex");
        }
        const scalar cr = (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0]);
        if (cr < scalar{0}) {
          throw std::runtime_error(
              "shape::polytope: non-convex polygon (a right turn); give its convex hull");
        }
      }
    }

    auto methods() {
      using env_t = decltype(self()->env());
      using scalar = typename env_t::scalar;
      using unbounded_matrix = typename env_t::template unbounded_matrix<scalar>;
      return collect(
          method("set_vertices", &interface<Handle>::set_vertices<unbounded_matrix>));
    }
  };
};

}  // namespace siconos::collision::shape
