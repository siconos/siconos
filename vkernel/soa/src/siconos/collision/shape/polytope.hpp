#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include <print>
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

    // convexity check: all consecutive edge cross products same sign.
    // tolerated (no-op) when nvertices() is still 0 at init time.
    void __init__() {
      const auto n = nvertices();
      if (n < 3) {
        return;
      }
      double sign = 0;
      for (auto i = 0u; i < n; ++i) {
        const auto& a = vertices()[i];
        const auto& b = vertices()[(i + 1) % n];
        const auto& c = vertices()[(i + 2) % n];
        const auto cr = (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0]);
        if (cr == 0) continue;
        if (sign == 0) sign = cr;
        if (cr * sign < 0) {
          throw std::runtime_error(
              "shape::polytope: non-convex polygon (inconsistent edge cross "
              "product sign); convex hull precomputation needed");
        }
      }
    }
  };
};

}  // namespace siconos::collision::shape
