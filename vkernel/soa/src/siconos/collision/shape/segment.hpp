#pragma once

#include "siconos/algebra/eigen.hpp"
#include "siconos/collision/collision.hpp"
#include "siconos/collision/collision_head.hpp"

namespace siconos::collision::shape {

struct segment : item<> {
  using attributes = gather<
      attribute<"points", some::vector<some::scalar, some::indice_value<6>>>,
      attribute<"normal", some::vector<some::scalar, some::indice_value<3>>>,
      attribute<"dp2p1", some::vector<some::scalar, some::indice_value<3>>>,
      attribute<"maxpoints", some::scalar>,
      attribute<"length_sq", some::scalar>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) p1()
    {
      using scalar = typename decltype(self()->env())::scalar;

      return algebra::matrix_view<algebra::vector<scalar, 3>>(
          attr<"points">(*self()).data());
    };
    decltype(auto) p2()
    {
      using scalar = typename decltype(self()->env())::scalar;

      return algebra::matrix_view<algebra::vector<scalar, 3>>(
          attr<"points">(*self()).data() + 3);
    };
    decltype(auto) x1() { return p1()[0]; };
    decltype(auto) y1() { return p1()[1]; };
    decltype(auto) x2() { return p2()[0]; };
    decltype(auto) y2() { return p2()[1]; };
    decltype(auto) dp2p1() { return attr<"dp2p1">(*self()); };
    decltype(auto) maxpoints() { return attr<"maxpoints">(*self()); };
    decltype(auto) length_sq() { return attr<"length_sq">(*self()); };
    decltype(auto) normal() { return attr<"normal">(*self()); };

    void compute_dp2p1() { dp2p1() = p2() - p1(); };
    void compute_length_sq()
    {
      const auto& v = dp2p1();
      length_sq() = algebra::dot(v, v);
    };

    void compute_normal()
    {
      const auto& v = dp2p1();
      normal()[0] = -v[1];
      normal()[1] = v[0];
    }

    void initialize()
    {
      compute_dp2p1();
      compute_length_sq();
      compute_normal();
    }

    decltype(auto) distance(match::vector auto& q)
    {
      /* dof 3 -> 2D + 1 (CompactNSearch) */
      auto qp = q;
      qp[2] = 0.;

      const auto t =
          fmax(0, fmin(1, algebra::dot(qp - p1(), dp2p1()) / length_sq()));
      const auto p = p1() + t * dp2p1();
      return collision::distance(qp, p);
    }

    decltype(auto) points_coords()
    {
      const auto p = p1();
      const auto pstep = 1. / maxpoints();
      const auto dir = dp2p1();
      return view::iota(0, maxpoints()) |
             view::transform([=](auto i) { return p + i * pstep * dir; });
    }

    auto methods()
    {
      return collect(method("initialize", &interface<Handle>::initialize));
    }
  };
};
}  // namespace siconos::collision::shape
