#pragma once

#include <print>

#include "siconos/algebra/eigen.hpp"
#include "siconos/collision/collision.hpp"
#include "siconos/collision/collision_head.hpp"

namespace siconos::collision::shape {

struct segment : item {
  using indice_t = std::size_t;
  template <typename T>
  using vector_t = algebra::vector<T, 6>;

  struct attributes {
    // fixed vector of size 1 => same interface for an unbounded vector
    // in the case of chained segment
    some::array<some::vector<some::scalar, some::indice_value<3>>,
                 some::indice_value<2>>
        nodes;
    some::array<some::vector<some::scalar, some::indice_value<3>>,
                 some::indice_value<1>>
        dp2p1;
    some::scalar maxpoints;
    some::array<some::scalar, some::indice_value<1>> length_sq;
  };

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) nodes() { return attr<"nodes">(*self()); }

    decltype(auto) p1(indice_t seg_index = 0)
    {
      return nodes()[seg_index * 2];
    };
    decltype(auto) p2(indice_t seg_index = 0)
    {
      return nodes()[seg_index * 2 + 1];
    };
    decltype(auto) x1(indice_t seg_index = 0) { return p1(seg_index)[0]; };
    decltype(auto) y1(indice_t seg_index = 0) { return p1(seg_index)[1]; };
    decltype(auto) x2(indice_t seg_index = 0) { return p2(seg_index)[0]; };
    decltype(auto) y2(indice_t seg_index = 0) { return p2(seg_index)[1]; };
    decltype(auto) dp2p1(indice_t seg_index = 0)
    {
      return attr<"dp2p1">(*self())[seg_index];
    };
    decltype(auto) maxpoints() { return attr<"maxpoints">(*self()); };
    decltype(auto) length_sq(indice_t seg_index = 0)
    {
      return attr<"length_sq">(*self())[seg_index];
    };

    void compute_dp2p1(indice_t seg_index = 0)
    {
      dp2p1(seg_index) = p2(seg_index) - p1(seg_index);
    };
    void compute_length_sq(indice_t seg_index = 0)
    {
      const auto& v = dp2p1(seg_index);
      length_sq() = algebra::dot(v, v);
    };

    void initialize(indice_t seg_index = 0)
    {
      compute_dp2p1(seg_index);
      compute_length_sq(seg_index);
    }

    decltype(auto) distance(match::vector auto& q, indice_t seg_index = 0)
    {
      /* dof 3 -> 2D + 1 (CompactNSearch) */
      auto qp = q;
      qp[2] = 0.;

      const auto t =
          fmax(0, fmin(1, algebra::dot(qp - p1(), dp2p1()) / length_sq()));
      const auto p = p1() + t * dp2p1();
      return collision::distance(qp, p);
    }

    decltype(auto) point_coord(indice_t point_index, indice_t seg_index = 0)
    {
      const auto p = p1(seg_index);
      const auto dir = dp2p1(seg_index);
      const auto pstep = 1. / maxpoints();

      // explicit return type to avoid dangling references
      decltype(p) return_value = p + point_index * pstep * dir;
      return return_value;
    }

    decltype(auto) points_coords(indice_t seg_index = 0)
    {
      const auto p = p1(seg_index);
      const auto dir = dp2p1(seg_index);
      const auto pstep = 1. / maxpoints();

      return view::iota(0, maxpoints()) |
             // return expression is ok as p and dir are copied into lambda
             // closure.
             view::transform([=](auto i) { return p + i * pstep * dir; });
    }

    void set_p1_p2(auto nodes_array, indice_t seg_index = 0)
    {
      p1(seg_index) = {nodes_array[0], nodes_array[1], nodes_array[2]};
      p2(seg_index) = {nodes_array[3], nodes_array[4], nodes_array[5]};
    }

    void set_maxpoints(indice_t mp) { maxpoints() = mp; }

    template <typename Scalar>
    auto insert(Scalar x, Scalar y, Scalar z = 0)
    {
      using env_t = decltype(self()->env());
      using indice = typename env_t::indice;

      using nodes_store_t = std::decay_t<decltype(nodes())>;

      if constexpr (match::push_back<nodes_store_t>) {
        nodes().push_back({x, y, 0.}); /* 2D */

        /* initialize must be called on even sizes */
        indice number_of_nodes = std::size(nodes());
        if (number_of_nodes % 2 == 0) self()->initialize(number_of_nodes);
      }
      else {
        // throw an exception ?
      };
    }

    auto methods()
    {
      using env_t = decltype(self()->env());
      using scalar = typename env_t::scalar;

      return collect(
          method("initialize", &interface<Handle>::initialize),
          method("set_maxpoints", &interface<Handle>::set_maxpoints),
          method("set_p1_p2",
                 &interface<Handle>::set_p1_p2<vector_t<scalar>>),
          method("insert", &interface<Handle>::insert<scalar>));
    }
  };
};


}  // namespace siconos::collision::shape
