#pragma once

#include "siconos/algebra/eigen.hpp"
#include "siconos/collision/collision.hpp"
#include "siconos/collision/collision_head.hpp"
#include "siconos/collision/shape/segment.hpp"

namespace siconos::collision::shape {

struct chained_segment : segment_base {
  // definition of segment attributes is over written
  // segment is able to manage the unbounded vector case

  using without_attributes_bindings = void;
  using without_attached_storages_bindings = void;

  struct attributes {
    some::unbounded_vector<some::vector<some::scalar, some::indice_value<3>>>
        points;
    some::unbounded_vector<some::vector<some::scalar, some::indice_value<3>>>
        dp2p1;
    some::scalar maxpoints;
    some::unbounded_vector<some::scalar> length_sq;
  };
};
}  // namespace siconos::collision::shape
