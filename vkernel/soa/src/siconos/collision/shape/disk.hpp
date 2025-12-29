#pragma once

namespace siconos::collision::shape {
struct disk : item {

  struct attributes {
    some::scalar radius;
    some::indice maxpoints;
  };

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    void __init__() { maxpoints() = 0; };

    decltype(auto) radius() { return attr<"radius">(*self()); };

    decltype(auto) maxpoints() { return attr<"maxpoints">(*self()); };

    decltype(auto) point_coord(auto point_index)
    {
      using env_t = decltype(self()->env());
      using scalar = typename env_t::scalar;
      using vector_t = typename env_t::template vector<scalar, 3>;

      if (maxpoints() < 3) {
        // not enough points, returns the center of the disk
        vector_t coord = {0., 0., 0.};
        return coord;
      }
      else {
        // returns a point on the circle.
        scalar t = 2 * std::numbers::pi_v<scalar> * point_index / maxpoints();
        vector_t coord = {radius() * cos(t), radius() * sin(t), 0.};

        return coord;
      }
    }
  };
};
}  // namespace siconos::collision::shape
