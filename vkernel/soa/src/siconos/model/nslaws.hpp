#pragma once

#include "siconos/model/model_head.hpp"

namespace siconos::model {

/// Represents an equality condition constraint (placeholder)
struct equality_condition {};

/// Represents a relay constraint (placeholder)
struct relay {};

/**
 * @brief Non-smooth law for Newton impact with friction
 *
 * Models a Newton impact law with Coulomb friction constraints.
 * Contains parameters for coefficient of restitution (e) and friction
 * coefficient (mu).
 */
struct newton_impact_friction : item {
  /// Degrees of freedom for this law (2: normal + tangential)
  static constexpr auto size = 2;

  /// Attributes defining the law's parameters
  struct attributes {
    some::scalar
        e;  ///< Coefficient of restitution (0 = inelastic, 1 = fully elastic)
    some::scalar mu;  ///< Friction coefficient (Coulomb friction)
  };

  /// Interface for accessing law parameters
  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    /// Access coefficient of restitution
    decltype(auto) e() { return attr<"e">(*self()); }

    /// Access friction coefficient
    decltype(auto) mu() { return attr<"mu">(*self()); }
  };
};

/**
 * @brief Newton impact law without friction
 *
 * Models a simple Newton impact law with only normal restitution.
 * Contains parameter for coefficient of restitution (e).
 */
struct newton_impact : item {
  /// Degrees of freedom for this law (1: normal direction only)
  static constexpr auto size = 1;

  /// Attributes defining the law's parameter
  struct attributes {
    some::scalar
        e;  ///< Coefficient of restitution (0 = inelastic, 1 = fully elastic)
  };

  /// Interface for accessing law parameter
  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    /// Access coefficient of restitution
    decltype(auto) e() { return attr<"e">(*self()); };
  };
};

}  // namespace siconos::model
