// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#pragma once

#include "siconos/vkernel/model/model_head.hpp"

namespace siconos::vkernel::model {

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
    some::scalar e;   ///< Coefficient of restitution (0 = inelastic, 1 = fully elastic)
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
 * @brief Frémond impact law with friction
 *
 * The Newton impact law in the normal direction, Coulomb friction written
 * on the theta-weighted velocity u(k + theta) instead of u(k + 1), see
 * Acary & Collins-Craft, On the Moreau-Jean scheme with the Frémond impact
 * law, JTCAM 2025. With 1/2 <= theta <= 1 / (1 + e) the Moreau-Jean scheme
 * dissipates energy at every contact; with theta = 1 it is the Newton
 * impact law. No tangential restitution: Frémond's law is et = 1 by
 * construction. Same parameters as newton_impact_friction, and a friction
 * law for everything that tests for it.
 */
struct fremond_impact_friction : newton_impact_friction {};

/**
 * @brief Newton impact law with Coulomb friction and rolling resistance
 *
 * The rows of a contact: normal, tangent, then the rolling row (in 2D the
 * relative angular velocity omega1 - omega2, written by the relations of
 * the disks and of the polytopes, one per manifold point, as the AOS kernel
 * Lagrangian2d3DR). The cone |r_t| <= mu r_n, |m_r| <= mu_r r_n, mu_r a
 * length (Acary & Bourrier 2021, as RollingFrictionContactProblem). The
 * Newton impact law on the normal row only, none on the rolling row. With
 * mu_r = 0 it is Coulomb friction, and so is a contact whose relation has
 * no rolling row (FEM self contact: no rotational dof).
 */
struct newton_impact_rolling_friction : newton_impact_friction {
  /// Degrees of freedom for this law (3: normal + tangential + rolling)
  static constexpr auto size = 3;
  /// The rolling rows, after those of newton_impact_friction
  static constexpr auto rolling_size = 1;

  struct attributes {
    some::scalar e;     ///< Coefficient of restitution
    some::scalar mu;    ///< Friction coefficient
    some::scalar mu_r;  ///< Rolling resistance, a length
  };

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) e() { return attr<"e">(*self()); }
    decltype(auto) mu() { return attr<"mu">(*self()); }
    /// Access rolling resistance
    decltype(auto) mu_r() { return attr<"mu_r">(*self()); }
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
    some::scalar e;  ///< Coefficient of restitution (0 = inelastic, 1 = fully elastic)
  };

  /// Interface for accessing law parameter
  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    /// Access coefficient of restitution
    decltype(auto) e() { return attr<"e">(*self()); };
  };
};

}  // namespace siconos::vkernel::model
