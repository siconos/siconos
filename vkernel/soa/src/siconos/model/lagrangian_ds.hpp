#pragma once

#include "siconos/model/model_head.hpp"

namespace siconos::model {

/**
 * @brief Represents a Lagrangian dynamical system.
 *
 * This struct defines a Lagrangian dynamical system with attributes related
 * to its configuration, velocity, mass matrix, and external forces.
 */
struct lagrangian_ds
    : item {
  using dof =
      some::indice_parameter<"dof">; /**< Degree of freedom parameter. */

  /**
   * @brief Defines the attributes of the Lagrangian dynamical system.
   *
   * The attributes include:
   * - q: position vector
   * - velocity: velocity vector
   * - mass_matrix: mass matrix of the system
   * - fext: external forces vector
   */
  using attributes =
      gather<attribute<"q", some::vector<some::scalar, dof>>,
             attribute<"velocity", some::vector<some::scalar, dof>>,
             attribute<"mass_matrix", some::matrix<some::scalar, dof, dof>>,
             attribute<"fext", some::vector<some::scalar, dof>>>;

  /**
   * @brief Interface for interacting with the Lagrangian dynamical system.
   *
   * This nested template struct provides methods to access the system's
   * attributes.
   *
   * @tparam Handle Type of the handle for accessing the interface.
   */
  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self; /**< Inherit the `self` pointer
                                              from the parent. */

    /**
     * @brief Access the mass matrix attribute.
     *
     * @return Reference to the mass matrix.
     */
    decltype(auto) mass_matrix() { return attr<"mass_matrix">(*self()); }

    /**
     * @brief Access the velocity attribute.
     *
     * @return Reference to the velocity vector.
     */
    decltype(auto) velocity() { return attr<"velocity">(*self()); }

    /**
     * @brief Access the position attribute.
     *
     * @return Reference to the position vector.
     */
    decltype(auto) q() { return attr<"q">(*self()); }

    /**
     * @brief Access the external forces attribute.
     *
     * @return Reference to the external forces vector.
     */
    decltype(auto) fext() { return attr<"fext">(*self()); }
  };
};

/**
 * @brief Represents a Lagrangian dynamical system with the number of
 *        degrees of freedom specified at runtime.
 *
 * This struct provides an alternative representation of a Lagrangian
 * dynamical system with unbounded attributes.
 */
struct rt_lagrangian_ds : item {
  using without_attached_storages_bindings =
      void; /**< No attached storage bindings. */

  /**
   * @brief Defines the unbounded attributes of the Lagrangian
   * dynamical system.
   *
   * The attributes include:
   * - q: position vector (unbounded)
   * - velocity: velocity vector (unbounded)
   * - mass_matrix: mass matrix of the system (unbounded)
   * - fext: external forces vector (unbounded)
   */
  using attributes =
      gather<attribute<"dof", some::scalar>,
             attribute<"q", some::unbounded_vector<some::scalar>>,
             attribute<"velocity", some::unbounded_vector<some::scalar>>,
             attribute<"mass_matrix", some::unbounded_matrix<some::scalar>>,
             attribute<"fext", some::unbounded_vector<some::scalar>>>;

  template <typename Handle>
  struct interface : lagrangian_ds::interface<Handle> {
    using lagrangian_ds::interface<Handle>::self;

    decltype(auto) dof() { return attr<"dof">(*self()); }
  };
};

}  // namespace siconos::model
