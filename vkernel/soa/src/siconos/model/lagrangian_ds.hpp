#pragma once

#include "siconos/model/model_head.hpp"

namespace siconos::model {
/**
 * @brief Represents a Lagrangian dynamical system with fixed-size attributes
 *
 * This struct defines a Lagrangian dynamical system with attributes related
 * to its configuration, velocity, mass matrix, and external forces. The
 * size of vectors and matrices is fixed at compile time through the "dof"
 * parameter.
 */
struct lagrangian_ds : item {
  /// @brief Degree of freedom parameter specifying system dimension
  using dof = some::indice_parameter<"dof">;

  /**
   * @brief Attributes defining the state and properties of the system
   *
   * Contains:
   * - q: Position vector (size = dof)
   * - velocity: Velocity vector (size = dof)
   * - mass_matrix: Mass matrix (size = dof x dof)
   * - fext: External forces vector (size = dof)
   */
  struct attributes {
    some::vector<some::scalar, dof> q;                 ///< Position vector
    some::vector<some::scalar, dof> velocity;          ///< Velocity vector
    some::matrix<some::scalar, dof, dof> mass_matrix;  ///< Mass matrix
    some::vector<some::scalar, dof> fext;  ///< External forces vector
  };

  /**
   * @brief Interface for accessing system attributes
   * @tparam Handle Type providing access to the system instance
   */
  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;  ///< Inherited self pointer

    /// @brief Access mass matrix attribute
    decltype(auto) mass_matrix() { return attr<"mass_matrix">(*self()); }

    /// @brief Access velocity attribute
    decltype(auto) velocity() { return attr<"velocity">(*self()); }

    /// @brief Access position attribute
    decltype(auto) q() { return attr<"q">(*self()); }

    /// @brief Access external forces attribute
    decltype(auto) fext() { return attr<"fext">(*self()); }
  };
};

/**
 * @brief Represents a Lagrangian dynamical system with runtime-sized
 * attributes
 *
 * Similar to lagrangian_ds but with unbounded vectors and matrices whose
 * sizes are determined at runtime. Includes an additional "dof" attribute
 * to specify the current system dimension.
 */
struct rt_lagrangian_ds : item {
  /// @brief Indicates no attached storage bindings (optimization hint)
  using without_attached_storages_bindings = void;

  /**
   * @brief Attributes defining the state and properties of the system
   *
   * Contains:
   * - dof: Current degree of freedom (system dimension)
   * - q: Position vector (unbounded size)
   * - velocity: Velocity vector (unbounded size)
   * - mass_matrix: Mass matrix (unbounded size)
   * - fext: External forces vector (unbounded size)
   */
  struct attributes {
    some::scalar dof;                        ///< Current degree of freedom
    some::unbounded_vector<some::scalar> q;  ///< Position vector
    some::unbounded_vector<some::scalar> velocity;     ///< Velocity vector
    some::unbounded_matrix<some::scalar> mass_matrix;  ///< Mass matrix
    some::unbounded_matrix<some::scalar> k_matrix;     // Rigidity matrix
    some::unbounded_vector<some::scalar> fext;  ///< External forces vector
  };

  /**
   * @brief Interface for accessing system attributes
   * @tparam Handle Type providing access to the system instance
   */
  template <typename Handle>
  struct interface : lagrangian_ds::interface<Handle> {
    using lagrangian_ds::interface<Handle>::self;  ///< Inherited self pointer

    /// @brief Access k_matrix
    decltype(auto) k_matrix() { return attr<"k_matrix">(*self()); }

    /// @brief Access degree of freedom attribute
    decltype(auto) dof() { return attr<"dof">(*self()); }
  };
};

// Check if a k_matrix member function is present.
static constexpr auto has_k_matrix =
    mp::is_valid([](auto &&x) -> decltype(x.k_matrix()) {});

}  // namespace siconos::model
