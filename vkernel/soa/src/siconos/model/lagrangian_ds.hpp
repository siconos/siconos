#pragma once

#include "siconos/model/model_head.hpp"
#include "siconos/algebra/eigen.hpp"

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
  using lagrangian_dynamical_system_t = void;

  /// @brief Degree of freedom parameter specifying system dimension
  using dof = some::indice_parameter<"dof">;

  using batch_capable = void;

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

struct elastic_lagrangian_ds : item {
  using dof = some::indice_parameter<"dof">;

  using without_attached_storages_bindings = void;

  struct attributes {
    some::vector<some::scalar, dof> q;                 ///< Position vector
    some::vector<some::scalar, dof> velocity;          ///< Velocity vector
    some::matrix<some::scalar, dof, dof> mass_matrix;  ///< Mass matrix
    some::matrix<some::scalar, dof, dof> k_matrix;     ///< Stiffness matrix
    some::vector<some::scalar, dof> fext;  ///< External forces vector
  };

  template <typename Handle>
  struct interface : lagrangian_ds::template interface<Handle> {
    using lagrangian_ds::template interface<Handle>::self;

    decltype(auto) k_matrix() { return attr<"k_matrix">(*self()); }
  };
};

static constexpr auto has_k_matrix = mp::is_valid(
    [](auto t) -> decltype(&(decltype(t)::attributes::k_matrix)) {});

template <typename Handle>
static constexpr bool runtime_dof()
{
  using q_t = decltype(std::declval<Handle>().q());
  return match::variable_size_vector<q_t>;
}

static_assert(has_k_matrix(elastic_lagrangian_ds{}));
static_assert(!has_k_matrix(lagrangian_ds{}));

}  // namespace siconos::model
