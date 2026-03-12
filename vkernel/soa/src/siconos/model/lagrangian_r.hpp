#pragma once

#include <functional>  // For std::mem_fn

#include "siconos/algebra/algebra.hpp"
#include "siconos/model/lagrangian_ds.hpp"
#include "siconos/model/model_head.hpp"

namespace siconos::storage::pattern::match {

/// Concept for linear relations
template <typename T>
concept linear_relation = match::handle<T, model::linear>;

/// Concept for relations with one argument
template <typename T>
concept relation1 = match::handle<T, model::relation1>;

/// Concept for relations with two arguments
template <typename T>
concept relation2 = match::handle<T, model::relation2>;

}  // namespace siconos::storage::pattern::match

namespace siconos::model {

/**
 * @brief Lagrangian relation template for non-smooth laws
 *
 * @tparam NSLSize Size parameter for the non-smooth law
 *
 * This template defines a Lagrangian relation that can serve as:
 * - A linear relation
 * - A relation with one argument (relation1)
 * - A relation with two arguments (relation2)
 * - Any Lagrangian relation
 *
 * It contains attributes for:
 * - Constant term (b)
 * - Constraint matrix (h_matrix)
 */
template <auto NSLSize>
struct lagrangian_r : item,
                      linear,
                      relation1,
                      relation2,
                      any_lagrangian_relation {
  /// Size of the non-smooth law
  using nslaw_size = some::param_val<NSLSize>;

  /// Degrees of freedom parameter
  using dof = some::indice_parameter<"dof">;

  /**
   * @brief Attributes structure holding relation data
   */
  struct attributes {
    some::scalar b;  ///< Constant term in the relation
    some::matrix<some::scalar, nslaw_size, dof>
        h_matrix;  ///< Constraint matrix
  };

  /**
   * @brief Interface for Lagrangian relation operations
   * @tparam Handle Type of the handle used to access the relation
   */
  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    /// Accessor for constraint matrix
    decltype(auto) h_matrix() { return attr<"h_matrix">(*self()); }

    /// Accessor for constant term
    decltype(auto) b() { return attr<"b">(*self()); }

    /**
     * @brief Compute Jacobian for two dynamical systems
     *
     * @param step Current simulation step
     * @param ds1 First dynamical system
     * @param ds2 Second dynamical system
     * @param h_matrix1 Output Jacobian for first system
     * @param h_matrix2 Output Jacobian for second system
     */
    decltype(auto) compute_jachq(auto step, auto& ds1, auto& ds2,
                                 auto& h_matrix1, auto& h_matrix2)
    {
      h_matrix1 = h_matrix();
      h_matrix2 = -h_matrix();
    }

    /**
     * @brief Compute Jacobian for a single dynamical system
     *
     * @param step Current simulation step
     * @param ds Dynamical system
     * @param h_matrix1 Output Jacobian
     */
    decltype(auto) compute_jachq(auto step, auto& ds, auto& h_matrix1)
    {
      h_matrix1 << -h_matrix();
    }
  };
};

}  // namespace siconos::model
