#pragma once

#include <functional>

#include "siconos/algebra/algebra.hpp"
#include "siconos/model/lagrangian_ds.hpp"
#include "siconos/model/model_head.hpp"

namespace siconos::storage::pattern::match {

/**
 * @brief Concept for linear relation handles
 * @tparam T Type to check
 *
 * Requires that T is a handle to a model::linear type.
 * Used to identify relations that have linear constraint matrices.
 */
template <typename T>
concept linear_relation = match::handle<T, model::linear>;

/**
 * @brief Concept for single-body relation handles
 * @tparam T Type to check
 *
 * Requires that T is a handle to a model::relation1 type.
 * Used for relations involving a single dynamical system (e.g., contact with
 * fixed obstacle).
 */
template <typename T>
concept relation1 = match::handle<T, model::relation1>;

/**
 * @brief Concept for two-body relation handles
 * @tparam T Type to check
 *
 * Requires that T is a handle to a model::relation2 type.
 * Used for relations involving two dynamical systems (e.g., contact between
 * two bodies).
 */
template <typename T>
concept relation2 = match::handle<T, model::relation2>;

}  // namespace siconos::storage::pattern::match

namespace siconos::model {

/**
 * @brief Lagrangian relation for non-smooth laws with compile-time constraint
 * size
 * @tparam NSLSize Number of constraints (compile-time constant)
 *
 * Represents a linear constraint relation of the form h(q) = H*q + b, where:
 * - H is the constraint matrix (h_matrix)
 * - b is the constant offset term
 * - q is the generalized coordinates
 *
 * This class supports both unilateral constraints (single body vs obstacle)
 * and bilateral constraints (two bodies). The constraint size is specified at
 * compile time via the template parameter NSLSize.
 *
 * Inherits from:
 * - item: Base class for storage pattern items
 * - linear: Marker for linear relations
 * - relation1: Single-body relation capability
 * - relation2: Two-body relation capability
 * - any_lagrangian_relation: General Lagrangian relation marker
 */
template <auto NSLSize>
struct lagrangian_r : item,
                      linear,
                      relation1,
                      relation2,
                      any_lagrangian_relation {
  /**
   * @brief Type representing the non-smooth law size (number of constraints)
   *
   * Wraps the compile-time constant NSLSize as a type for use in attribute
   * definitions.
   */
  using nslaw_size = some::param_val<NSLSize>;

  /**
   * @brief Type representing the degrees of freedom parameter
   *
   * Used to specify the dimension of the configuration space for the
   * constraint matrix. The actual value is provided by the environment
   * parameters.
   */
  using dof = some::indice_parameter<"dof">;

  /**
   * @brief Attributes of the Lagrangian relation
   *
   * Defines the storage layout for the constant term b and the constraint
   * matrix h_matrix. Both are sized according to the compile-time constraint
   * dimension NSLSize.
   */
  struct attributes {
    /**
     * @brief Constant term b in the constraint equation h(q) = H*q + b
     *
     * Vector of size NSLSize containing the constant offset of the
     * constraint. For contact constraints, this typically represents the
     * initial gap or penetration.
     */
    some::vector<some::scalar, nslaw_size> b;

    /**
     * @brief Constraint matrix H in the constraint equation h(q) = H*q + b
     *
     * Matrix of size NSLSize x dof containing the linear coefficients
     * relating generalized coordinates to constraint values. Each row
     * corresponds to one constraint equation.
     */
    some::matrix<some::scalar, nslaw_size, dof> h_matrix;
  };

  /**
   * @brief Interface providing methods to access and manipulate the relation
   * @tparam Handle The handle type used to access this item in storage
   *
   * Provides methods to get/set the constraint matrix and constant term,
   * as well as compute Jacobians for the non-smooth solver.
   */
  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    /**
     * @brief Get the constraint matrix H
     * @return Reference to the constraint matrix h_matrix
     *
     * Returns a reference to the matrix H such that h(q) = H*q + b.
     * The matrix has dimensions NSLSize x dof.
     */
    decltype(auto) h_matrix() { return attr<"h_matrix">(*self()); }

    /**
     * @brief Get the constant term b
     * @return Reference to the constant vector b
     *
     * Returns a reference to the vector b such that h(q) = H*q + b.
     * The vector has length NSLSize.
     */
    decltype(auto) b() { return attr<"b">(*self()); }

    /**
     * @brief Compute Jacobians for bilateral constraint between two bodies
     * @tparam StepType Type of the simulation step index
     * @tparam DS1Type Type of the first dynamical system handle
     * @tparam DS2Type Type of the second dynamical system handle
     * @tparam MatrixType Type of the output Jacobian matrices
     * @param[in] step Current simulation step index (unused but required by
     * interface)
     * @param[in] ds1 Handle to the first dynamical system
     * @param[in] ds2 Handle to the second dynamical system
     * @param[out] h_matrix1 Jacobian matrix for the first system (filled with
     * H)
     * @param[out] h_matrix2 Jacobian matrix for the second system (filled
     * with -H)
     *
     * Computes the constraint Jacobians for a bilateral constraint of the
     * form: h(q1, q2) = H*q1 - H*q2 + b = 0
     *
     * The Jacobian with respect to q1 is H, and with respect to q2 is -H.
     * This is used for contact constraints between two moving bodies.
     */
    decltype(auto) compute_jachq(auto step, auto& ds1, auto& ds2,
                                 auto& h_matrix1, auto& h_matrix2)
    {
      h_matrix1 = h_matrix();
      h_matrix2 = -h_matrix();
    }

    /**
     * @brief Compute Jacobian for unilateral constraint with fixed obstacle
     * @tparam StepType Type of the simulation step index
     * @tparam DSType Type of the dynamical system handle
     * @tparam MatrixType Type of the output Jacobian matrix
     * @param[in] step Current simulation step index (unused but required by
     * interface)
     * @param[in] ds Handle to the dynamical system
     * @param[out] h_matrix1 Jacobian matrix (filled with -H)
     *
     * Computes the constraint Jacobian for a unilateral constraint of the
     * form: h(q) = -H*q + b >= 0
     *
     * The negative sign indicates that the constraint is typically formulated
     * as a gap function that must remain non-negative (e.g., distance to
     * obstacle). This is used for contact constraints between a body and a
     * fixed obstacle.
     */
    decltype(auto) compute_jachq(auto step, auto& ds, auto& h_matrix1)
    {
      h_matrix1 << -h_matrix();
    }
  };
};

/**
 * @brief Runtime-sized Lagrangian relation for non-smooth laws
 *
 * Represents a linear constraint relation h(q) = H*q + b where the constraint
 * size is determined at runtime rather than compile time. This is necessary
 * for applications with adaptive discretizations or FEM where the number of
 * active constraints varies during simulation.
 *
 * Unlike lagrangian_r, this class uses unbounded (dynamically sized) storage
 * for the constraint matrix and constant vector, allowing the constraint
 * dimension to change at runtime.
 *
 * Currently only supports single-body relations (relation1).
 */
struct rt_lagrangian_r : item, model::linear, model::relation1 {
  /**
   * @brief Attributes of the runtime Lagrangian relation
   *
   * Defines the storage layout for dynamically sized constraint data.
   * The sizes are determined at runtime based on the specific constraint
   * configuration.
   */
  struct attributes {
    /**
     * @brief Constraint matrix H with runtime-determined dimensions
     *
     * Unbounded matrix storing the linear coefficients H such that h(q) = H*q
     * + b. Dimensions are NSLSize x dof, determined at runtime.
     */
    some::unbounded_matrix<some::scalar> h_matrix;

    /**
     * @brief Constant term b with runtime-determined size
     *
     * Unbounded vector storing the constant offset b such that h(q) = H*q +
     * b. Length is NSLSize, determined at runtime.
     */
    some::unbounded_vector<some::scalar> b;
  };

  /**
   * @brief Interface providing methods to access and manipulate the runtime
   * relation
   * @tparam Handle The handle type used to access this item in storage
   *
   * Provides methods to get/set the constraint matrix and constant term,
   * compute constraint values, and compute Jacobians.
   */
  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    /**
     * @brief Get the constraint matrix H
     * @return Reference to the unbounded constraint matrix h_matrix
     *
     * Returns a reference to the runtime-sized matrix H.
     * The matrix dimensions are determined by the current constraint
     * configuration.
     */
    decltype(auto) h_matrix() { return attr<"h_matrix">(*self()); }

    /**
     * @brief Get the constant term b
     * @return Reference to the unbounded constant vector b
     *
     * Returns a reference to the runtime-sized vector b.
     * The vector length is determined by the current constraint
     * configuration.
     */
    decltype(auto) b() { return attr<"b">(*self()); }

    /**
     * @brief Compute the constraint value h(q) = H*q + b
     * @tparam DSType Type of the dynamical system handle
     * @param[in] ds Handle to the dynamical system
     * @return Vector containing the constraint values h(q)
     *
     * Evaluates the constraint function for the given system configuration.
     * Returns a vector where each element corresponds to one constraint
     * equation. Positive values typically indicate penetration (for contact
     * constraints).
     */
    auto compute_h(auto& ds)
    {
      return h_matrix() * storage::attr<"q">(ds) + b();
    }

    /**
     * @brief Compute the constraint Jacobian dh/dq = H
     * @tparam StepType Type of the simulation step index
     * @tparam DSType Type of the dynamical system handle
     * @tparam MatrixType Type of the output Jacobian matrix
     * @param[in] step Current simulation step index (unused but required by
     * interface)
     * @param[in] ds Handle to the dynamical system
     * @param[out] h_mat Output matrix filled with the constraint Jacobian H
     *
     * Computes the Jacobian of the constraint function with respect to the
     * generalized coordinates. For linear constraints, this is simply the
     * constraint matrix H.
     */
    void compute_jachq(auto step, auto& ds, auto& h_mat)
    {
      h_mat = h_matrix();
    }
  };
};

}  // namespace siconos::model
