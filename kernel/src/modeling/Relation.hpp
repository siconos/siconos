/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*! \file Relation.hpp
  \brief General interface for relations.
*/

#ifndef RELATION_H
#define RELATION_H
#include <memory>
#include <string>

#include "FunctionTypes.hpp"
#include "Interaction.hpp"
#include "RelationVisitor.hpp"
#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosSerialization.hpp"  // For ACCEPT_SERIALIZATION

namespace siconos::modeling {

/** @enum RelationType
 *  @brief List of authorized Relation types
 */
enum class RelationType {
  /** First Order */
  FirstOrder,
  /** Lagrangian (2nd order) with dense storage*/
  Lagrangian,
  /** Lagrangian (2nd order) with sparse storage */
  LagrangianSparse,
  /** Newton-Euler (2nd order) */
  NewtonEuler
};

/** @enum RelationSubType
 *  @brief List of authorized Relation subtypes
 */
enum class RelationSubType {
  /** non linear */
  NonLinearR,
  /** linear */
  LinearR,
  /** Linear and time invariant */
  LinearTIR,
  /** Scleronomous: depends on state (lagrangian only) */
  ScleronomousR,
  /** Rheonomous: depends on time and state (lagrangian only) */
  RheonomousR,
  /** Compliant: depends on state and \f$\lambda\f$ (lagrangian only) */
  CompliantR,
  /** CompliantLinearTIR: compliant and linear time invariant (lagrangian only) */
  CompliantLinearTIR,
  /** */
  Type1R,
  /** */
  Type2R,
  /**  */
  StressLinearTIR
};

/**
   @brief General Non Linear Relation (Abstract Base class for Relations).

   The present class is an interface to all relations and provides
   tools to define and describe them.

   A relation is a link between global variables of one or two dynamical
   systems and local variables (e.g. at contact), named y and \f$\lambda\f$
   belonging to an Interaction.

   Notations:

   \f[
   output &= y = h(...) \\
   input &= R = g(\lambda,...)
   \f]

   with \f$ y \in  \mathbb{R}^{M}, \lambda \in  \mathbb{R}^{M} \f$.

  \f$ R \in \mathbb{R}^{N} \f$ is the nonsmooth "input" to the dynamics (see r in
  DynamicalSystem). It is composed with the concatenation of the r components of the dynamical
  systems (one or two) involved in the relation.

  \f[
    J^h_x(...)
    = \frac{\partial h}{\partial x} \in \mathbb{R}^{M \times N}
    =
   \left[
   \frac{\partial h_i}{\partial x_j}
   \right]_{1\le i\le M,\;1\le j\le N}.
  \f]

   with \f$ x \in \mathbb{R}^{N} \f$

  - All relations are specified by their type (First order, Lagrangian, NewtonEuler, see
  RelationType) accessed by getType() and their sub-type (linear, scleronomous ..., see
  RelationSubType), returned by getSubType().

  - The variables to be taken into account for g and h, as well as their Jacobians, depend on
  the specific type of relationship.

  - A relation provides functions to compute:

      - the "output" (function computeOutput()) to update y using dynamical systems global
        variables

      - the "input" to the dynamics, (function computeInput()) to update the non-smooth
  dynamical systems parts, (e.g. r or p).

      - the jacobians of the output and of the input.

*/
class Relation {
 protected:
  ACCEPT_SERIALIZATION(Relation);

  /** type of the Relation: FirstOrder or Lagrangian */
  RelationType _relationType;

  /** sub-type of the Relation (exple: LinearTIR or ScleronomousR ...) */
  RelationSubType _subType;

  /** base and only constructor
   *
   *  @param type type of the relation
   *  @param subtype subtype of the relation
   */
  Relation(RelationType type, RelationSubType subtype)
      : _relationType(type), _subType(subtype) {};

 private:
  // Rule of five ...
  Relation() = delete;
  Relation(const Relation&) = delete;
  Relation(Relation&&) = delete;
  Relation& operator=(const Relation&) = delete;
  Relation& operator=(Relation&&) = delete;

 public:
  /** destructor */
  virtual ~Relation() noexcept = default;

  /** @return the type of the Relation */
  inline RelationType getType() const { return _relationType; }

  /** @return the subType of the Relation */
  inline RelationSubType getSubType() const { return _subType; }

  /** initialize the relation (check sizes, memory allocation ...)
   *
   *  @param inter the interaction using this relation
   */
  virtual void initialize(Interaction& inter) = 0;

  /** check sizes of the relation specific operators.
   *
   *  @param inter an Interaction using this relation
   */
  virtual void checkSize(const Interaction& inter) const = 0;

  /** compute all the Jacobian of the output (h(...)) for the current relation
   *
   *  @param time the current time
   *  @param inter the interaction using this relation
   */
  virtual void computeJach(double time, Interaction& inter) = 0;

  /** compute all the Jacobian of the input (r(...) or p(...)) for the current relation
   *
   *  @param time the current time
   *  @param inter the interaction using this relation
   */
  virtual void computeJacg(double time, Interaction& inter) {
    /*Does nothing by default. Reimplement if required*/
  };

  /** default function to compute y or one of its derivative
   *
   *  @param time the current time
   *  @param inter the interaction using this relation
   *  @param derivativeNumber number of the derivative to compute (optional,
   *  default = 0)
   */
  virtual void computeOutput(double time, Interaction& inter,
                             siconos::algebra::blocks::size_type derivativeNumber = 0) = 0;
  /** default function to compute r or one of its derivative
   *
   *  @param time the current time
   *  @param inter the interaction using this relation
   *  @param level the input "derivative" order of lambda used to compute input
   */
  virtual void computeInput(double time, Interaction& inter,
                            siconos::algebra::blocks::size_type level = 0) = 0;

  /** @return true if the relation is linear */
  virtual bool isLinear() const { return false; }

  /** @return true if the relation requires the computation of residu */
  virtual bool requireResidu() { return false; }

  /** main relation members display */
  virtual void display() const = 0;

  /** @return True if \f$ J^h_\lambda h(...) \f$ is taken into account */
  virtual bool hasJacobianhOver_lambda() const { return false; }

  /** @return a read - only view on the matrix \f$ J^h_{\lambda} h(...) \f$

    warning: use hasJacobianhOver_lambda before any call to ensure the jacobian is defined and
    has sense */
  virtual const siconos::algebra::ConstMapType jacobianhOver_lambda() const {
    // To be overriden in classes where jacobianhOver_lambda has sense
    THROW_EXCEPTION("jacobian h over lambda is not defined for this kind of relation.");
  }

  /** @brief allocation of memory space for relations in the graph
   *      Warning: internal use only (called from Topology)
   *  @param dslink a container of vectors (pointers), from the parent interaction
   *  @param ds1 first ds concerned by the relation
   *  @param ds2 second ds concerned by the relation
   */
  virtual void allocate_read_dynamical_systems_var_vectors(
      std::vector<std::shared_ptr<siconos::algebra::BlockVector>>& ds_vars,
      modeling::DynamicalSystem& ds1, modeling::DynamicalSystem& ds2) const = 0;

  /**
   * @brief Accept a visitor.
   *
   * Dispatch of the Visitor: forwards this object to the appropriate
   * `relations::Visitor::visit()` (corresponding to the concrete relation type)
   *
   * @param visitor Visitor receiving this relation.
   */
  virtual void accept(relations::Visitor&) const {
    throw std::logic_error("accept (relation): no visitor defined");
  }
};
}  // namespace siconos::modeling

#endif  // RELATION_H
