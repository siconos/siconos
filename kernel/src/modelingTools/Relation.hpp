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
#include "RelationVisitor.hpp"
#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosSerialization.hpp"  // For ACCEPT_SERIALIZATION
#include "SiconosVector.hpp"
#include "Interaction.hpp"

namespace siconos::modeling {

// class Interaction;

/** List of possible Relations types*/
enum class RelationType {
  /** First Order */
  FirstOrder,
  /** Lagrangian */
  Lagrangian,
  /** Lagrangian */
  NewtonEuler
};

/** List of possible Relations subtypes*/
enum class RelationSubType {
  /** non linear */
  NonLinearR,
  /** linear */
  LinearR,
  /** Linear and time invariant */
  LinearTIR,
  /** Scleronomous (lagrangian only) */
  ScleronomousR,
  /** Rheonomous (lagrangian only) */
  RheonomousR,
  /** Compliant (lagrangian only) */
  CompliantR,
  /** CompliantLinearTIR (lagrangian only) */
  CompliantLinearTIR,
  /** */
  Type1R,
  /** */
  Type2R
};

/**
   General Non Linear Relation (Abstract Base class for Relations).

   The present class is an interface to all relations and provides
   tools to define and describe them.

   A relation is a link between global variables of one or two dynamical
   systems and local variables (e.g. at contact), named y and lambda
   belonging to an interaction.

   All relations are specified by their type (First order or Lagrangian)
   accessed by getType() and their sub-type (linear, scleronomous ...),
   returned by getSubType().

   A relation provides functions to compute:

   - the "output" (function computeOutput()) to update y using dynamical systems global
   variables

   - the "input" (function computeInput()) to update the non-smooth dynamical systems parts
   (e.g. r or p) using  \f$ \lambda \f$ .

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
   *  \param type type of the relation
   *  \param subtype subtype of the relation
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

  /** \return the type of the Relation (FirstOrder or Lagrangian) */
  inline RelationType getType() const { return _relationType; }

  /** \return the subType of the Relation */
  inline RelationSubType getSubType() const { return _subType; }

  /** initialize the relation (check sizes, memory allocation ...)
   *
   *  \param inter the interaction using this relation
   */
  virtual void initialize(Interaction& inter) = 0;

  /** check sizes of the relation specific operators.
   *
   *  \param inter an Interaction using this relation
   */
  virtual void checkSize(const Interaction& inter) const = 0;

  /** compute all the H Jacobian
   *
   *  \param time the current time
   *  \param inter the interaction using this relation
   */
  virtual void computeJach(double time, Interaction& inter) = 0;

  /** compute all the G Jacobian
   *
   *  \param time the current time
   *  \param inter the interaction using this relation
   *  \param interProp
   */
  virtual void computeJacg(
      double time,
      Interaction& inter) { /*Does nothing by default. Reimplement if required*/ };

  /** default function to compute y
   *
   *  \param time the current time
   *  \param inter the interaction using this relation
   *  \param derivativeNumber number of the derivative to compute (optional,
   *  default = 0)
   */
  virtual void computeOutput(double time, Interaction& inter,
                             unsigned int derivativeNumber = 0) = 0;
  /** default function to compute r
   *
   *  \param time the current time
   *  \param inter the interaction using this relation
   *  \param level the input "derivative" order of lambda used to compute input
   */
  virtual void computeInput(double time, Interaction& inter, unsigned int level = 0) = 0;

  /** \return true if the relation is linear */
  virtual bool isLinear() const { return false; }

  /** \return true if the relation requires the computation of residu */
  virtual bool requireResidu() { return false; }

  /** main relation members display */
  virtual void display() const = 0;

  /** \return True if \f$ \nabla_x h(x,t,\lambda) \f$ is taken into account */
  virtual bool hasJacobianhOver_lambda() const { return false; }

  /*  \return a read-only view on the matrix \f$ \nabla^\top_{\lambda}h(q,\lambda) \f$

      warning: use hasJacobianhOver_lambda before any call to ensure the jacobian is
      defined and has sense
  */
  virtual const siconos::algebra::ConstMapType jacobianhOver_lambda() const {
    // To be overriden in classes where jacobianhOver_lambda has sense
    THROW_EXCEPTION("jacobian h over lambda is not defined for this kind of relation.");
    // handle return value to avoid compiler warning.
    // static const siconos::algebra::SiconosMatrix empty_matrix;
    // static const siconos::algebra::ConstMapType empty_map(empty_matrix.data(), 0, 0);
    // return empty_map;
  }

  virtual void accept(relations::Visitor&) const {
    throw std::logic_error("accept (relation): no visitor defined");
  }
};
}  // namespace siconos::modeling

#endif  // RELATION_H
