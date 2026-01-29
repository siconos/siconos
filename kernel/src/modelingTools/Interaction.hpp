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

/*! \file Interaction.hpp
  \brief Interaction class and related typedef
*/

#ifndef INTERACTION_H
#define INTERACTION_H

#include <assert.h>

#include <memory>
#include <vector>

#include "SiconosMatrix.hpp"
#include "SiconosSerialization.hpp"
#include "SiconosVector.hpp"

namespace siconos::algebra {
class SiconosMemory;
class BlockVector;
}  // namespace siconos::algebra

namespace siconos::modeling {

class NonSmoothLaw;
class Relation;
class DynamicalSystem;

/**
    Description of a non-smooth interaction.
    The object Interaction is used to defined a "link" between one or two DynamicalSystem,
    like unilateral constraints and some nonsmooth law (e.g. complementarity).

    It holds two vectors of "local" variables, \f$ y \f$ and \f$ \lambda \f$
    and their derivatives, which are related to the state variables and the inputs of the
    DynamicalSystem (x,R) through constraints defined in a Relation and completed with
    a NonSmoothLaw involving those variables.

    Remarks:

    - one and only one Relation (access: relation()) per Interaction
    - one and only one NonSmoothLaw (access: nonSmoothLaw()) per Interaction
    - dimension() is the size of the interaction and so the size of vectors y, lambda
    and their derivatives.
    - output: y(i), to get derivative i of y
    - input: lambda(i), to get derivative i of lambda

 */
class Interaction : public std::enable_shared_from_this<Interaction> {
 private:
  ACCEPT_SERIALIZATION(Interaction);

  /* internal counter used to set interaction number */
  static size_t count_;

  /** Interaction id */
  size_t _number{0};

  /** Minimum required 'level' for output y
   *  y will be initialized from
   *  y[_lowerLevelForOutput] to y[_upperLevelForOutput]
   */
  unsigned int _lowerLevelForOutput = 0;

  /** Maximum required 'level' for output y
   *  y will be initialized from
   *  y[_lowerLevelForOutput] to y[_upperLevelForOutput]
   */
  unsigned int _upperLevelForOutput = 0;

  /** Minimum required 'level' for input lambda
   *  lambda will be initialized from
   *  lambda[_lowerLevelForIntput] to lambda[_upperLevelForInput]
   */
  unsigned int _lowerLevelForInput = 0;

  /** Maximum required 'level' for input lambda
   *  lambda will be initialized from
   *  lambda[_lowerLevelForIntput] to lambda[_upperLevelForInput]
   */
  unsigned int _upperLevelForInput = 0;

  /** size of the interaction, ie size of y[i] and _lambda[i] */
  siconos::algebra::Index _interactionSize = 0;

  /** sum of all DS sizes, for DS involved in the interaction */
  siconos::algebra::Index _sizeOfDS = 0;

  /** Bool to check the number of DS concerned by this interaction
      (1 or 2 indeed)
      True if 2 DS.
      Note FP : usefull in NewtonEuler jacobians computation.
  */
  bool _has2Bodies = false;

  /** relation between constrained variables and states variables
   * vector of output derivatives
   * y[0] is y, y[1] is yDot and so on
   */
  std::vector<std::shared_ptr<siconos::algebra::SiconosVector>> _y = {};

  /** memory of previous coordinates of the system */
  std::vector<siconos::algebra::SiconosMemory> _yMemory;

  /** result of the computeInput function */
  std::vector<std::shared_ptr<siconos::algebra::SiconosVector>> _lambda = {};

  /** memory of previous coordinates of the system */
  std::vector<siconos::algebra::SiconosMemory> _lambdaMemory;

  /** the Non-smooth Law of the interaction*/
  std::shared_ptr<NonSmoothLaw> _nslaw{nullptr};

  /** the type of Relation of the interaction */
  std::shared_ptr<Relation> _relation{nullptr};

  /** pointer links to DS variables needed for computation,
   *  mostly used in Relations (computeOutput and computeInput)
   * and OneStepIntegrator classes. */
  std::vector<std::shared_ptr<siconos::algebra::BlockVector>> _linkToDSVariables = {};

  std::vector<std::shared_ptr<siconos::algebra::SiconosVector>> _relationVectors = {};

  // internal struct used to handle visitors process to set Interaction levels
  // depending on the nslaw and the relation.
  struct SetLevels;

  // === PRIVATE FUNCTIONS ===

  /* forbid default constructor, copy and assignment */
  Interaction(const Interaction& inter) = delete;
  Interaction& operator=(const Interaction&) = delete;

 protected:
  Interaction() = default; /* for serialization only */

 public:
  /**
     Interaction constructor

     \param NSL pointer object describing the nonsmooth law;
     the interaction size if infered from the size of this law.
     \param rel a pointer object describing the functions used to compute the constraints
  */
  Interaction(std::shared_ptr<NonSmoothLaw> NSL, std::shared_ptr<Relation> rel);

  /** destructor  */
  ~Interaction() noexcept = default;

  /**
     Update interactions attributes.
     Must be called when levels have been modified.
  */
  void reset();

  /**
      set the links  between the interaction and the DynamicalSystem(s) members.

      \param ds1 first ds linked to this Interaction (i.e IG->vertex.source)
      \param ds2 second ds linked to this Interaction (i.e IG->vertex.target) ds1 == ds2 is
     allowed
  */
  void initializeLinkToDsVariables(DynamicalSystem& ds1, DynamicalSystem& ds2);

  /** set all lambda to zero */
  void resetAllLambda();

  /** set lambda to zero for a given level
   *
   *  \param level
   */
  void resetLambda(unsigned int level);

  /** build memories vectors for y and \f$ \lambda \f$
   *
   *  \param computeResiduY true if interaction should compute extra residu value
   *  \param steps number of required memories (depends on the OSI)
   */
  void initializeMemory(unsigned int steps);

  // === GETTERS/SETTERS ===
  /** \return the id of the interaction */
  inline size_t number() const { return _number; }

  /** Set the lower level for output y.
   *
   *  \param newVal new level
   */
  inline void setLowerLevelForOutput(const unsigned int newVal) {
    _lowerLevelForOutput = newVal;
  };

  /** Set the upper level for output y.
   *
   *  \param newVal :new level
   */
  inline void setUpperLevelForOutput(const unsigned int newVal) {
    _upperLevelForOutput = newVal;
  };

  /** \return the lower level for output y. */
  inline unsigned int lowerLevelForOutput() { return _lowerLevelForOutput; };

  /**\return the upper level for output y. */
  inline unsigned int upperLevelForOutput() { return _upperLevelForOutput; };

  /** Set the lower level for input Lambda.
   *
   *  \param newVal : new level
   */
  inline void setLowerLevelForInput(const unsigned int newVal) {
    _lowerLevelForInput = newVal;
  };

  /** Set the upper level for input Lambda.
   *
   *  \param newVal new level
   */
  inline void setUpperLevelForInput(const unsigned int newVal) {
    _upperLevelForInput = newVal;
  };

  /**\return the lower level for input Lambda */
  inline unsigned int lowerLevelForInput() { return _lowerLevelForInput; };

  /**\return the upper level for input Lambda  */
  inline unsigned int upperLevelForInput() { return _upperLevelForInput; };

  /** returns dimension (i.e. nslaw size == y and lambda size) */
  inline auto dimension() const { return _interactionSize; }

  /** \return the sum of DS sizes, for DS involved in interaction */
  inline auto getSizeOfDS() const { return _sizeOfDS; }

  /**
     Set the number of dynamical systems concerned by
     this interaction. Warning FP: this function is supposed
     to be called only during topology->link(inter, ds1, ds2) call.

     \param val : true if two ds, else false
   */
  void setHas2Bodies(bool val) { _has2Bodies = val; }

  /**
     Check the number of dynamical systems concerned by
     this interaction.

     \return bool : true if two ds, else false
   */
  bool has2Bodies() const { return _has2Bodies; }

  // -- y --

  /** Get y[i], derivative number i of output
   *
   *  \param i : the derivative number
   *  \return BlockVector
   */
  const siconos::algebra::SiconosVector getCopyOfy(const unsigned int i) const;

  /** get vector of output derivatives
   *
   *  \return a std::vector<std::shared_ptr<siconos::algebra::SiconosVector>>
   */
  inline const std::vector<std::shared_ptr<siconos::algebra::SiconosVector>> y() const {
    return _y;
  }

  /** get y[i], derivative number i of output
   *
   *  \param i derivative number i of output
   *  \return pointer on a SiconosVector
   */
  inline std::shared_ptr<siconos::algebra::SiconosVector> y(const unsigned int i) const {
    assert(_y[i]);
    return _y[i];
  }

  /** get y[i], derivative number i of output
   *
   *  \param i derivative number i of output
   *  \return reference on a SiconosVector
   */
  inline siconos::algebra::SiconosVector& y_python(const unsigned int i) const {
    assert(_y[i]);
    return *(_y[i]);
  }

  /** set the output vector y to newVector with copy of the y[i] (ie
   *  memory allocation)
   *
   *  \param v std::vector<std::shared_ptr<siconos::algebra::SiconosVector>>
   */
  void setY(const std::vector<std::shared_ptr<siconos::algebra::SiconosVector>>& v);

  /** set the output vector y to newVector with direct pointer
   *  equality for the y[i]
   *
   *  \param v std::vector<std::shared_ptr<siconos::algebra::SiconosVector>>
   */
  void setYPtr(const std::vector<std::shared_ptr<siconos::algebra::SiconosVector>>& v);

  /** set y[i]
   *
   *  \param i derivative number i of output
   *  \param v the new yi value
   */
  void setY(const unsigned int i, const siconos::algebra::SiconosVector& v);

  /** set y[i] to pointer newPtr
   *
   *  \param i derivative number i of output
   *  \param v the new value for yi
   */
  void setYPtr(const unsigned int i, std::shared_ptr<siconos::algebra::SiconosVector> v);

  /** get all the values of the output y stored in memory
   *
   *  \param level
   *  \return a memory
   */
  siconos::algebra::SiconosMemory& yMemory(unsigned int level);

  /** get the last value of the  output y stored in memory
   *
   *  \param level
   *  \return a SiconosVector reference
   */
  const siconos::algebra::SiconosVector& y_k(const unsigned int i) const;

  // -- _lambda --

  /** get vector of input derivatives
   *
   *  \return a std::vector<std::shared_ptr<siconos::algebra::SiconosVector>>
   */
  inline const std::vector<std::shared_ptr<siconos::algebra::SiconosVector>> getLambda()
      const {
    return _lambda;
  }

  /** get _lambda[i], derivative number i of input
   *
   *  \param i derivative number i of output
   *  \return SiconosVector
   */
  const siconos::algebra::SiconosVector getLambda(const unsigned int i) const;

  /** get _lambda[i], derivative number i of input
   *
   *  \param i derivative number i of output
   *  \return pointer on a SiconosVector
   */
  inline std::shared_ptr<siconos::algebra::SiconosVector> lambda(const unsigned int i) const {
    assert(_lambda[i]);
    return _lambda[i];
  }

  /** get _lambda[i], derivative number i of input
   *
   *  \param i derivative number i of output
   *  \return pointer on a SiconosVector
   */
  inline siconos::algebra::SiconosVector& lambda_python(const unsigned int i) const {
    assert(_lambda[i]);
    return *(_lambda[i]);
  }

  /** get all the values of the multiplier lambda stored in memory
   *
   *  \param level
   *  \return a memory
   */
  siconos::algebra::SiconosMemory& lambdaMemory(unsigned int level);

  /** get the last value of the multiplier lambda stored in memory
   *
   *  \param level
   *  \return a SiconosVector reference
   */
  const siconos::algebra::SiconosVector& lambda_k(const unsigned int i) const;

  /** set the input vector _lambda to newVector
   *
   *  \param v std::vector<std::shared_ptr<siconos::algebra::SiconosVector>>
   */
  void setLambda(const std::vector<std::shared_ptr<siconos::algebra::SiconosVector>>& v);

  /** set vector _lambda to newVector with direct pointer equality for the _lambda[i]
   *
   *  \param v std::vector<std::shared_ptr<siconos::algebra::SiconosVector>>
   */
  void setLambdaPtr(const std::vector<std::shared_ptr<siconos::algebra::SiconosVector>>& v);

  /** set _lambda[i] to newValue
   *
   *  \param i derivative number i of output
   *  \param newValue a SiconosVector
   */
  void setLambda(const unsigned int i, const siconos::algebra::SiconosVector& newValue);

  /** set _lambda[i] to pointer newPtr
   *
   *  \param i derivative number i of output
   *  \param newPtr a std::shared_ptr<siconos::algebra::SiconosVector>
   */
  void setLambdaPtr(const unsigned int i,
                    std::shared_ptr<siconos::algebra::SiconosVector> newPtr);

  /** get the Relation of this Interaction
   *
   *  \return a pointer on this Relation
   */
  inline std::shared_ptr<Relation> relation() const { return _relation; }

  /** get the NonSmoothLaw of this Interaction
   *
   *  \return a pointer on this NonSmoothLaw
   */
  inline std::shared_ptr<NonSmoothLaw> nonSmoothLaw() const { return _nslaw; }

  inline std::vector<std::shared_ptr<siconos::algebra::BlockVector>>& linkToDSVariables() {
    return _linkToDSVariables;
  };

  inline std::vector<std::shared_ptr<siconos::algebra::SiconosVector>>& relationVectors() {
    return _relationVectors;
  };

  // --- OTHER FUNCTIONS ---

  /** set interaction 'ds-dimension', i.e. sum of all sizes of the dynamical systems linked
   *  by the current interaction. This must be done by topology during call to link(inter, ds,
   * ...).
   *
   *  \param s1 int sum of ds sizes
   */
  inline void setDSSizes(siconos::algebra::Index s1) { _sizeOfDS = s1; }

  /** Must be call to fill the memory. (after convergence of the Newton iterations)
   */
  void swapInMemory();

  /** print the data to the screen
   */
  void display(bool brief = true) const;

  /** reset the global Interaction counter (for ids)
   *
   *  \return the previous value of count
   */
  static inline size_t resetCount(size_t new_count = 0) {
    size_t old_count = count_;
    count_ = new_count;
    return old_count;
  };

  /** Computes output y, depends on the relation type.
   *
   *  \param time current time
   *  \param derivativeNumber number of the derivative to compute, optional, default = 0
   */
  void computeOutput(double time, unsigned int derivativeNumber = 0);

  /** Compute input r of all Dynamical Systems involved in the present
   *  Interaction.
   *
   *  \param time current time
   *  \param level order of _lambda used to compute input.
   */
  void computeInput(double time, unsigned int level = 0);

  /** \return a read-only view on the matrix 'left' that must be used in interactionBlock
   * computation, (left * W * right). It depends on the relation type.
   */
  const siconos::algebra::ConstMapType getLeftInteractionBlock() const;

  /** gets the matrix used in interactionBlock computation

   *  (left * W * right), depends on the relation type (ex, LinearTIR, left = C, right = B).
   *   We get only the part corresponding to one ds.
   *
   *  \param pos int, relative position of the beginning of the required block in relation
   matrix.
   *  \param size int, size(0) of the block
   *  \param sizeDS int, cols() of the block
   *  \return InteractionBlock a pointer to SiconosMatrix (in-out parameter): the resulting
   interactionBlock matrix
   */
  std::shared_ptr<siconos::algebra::SiconosMatrix> getLeftInteractionBlockForDS(
      siconos::algebra::Index pos, siconos::algebra::Index size,
      siconos::algebra::Index sizeDS) const;

  /** gets the matrix used in interactionBlock computation. Used only for the formulation
   * projecting on the constraints. We get only the part corresponding to ds.
   *
   *  \param pos int, relative position of the beginning of the required block in relation
   * matrix. \param InteractionBlock a pointer to SiconosMatrix (in-out parameter): the
   * resulting interactionBlock matrix
   */
  void getLeftInteractionBlockForDSProjectOnConstraints(
      siconos::algebra::Index pos,
      std::shared_ptr<siconos::algebra::SiconosMatrix> InteractionBlock) const;

  /** gets the matrix used in interactionBlock computation
   *
   *  (left * W * rigth), depends on the relation type (ex, LinearTIR, left = C, right = B).
   *  We get only the part corresponding to ds.
   *
   *  \param pos relative position of the beginning of the required block in relation matrix.
   *  \param sizeDS number of rows of the block
   *  \param size number of columns of the block
   *  \return the resulting interactionBlock matrix
   */
  std::shared_ptr<siconos::algebra::SiconosMatrix> getRightInteractionBlockForDS(
      siconos::algebra::Index pos, siconos::algebra::Index sizeDS,
      siconos::algebra::Index size) const;

  /** gets extra interactionBlock corresponding to the present Interaction
   *
   *  \param[in,out] InteractionBlock std::shared_ptr<siconos::algebra::SiconosMatrix>
   */
  void getExtraInteractionBlock(
      std::shared_ptr<siconos::algebra::SiconosMatrix> InteractionBlock) const;
};

}  // namespace siconos::modeling

#endif  // INTERACTION_H
