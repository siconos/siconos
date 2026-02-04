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
/*! \file
 */
#ifndef QP_H
#define QP_H

#include "OneStepNSProblem.hpp"

namespace siconos::nonsmooth_formulations {

/** Quadratic Problem
 *
 */
class QP : public OneStepNSProblem {
 private:
  ACCEPT_SERIALIZATION(QP);

  /** contains the Q matrix of a QP problem */
  std::shared_ptr<siconos::algebra::SiconosMatrix> _Q{nullptr};

  /** contains the p vector of a QP problem */
  std::shared_ptr<siconos::algebra::SiconosVector> _p{nullptr};

  //  /** contains the data of the QP, according to siconos/numerics */
  //  QPStructure QPMethod;

 public:
  /* default constructor */
  QP() = default;

  /** Destructor */
  ~QP() noexcept = default;

  /** get Q
   *
   *  \return pointer on a SiconosMatrix
   */
  inline std::shared_ptr<siconos::algebra::SiconosMatrix> q() const { return _Q; }

  /** set the value of Q to newValue
   *
   *  \param newValue SiconosMatrix
   */
  void setQ(const siconos::algebra::SiconosMatrix& newValue);

  /** set Q to pointer newPtr
   *
   *  \param newPtr the new matrix
   */
  inline void setQPtr(std::shared_ptr<siconos::algebra::SiconosMatrix> newPtr) { _Q = newPtr; }

  /** get p, the initial state of the DynamicalSystem
   *
   *  \return pointer on a siconos::algebra::SiconosVector
   */
  inline std::shared_ptr<siconos::algebra::SiconosVector> p() const { return _p; }

  /** set the value of p to newValue
   *
   *  \param newValue siconos::algebra::SiconosVector
   */
  void setP(const siconos::algebra::SiconosVector& newValue);

  /** set p to pointer newPtr
   *
   *  \param newPtr siconos::algebra::SiconosVector *
   */
  inline void setPPtr(std::shared_ptr<siconos::algebra::SiconosVector> newPtr) { _p = newPtr; }

  // --- OTHER FUNCTIONS ---

  /** To run the solver for ns problem
   *
   *  \param time current time
   *  \return int, information about the solver convergence.
   */
  int compute(double time);

  /** print the data to the screen
   */
  void display() const;

  /* pure virtual in OneStepNSProblem.hpp */
  void computeInteractionBlock(const siconos::graphs::InteractionsGraph::EDescriptor&)
  {
    assert(false);
  }

  void computeDiagonalInteractionBlock(const siconos::graphs::InteractionsGraph::VDescriptor&)
  {
    assert(false);
  }

  virtual bool preCompute(double time)
  {
    assert(false);
    return false;
  }

  void postCompute() { assert(false); }
};
}  // namespace siconos::nonsmooth_formulations

#endif  // QP_H
