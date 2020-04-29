/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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
#include "SimpleMatrix.hpp"
#include "SiconosVector.hpp"


/** Quadratic Problem
 *
 */
class QP : public OneStepNSProblem
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(QP);


  /** contains the Q matrix of a QP problem */
  SP::SiconosMatrix _Q;

  /** contains the p vector of a QP problem */
  SP::SiconosVector _p;

  //  /** contains the data of the QP, according to siconos/numerics */
  //  QPStructure QPMethod;

public:

  /* default constructor */
  QP() = default;

  /** Destructor */
  ~QP(){};

  // --- GETTERS/SETTERS ---

  // --- Q ---
  /** get the value of Q
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getQ() const
  {
    return *_Q;
  }

  /** get Q
   *  \return pointer on a SiconosMatrix
   */
  inline SP::SiconosMatrix q() const
  {
    return _Q;
  }

  /** set the value of Q to newValue
   *  \param newValue SiconosMatrix 
   */
  inline void setQ(const SiconosMatrix& newValue)
  {
    *_Q = newValue;
  }

  /** set Q to pointer newPtr
   *  \param newPtr the new matrix
   */
  inline void setQPtr(SP::SiconosMatrix newPtr)
  {
    _Q = newPtr;
  }

  // --- P ---
  /** get the value of p, the initial state of the DynamicalSystem
   *  \return SiconosVector
   *  \warning: SiconosVector is an abstract class => can not be an lvalue => return SiconosVector
   */
  inline const SiconosVector getP() const
  {
    return *_p;
  }

  /** get p, the initial state of the DynamicalSystem
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector p() const
  {
    return _p;
  }

  /** set the value of p to newValue
   *  \param newValue SiconosVector 
   */
  inline void setP(const SiconosVector& newValue)
  {
    *_p = newValue;
  }

  /** set p to pointer newPtr
   *  \param newPtr SiconosVector * 
   */
  inline void setPPtr(SP::SiconosVector newPtr)
  {
    _p = newPtr;
  }

  // --- OTHER FUNCTIONS ---

  /** To run the solver for ns problem
   *   \param time current time
   *  \return int, information about the solver convergence.
   */
  int compute(double time);

  /** print the data to the screen
   */
  void display() const;

  /* pure virtual in OneStepNSProblem.hpp */
  void computeInteractionBlock(const InteractionsGraph::EDescriptor&)
  {
    assert(false);
  }

  void computeDiagonalInteractionBlock(const InteractionsGraph::VDescriptor&)
  {
    assert(false);
  }

  virtual bool preCompute(double time)
  {
    assert(false);
    return false;
  }

  void postCompute()
  {
    assert(false);
  }

};

#endif // QP_H
