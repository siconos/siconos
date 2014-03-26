/* Siconos-Kernel, Copyright INRIA 2005-2012.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
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
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Apr 26, 2004
 *
 *
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

  /* default constructor */
  QP() {};

public:

  /** Destructor */
  ~QP();

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
   *  \param SiconosMatrix newValue
   */
  inline void setQ(const SiconosMatrix& newValue)
  {
    *_Q = newValue;
  }

  /** set Q to pointer newPtr
   *  \param SP::SiconosMatrix  newPtr
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
   *  \param SiconosVector newValue
   */
  inline void setP(const SiconosVector& newValue)
  {
    *_p = newValue;
  }

  /** set p to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  inline void setPPtr(SP::SiconosVector newPtr)
  {
    _p = newPtr;
  }

  // --- OTHER FUNCTIONS ---

  /** To run the solver for ns problem
   *   \param double : current time
   *  \return int, information about the solver convergence.
   */
  int compute(double);

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
