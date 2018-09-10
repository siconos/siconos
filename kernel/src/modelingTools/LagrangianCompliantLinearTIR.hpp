/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
/*! \file LagrangianCompliantLinearTIR.hpp
 */
#ifndef LAGRANGIANCOMPLIANTLINEARRELATION_H
#define LAGRANGIANCOMPLIANTLINEARRELATION_H

#include "LagrangianR.hpp"


/**  Lagrangian Linear Relation.


Lagrangian Relation with:

\f$
y= Cq + e + D\lambda + Fz
\f$

\f$
p = C^t \lambda
\f$

C is the only required input to built a LagrangianCompliantLinearTIR.
D is mandatory and may represent a stiffness in the relation




 */
class LagrangianCompliantLinearTIR : public LagrangianR
{

protected:
  /* serialization hooks */
  ACCEPT_SERIALIZATION(LagrangianCompliantLinearTIR);

  /** F matrix, coefficient of z */
  SP::SimpleMatrix _F;

  /** e*/
  SP::SiconosVector _e;

public:

  /** Default constructor
  */
  LagrangianCompliantLinearTIR() : LagrangianR(RELATION::CompliantLinearTIR) {};

  /** create the Relation from a set of data
  *  \param C the matrix C
  *  \param D the matrix D
  */
  LagrangianCompliantLinearTIR(SP::SimpleMatrix C, SP::SimpleMatrix D);

  /** create the Relation from a set of data
  *  \param C the matrix C
  *  \param D the matrix D
  *  \param e the vector e
  */
  LagrangianCompliantLinearTIR(SP::SimpleMatrix C, SP::SimpleMatrix D, SP::SiconosVector e);

  /** create the Relation from a set of data
  *  \param C the matrix C
  *  \param D the matrix D
  *  \param F the matrix F
  *  \param e the vector e
  */
  LagrangianCompliantLinearTIR(SP::SimpleMatrix C, SP::SimpleMatrix D, SP::SimpleMatrix F, SP::SiconosVector e);

  /** destructor
  */
  virtual ~LagrangianCompliantLinearTIR() {};

  /** initialize LagrangianCompliantLinearTIR specific operators.
   * \param inter an Interaction using this relation
   */
  void initialize(Interaction& inter);

  /** check sizes LagrangianCompliantLinearTIR specific operators.
   * \param inter an Interaction using this relation
   */
  void checkSize(Interaction& inter);

  /** default function to compute y
  *  \param time not used
  *  \param inter the Interaction we want to update
  *  \param derivativeNumber the derivative of y we want to compute
  */
  void computeOutput(double time, Interaction& inter, unsigned int derivativeNumber = 0);

  /** default function to compute r
  *  \param time not used
  *  \param inter the Interaction we want to update
  *  \param level the derivative of lambda we want to compute
  */
  void computeInput(double time, Interaction& inter, unsigned int level = 0);
  /* compute all the H Jacobian
   *  \param time not used
   *  \param inter the Interaction we want to update
   *  \param interProp interaction properties
   */
  void computeJach(double time, Interaction& inter)
  {
    ;
  }


  /* compute all the G Jacobian
   *  \param time not used
   *  \param inter the Interaction we want to update
   *  \param interProp interaction properties
   */
  void computeJacg(double time, Interaction& inter)
  {
    ;
  }





  // GETTERS/SETTERS

  // -- C --
  /** get C
   *  \return pointer on a plugged matrix
   */
  inline SP::SimpleMatrix C() const
  {
    return _jachq;
  }


  /** set C to pointer newPtr
   *  \param newPtr a SP to plugged matrix
   */
  inline void setCPtr(SP::SimpleMatrix newPtr)
  {
    _jachq = newPtr;
  }

  // -- D --

  /** get D
   *  \return pointer on a plugged matrix
   */
  inline SP::SimpleMatrix D() const
  {
    return _jachlambda;
  }


  /** set D to pointer newPtr
   * \param newPtr a SP to plugged matrix
   */
  inline void setDPtr(SP::SimpleMatrix newPtr)
  {
    _jachlambda = newPtr;
  }

  // -- F --

  /** get F
  *  \return pointer on a plugged matrix
  */
  inline SP::SimpleMatrix F() const
  {
    return _F;
  }

  /** set F to pointer newPtr
   * \param newPtr a SP to plugged matrix
   */
  inline void setFPtr(SP::SimpleMatrix newPtr)
  {
    _F = newPtr;
  }

  // -- e --

  /** get e
   *  \return pointer on a plugged vector
   */
  inline SP::SiconosVector e() const
  {
    return _e;
  }

  /** set e to pointer newPtr
   *  \param newPtr a SP to plugged vector
   */
  inline void setEPtr(SP::SiconosVector newPtr)
  {
    _e = newPtr;
  }

  /** get a pointer on matrix Jach[index]
   *  \return a pointer on a SimpleMatrix
   */

  /** print the data to the screen
   */
  void display() const;

  /**
   * \return true if the relation is linear.
   */

  virtual bool isLinear()
  {
    return true;
  }
  ACCEPT_STD_VISITORS();

};

TYPEDEF_SPTR(LagrangianCompliantLinearTIR)

#endif
