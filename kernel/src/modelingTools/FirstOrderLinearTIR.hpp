/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
/*! \file FirstOrderLinearTIR.hpp

 */
#ifndef FirstOrderLinearTIR_H
#define FirstOrderLinearTIR_H

#include "FirstOrderR.hpp"

/** Linear Time Invariant Relation, derived from class FirstOrderR

\author SICONOS Development Team - copyright INRIA
\version 3.0.0.
\date Apr 15, 2007

Linear Relation for First Order Dynamical Systems:

\f{eqnarray}
y &=& Cx(t) + Fz + D\lambda + e\\

R &=& B\lambda
\f}

 */
class FirstOrderLinearTIR : public FirstOrderR
{

protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(FirstOrderLinearTIR);

  /** initialize the relation (check sizes, memory allocation ...)
   * \param inter the interaction that owns this relation
   * \param DSlink
   * \param workV
   * \param workM
   */
  virtual void initComponents(Interaction& inter, VectorOfBlockVectors& DSlink,
                              VectorOfVectors& workV, VectorOfSMatrices& workM);

  SP::SiconosVector _e;

public:

  /** default constructor, protected
  */
  FirstOrderLinearTIR();

  /** create the Relation from a set of data
  *  \param C the matrix C
  *  \param B the matrix B
  */
  FirstOrderLinearTIR(SP::SimpleMatrix C, SP::SimpleMatrix B);

  /** create the Relation from a set of data
   *  \param C the C matrix
   *  \param D the D matrix
   *  \param F the F matrix
   *  \param e the e matrix
   *  \param B the B matrix
   */
  FirstOrderLinearTIR(SP::SimpleMatrix C, SP::SimpleMatrix D, SP::SimpleMatrix F, SP::SiconosVector e, SP::SimpleMatrix B);

  /** destructor
   */
  virtual ~FirstOrderLinearTIR() {};

  // GETTERS/SETTERS

  /** default function to compute h
   *  \param x
   *  \param lambda
   *  \param z
   *  \param y value of h
   */
  void computeh(BlockVector& x, SiconosVector& lambda, BlockVector& z, SiconosVector& y);

  /** default function to compute g
   *  \param lambda
   *  \param r non-smooth input
   */
  void computeg(SiconosVector& lambda, BlockVector& r);

  /** default function to compute y
   *  \param time current time
   *  \param inter Interaction using this Relation
   *  \param interProp
   *  \param level not used
   */
  virtual void computeOutput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int level = 0);

  /** default function to compute r
   *  \param time current time
   *  \param inter Interaction using this Relation
   *  \param interProp
   *  \param level not used
   */
  virtual void computeInput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int level = 0);

  /** print the data to the screen
   */
  void display() const;

  /** compute the jacobian of h: nothing to be done here
   *  \param time current time
   *  \param inter Interaction using this Relation
   *  \param interProp
   */
  virtual void computeJach(double time, Interaction& inter, InteractionProperties& interProp) {};

  /** compute the jacobian of g: nothing to be done here
   *  \param time current time
   *  \param inter Interaction using this Relation
   *  \param interProp
   */
  virtual void computeJacg(double time, Interaction& inter, InteractionProperties& interProp) {};
 

  /** set e
   *  \param  newe the new value of e
   */
  inline void setePtr(SP::SiconosVector newe)
  {
    _e = newe;
  }

  /** get e
   *  \return e matrix
   */
  inline SP::SiconosVector e() const
  {
    return _e;
  }

  /**
   * \return true if the relation is linear.
   */
  virtual bool isLinear()
  {
    return true;
  }

  ACCEPT_STD_VISITORS();

};

TYPEDEF_SPTR(FirstOrderLinearTIR)

#endif
