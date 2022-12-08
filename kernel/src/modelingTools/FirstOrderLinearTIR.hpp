/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

/** 
    Linear Time Invariant Relation, derived from class FirstOrderR

    Linear Relation for First Order Dynamical Systems with time-independant
    operators
    
    \f[
    y &=& Cx(t) + Fz + D\lambda + e \\      
    R &=& B\lambda 
    \f]
    
*/
class FirstOrderLinearTIR : public FirstOrderR {

protected:

  ACCEPT_SERIALIZATION(FirstOrderLinearTIR);

  /** e operator (constant vector) */
  SP::SiconosVector _e{nullptr};

  /** initialize the relation (check sizes, memory allocation ...)
   *
   *  \param inter the interaction that owns this relation
   */
  void initialize(Interaction &inter) override;

  /** check sizes of the relation specific operators.
   *
   *  \param inter an Interaction using this relation
   */
  void checkSize(Interaction &inter) override;

public:
  /** create the Relation from a set of data
   *
   *  \param C the matrix C
   *  \param B the matrix B
   */
  FirstOrderLinearTIR(SP::SimpleMatrix C, SP::SimpleMatrix B);

  /** create the Relation from a set of data
   *
   *  \param C the C matrix
   *  \param D the D matrix
   *  \param F the F matrix
   *  \param e the e matrix
   *  \param B the B matrix
   */
  FirstOrderLinearTIR(SP::SimpleMatrix C, SP::SimpleMatrix D,
                      SP::SimpleMatrix F, SP::SiconosVector e,
                      SP::SimpleMatrix B);

  /** destructor
   */
  virtual ~FirstOrderLinearTIR() noexcept = default;

  /** default function to compute h = y = Cx(t) + Fz + Dlambda + e
   *
   *  \param x
   *  \param lambda
   *  \param z
   *  \param y the resulting vector
   */
  void computeh(const BlockVector &x, const SiconosVector &lambda, BlockVector &z,
                SiconosVector &y);

  /** default function to compute g = Blambda
   *
   *  \param lambda
   *  \param r non-smooth input
   */
  void computeg(const SiconosVector &lambda, BlockVector &r);

  /** default function to compute y
   *
   *  \param time current time
   *  \param inter Interaction using this Relation
   *  \param level
   */
  void computeOutput(double time, Interaction &inter,
                     unsigned int level = 0) override;

  /** default function to compute r
   *
   *  \param time current time
   *  \param inter Interaction using this Relation
   *  \param level
   */
  void computeInput(double time, Interaction &inter,
                    unsigned int level = 0) override;

  /** print the data to the screen
   */
  void display() const override;

  /** set e
   *
   *  \param  newe the new value of e
   */
  inline void setePtr(SP::SiconosVector newe) { _e = newe; }

  /** get e
   *
   *  \return e matrix
   */
  inline SP::SiconosVector e() const { return _e; }

  /** determine if the Relation is linear
   *
   *  \return true if the relation is linear.
   */
  inline bool isLinear() override { return true; }

  // Jacobians: required to fullfill base abstract class API but do nothing.
  // Note FP: final would be better than override but swig cannot handle it.
  void computeJach(double time, Interaction &inter) override{};
  void computeJacg(double time, Interaction &inter) override{};

  ACCEPT_STD_VISITORS();
};

TYPEDEF_SPTR(FirstOrderLinearTIR)

#endif
