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
/*! \file StressLinearTIR.hpp

 */
#ifndef STRESSLINEARRELATION_H
#define STRESSLINEARRELATION_H

#include "LagrangianLinearTIR.hpp"
namespace siconos::modeling {
/**
    Stress Linear Relation.

    Stress Relation with:

    \f$ y= [0 C] [v ; \sigma] + e + Fz \f$

    \f$ p = [0 C]^t \lambda = [0 C^t \lambda] \f$

    C will typically be identity for a StressLinearTIR. \lambda will be intepreted as the plastic rate.

*/
class StressLinearTIR : public LagrangianLinearTIR {
 protected:
  ACCEPT_SERIALIZATION(StressLinearTIR);

 public:
  enum SolidLinearDS { sigma, dotEpsilon, solidDSlinkSize };


  /** Default constructor
   */
  StressLinearTIR() : LagrangianLinearTIR(){};

  /** create the Relation from a set of data
   *
   *  \param C the matrix C
   */
  StressLinearTIR(std::shared_ptr<siconos::algebra::SimpleMatrix> C);

  /** create the Relation from a set of data
   *
   *  \param C the matrix C
   *  \param F the matrix F
   *  \param e the vector e
   */
  StressLinearTIR(std::shared_ptr<siconos::algebra::SimpleMatrix> C,
                      std::shared_ptr<siconos::algebra::SimpleMatrix> F,
                      std::shared_ptr<siconos::algebra::SiconosVector> e);

  /** create the Relation from a set of data
   *
   *  \param C the matrix C
   *  \param e the vector e
   */
  StressLinearTIR(std::shared_ptr<siconos::algebra::SimpleMatrix> C,
                      std::shared_ptr<siconos::algebra::SiconosVector> e);

  /** destructor
   */
  virtual ~StressLinearTIR() noexcept = default;

  /** check sizes of the relation specific operators.
   *
   *  \param inter an Interaction using this relation
   */
  void checkSize(Interaction &inter) override;
  ;

  /** default function to compute y
   *
   *  \param time not used
   *  \param inter the Interaction we want to update
   *  \param derivativeNumber the derivative of y we want to compute
   */
  void computeOutput(double time, Interaction &inter,
                     unsigned int derivativeNumber = 0) override;
  ;

  /** default function to compute r
   *
   *  \param time not used
   *  \param inter the Interaction we want to update
   *  \param level the derivative of lambda we want to compute
   */
  void computeInput(double time, Interaction &inter, unsigned int level = 0) override;

  /** compute all the H Jacobian
   *
   *  \param time not used
   *  \param inter the Interaction we want to update
   *  \param interProp interaction properties
   */
  void computeJach(double time, Interaction &inter) override {}

  /** compute all the G Jacobian
   *
   *  \param time not used
   *  \param inter the Interaction we want to update
   *  \param interProp interaction properties
   */
  void computeJacg(double time, Interaction &inter) override {}

  // GETTERS/SETTERS

  // -- C --
  /** \return pointer on a plugged matrix
   */
  inline std::shared_ptr<siconos::algebra::SimpleMatrix> C() const override { return _jachq; }

  /** set C to pointer newPtr
   *
   *  \param newPtr a SP to plugged matrix
   */
  inline void setCPtr(std::shared_ptr<siconos::algebra::SimpleMatrix> newPtr)
  {
    _jachq = newPtr;
  }

  // -- D --

  /** \return pointer on a plugged matrix
   */
  inline std::shared_ptr<siconos::algebra::SimpleMatrix> D() const { return _jachlambda; }

  /** set D to pointer newPtr
   *
   *  \param newPtr a SP to plugged matrix
   */
  inline void setDPtr(std::shared_ptr<siconos::algebra::SimpleMatrix> newPtr)
  {
    _jachlambda = newPtr;
  }

  // -- F --

  /** \return pointer on a plugged matrix
   */
  inline std::shared_ptr<siconos::algebra::SimpleMatrix> F() const { return _F; }

  /** set F to pointer newPtr
   *
   *  \param newPtr a SP to plugged matrix
   */
  inline void setFPtr(std::shared_ptr<siconos::algebra::SimpleMatrix> newPtr) { _F = newPtr; }

  // -- e --

  /** \return pointer on a plugged vector
   */
  inline std::shared_ptr<siconos::algebra::SiconosVector> e() const { return _e; }

  /** set e to pointer newPtr
   *
   *  \param newPtr a SP to plugged vector
   */
  inline void setEPtr(std::shared_ptr<siconos::algebra::SiconosVector> newPtr) { _e = newPtr; }

  /** print the data to the screen
   */
  void display() const override;

  /** \return true if the relation is linear.
   */

  bool isLinear() override  // final would be better but swig does not like it
  {
    return true;
  }
};

}  // namespace siconos::modeling
#endif
