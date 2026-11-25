// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file LagrangianLinearTIR.hpp

 */
#ifndef LAGRANGIANLINEARRELATION_H
#define LAGRANGIANLINEARRELATION_H

#include "LagrangianR.hpp"

namespace siconos::modeling {
/**
    Lagrangian Linear Relation.

    Lagrangian Relation with:

    \f$ y = H q + e \f$

    \f$ p = H^\top \lambda \f$

    with a constant contact configuration matrix \f$ H \f$ (contactMatrix(), called C in the
    constructors and formerly returned by CMatrix(), deprecated). \f$ H \f$ is the only
    required input to build a LagrangianLinearTIR.

*/
class LagrangianLinearTIR : public LagrangianR {
 protected:
  ACCEPT_SERIALIZATION(LagrangianLinearTIR);

  /** e*/
  std::shared_ptr<siconos::algebra::MapVectorType> eVector_view_{nullptr};

 public:
  /** Build a relation with only C matrix defined (i.e. e=0)
   *
   *  \param C matrix operator (constant) C
   */
  LagrangianLinearTIR(Eigen::Ref<siconos::algebra::SiconosMatrix> C);

  /** Build a complete relation (C and e). Warning: shared memory between input and class
   * attributes.
   *
   *  \param C matrix operator (constant) C
   *  \param e vector operator (constant) e
   */
  LagrangianLinearTIR(Eigen::Ref<siconos::algebra::SiconosMatrix> C,
                      Eigen::Ref<siconos::algebra::SiconosVector> e);

  /** destructor */
  virtual ~LagrangianLinearTIR() noexcept = default;

  /** check sizes of the relation specific operators.
   *
   *  \param inter an Interaction using this relation
   */
  void checkSize(const Interaction& inter) const override;

  /** \deprecated use contactMatrix() */
  [[deprecated("use contactMatrix()")]] inline auto CMatrix() const { return contactMatrix(); }

  /** \return  a read-only view on e vector */
  inline auto eVector() const {
    return siconos::algebra::ConstMapVectorType(eVector_view_->data(), eVector_view_->size());
  }

  /** True if e is defined */
  bool haseVector() const { return eVector_view_ != nullptr; }

  /** default function to compute y
   *
   *  \param time not used
   *  \param inter the Interaction we want to update
   *  \param derivativeNumber the derivative of y we want to compute
   */
  void computeOutput(double time, Interaction& inter,
                     siconos::algebra::blocks::size_type derivativeNumber = 0) override;

  /** default function to compute r
   *
   *  \param time not used
   *  \param inter the Interaction we want to update
   *  \param level the derivative of lambda we want to compute
   */
  void computeInput(double time, Interaction& inter,
                    siconos::algebra::blocks::size_type level = 0) override;

  /** compute all the H Jacobian
   *
   *  \param time not used
   *  \param inter the Interaction we want to update
   *  \param interProp interaction properties
   */
  void computeJach(double time, Interaction& inter) override {}

  /** print the data to the screen
   */
  void display() const override;

  /** \return true if the relation is linear.
   */

  bool isLinear() const final { return true; }
};

}  // namespace siconos::modeling
#endif
