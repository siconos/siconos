// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file LagrangianCompliantLinearTIR.hpp
 */
#ifndef LAGRANGIANCOMPLIANTLINEARRELATION_H
#define LAGRANGIANCOMPLIANTLINEARRELATION_H

#include "LagrangianR.hpp"

namespace siconos::modeling {
/**  Lagrangian Compliant Linear Relation.

This class provides tools to describe Lagrangian (2nd order) linear relations, with:

\f$
y = H q + e + D\lambda
\f$

\f$
p = H^\top \lambda
\f$

Minimal setup requires the H (called C in the constructors) and D matrices.

With the notations of LagrangianR and LagrangianCompliantR,
\f$ H = \nabla_q h(q,\lambda) \f$ (contactMatrix()) and \f$ D = \nabla_{\lambda}h(q,\lambda)
\f$ (D()). CMatrix() and DMatrix() are the former (deprecated) names of contactMatrix() and
D().

*/
class LagrangianCompliantLinearTIR : public LagrangianR {
 protected:
  ACCEPT_SERIALIZATION(LagrangianCompliantLinearTIR);

  /** D matrix */
  std::shared_ptr<siconos::algebra::MapType> D_view_{nullptr};

  /** e*/
  std::shared_ptr<siconos::algebra::MapVectorType> eVector_view_{nullptr};

 public:
  /** create the Relation from a set of data
   *  \param C the matrix C
   *  \param D the matrix D
   */
  LagrangianCompliantLinearTIR(Eigen::Ref<siconos::algebra::SiconosMatrix> C,
                               Eigen::Ref<siconos::algebra::SiconosMatrix> D);

  /** create a complete Relation from a set of data
   *  \param C matrix operator (constant) C
   *  \param D matrix operator (constant) D
   *  \param e vector operator (constant) e
   */
  LagrangianCompliantLinearTIR(Eigen::Ref<siconos::algebra::SiconosMatrix> C,
                               Eigen::Ref<siconos::algebra::SiconosMatrix> D,
                               Eigen::Ref<siconos::algebra::SiconosVector> e);

  /** destructor */
  virtual ~LagrangianCompliantLinearTIR() noexcept = default;

  /** initialize LagrangianCompliantLinearTIR specific operators.
   * \param inter an Interaction using this relation
   */
  void initialize(Interaction& inter) override;

  /** check sizes LagrangianCompliantLinearTIR specific operators.
   * \param inter an Interaction using this relation
   */
  void checkSize(const Interaction& inter) const override;

  /** \deprecated use contactMatrix() */
  [[deprecated("use contactMatrix()")]] inline auto CMatrix() const { return contactMatrix(); }

  /** \return a read-only view on \f$ D = \nabla^\top_{\lambda}h(q,\lambda) \f$ */
  inline auto D() const {
    return siconos::algebra::ConstMapType(D_view_->data(), D_view_->rows(), D_view_->cols());
  }

  /** \deprecated use D() */
  [[deprecated("use D()")]] inline auto DMatrix() const { return D(); }

  /** \return a read-only view on \f$ D \f$ (the Relation interface, same as D()) */
  inline const siconos::algebra::ConstMapType jacobianhOver_lambda() const override {
    return D();
  }

  /** \return  a read-only view on e vector */
  inline auto eVector() const {
    return siconos::algebra::ConstMapVectorType(eVector_view_->data(), eVector_view_->size());
  }

  virtual bool hasJacobianhOver_lambda() const override { return D_view_ != nullptr; }

  /** True if e is defined */
  bool haseVector() const { return eVector_view_ != nullptr; }

  /** default function to compute y
   *  \param time dummy parameter for this kind of relation
   *  \param inter the Interaction we want to update
   *  \param derivativeNumber the derivative of y we want to compute
   */
  void computeOutput(double time, Interaction& inter,
                     siconos::algebra::blocks::size_type derivativeNumber = 0) override;

  /** default function to compute r
   *  \param time dummy parameter for this kind of relation
   *  \param inter the Interaction we want to update
   *  \param level the derivative of lambda we want to compute
   */
  void computeInput(double time, Interaction& inter,
                    siconos::algebra::blocks::size_type level = 0) override;

  /* compute all the H Jacobian
   *  \param time dummy parameter for this kind of relation
   *  \param inter the Interaction we want to update
   *  \param interProp interaction properties
   */
  void computeJach(double time, Interaction& inter) override {}

  /** print the data to the screen
   */
  void display() const override;

  /**
   * \return true if the relation is linear.
   */

  bool isLinear() const override { return true; }
};
}  // namespace siconos::modeling
#endif
