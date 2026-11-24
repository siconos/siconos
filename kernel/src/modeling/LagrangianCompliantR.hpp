// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file LagrangianCompliantR.hpp

 */
#ifndef LagrangianCompliantR_H
#define LagrangianCompliantR_H

#include "LagrangianR.hpp"

namespace siconos::modeling {
/**
    @brief Lagrangian (2nd order) compliant relations, with \f$ y = h(q, \lambda) \f$
    (dense storage)

    \f[
    \begin{aligned}
      y &= h(q, \lambda), \\
      \dot y &= H(q, \lambda)\, \dot q + D(q, \lambda)\, \dot\lambda, \\
      p &= H(q, \lambda)^\top \lambda,
    \end{aligned}
    \qquad
    H = \nabla_q^\top h(q, \lambda), \quad D = \nabla_\lambda^\top h(q, \lambda).
    \f]

    \f$ H \f$ is the contact configuration matrix (see LagrangianR), and \f$ D \f$ has the
    same meaning as in the first-order relations. The following operators can be set by
    user-defined functions, or computed by overriding the virtual functions:

    - \f$ h(q, \lambda) \f$: setComputehFunction, or override computeh;
    - \f$ H(q, \lambda) \f$: setComputeContactMatrixFunction, or override computeContactMatrix;
    - \f$ D(q, \lambda) \f$: setComputeDFunction, or override computeD.

    <b>Former names</b> (deprecated, kept for one release): setComputeJacobianhOver_qFunction,
    computeJacobianhOver_q, setComputeJacobianhOver_lambdaFunction,
    computeJacobianhOver_lambda (a class overriding them is still called).
    jacobianhOver_lambda() (the Relation interface) returns D().
*/
class LagrangianCompliantR : public LagrangianR {
 protected:
  /** function wrapper used to compute \f$ h(q, \lambda) \f$ */
  siconos::modeling::func_prototypes::FunctionBVV_V computeh_{nullptr};

  /** function wrapper used to compute \f$ H = \nabla^\top_q h(q, \lambda) \f$ */
  siconos::modeling::func_prototypes::FunctionBVV_M computeContactMatrix_{nullptr};

  /** \f$ D = \nabla^\top_{\lambda}h(q, \lambda) \f$  */
  std::shared_ptr<siconos::algebra::SiconosMatrix> D_{nullptr};

  /** function wrapper used to compute \f$ D = \nabla^\top_{\lambda} h(q, \lambda) \f$ */
  siconos::modeling::func_prototypes::FunctionBVV_M computeD_{nullptr};

  ACCEPT_SERIALIZATION(LagrangianCompliantR);

  /** initialize H and D matrices or components specific to derived classes
   *
   *  \param inter : the Interaction
   */
  void initialize(Interaction& inter) override;

 private:
  /** flags of the dispatch to the deprecated virtual names (see RelationCompatibility.hpp) */
  bool computeContactMatrixLegacy_{false};
  bool computeDLegacy_{false};

 public:
  /** default and only constructor */
  LagrangianCompliantR() : LagrangianR(RelationSubType::CompliantR) {};

  /** destructor */
  virtual ~LagrangianCompliantR() noexcept = default;

  /** set a user-defined function to compute \f$ h(q, \lambda) \f$
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputehFunction(const siconos::modeling::func_prototypes::FunctionBVV_V& fct);

  /**
    to compute the output y = h(q, \lambda) of the Relation

    \param q coordinates of the dynamical systems involved in the relation
    \param lambda interaction  \f$ \lambda \f$  vector
    \param y the resulting vector
  */
  virtual void computeh(const siconos::algebra::BlockVector& q,
                        const Eigen::Ref<const siconos::algebra::SiconosVector>& lambda,
                        Eigen::Ref<siconos::algebra::SiconosVector> y);

  /** set a user-defined function to compute \f$ H = \nabla^\top_q h(q, \lambda) \f$
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeContactMatrixFunction(
      const siconos::modeling::func_prototypes::FunctionBVV_M& fct);

  /** Computes \f$ H = \nabla^\top_q h(q, \lambda) \f$. The default calls the user-defined
   *  function, if any. Override this function to compute \f$ H \f$ in a derived class.
   *
   *  \param q coordinates of the dynamical systems involved in the relation
   *  \param lambda interaction  \f$ \lambda \f$  vector
   */
  virtual void computeContactMatrix(
      const siconos::algebra::BlockVector& q,
      const Eigen::Ref<const siconos::algebra::SiconosVector>& lambda);

  /** set a user-defined function to compute \f$ D = \nabla^\top_{\lambda} h(q, \lambda) \f$
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeDFunction(const siconos::modeling::func_prototypes::FunctionBVV_M& fct);

  /** Computes \f$ D = \nabla^\top_{\lambda} h(q, \lambda) \f$. The default calls the
   *  user-defined function, if any. Override this function to compute \f$ D \f$.
   *
   *  \param q coordinates of the dynamical systems involved in the relation
   *  \param lambda interaction  \f$ \lambda \f$  vector
   */
  virtual void computeD(const siconos::algebra::BlockVector& q,
                        const Eigen::Ref<const siconos::algebra::SiconosVector>& lambda);

  /** \return a read-only view on \f$ D = \nabla^\top_{\lambda}h(q,\lambda) \f$ */
  inline auto D() const {
    return siconos::algebra::ConstMapType(D_->data(), D_->rows(), D_->cols());
  }

  /** \return a read-only view on \f$ D = \nabla^\top_{\lambda}h(q,\lambda) \f$ (the
   *  Relation interface, same as D()) */
  inline const siconos::algebra::ConstMapType jacobianhOver_lambda() const override {
    return D();
  }

  /** to compute output
   *
   *  \param time the current time
   *  \param inter the Interaction owning y
   *  \param derivativeNumber the number of the derivative to compute, optional, default = 0.
   */
  void computeOutput(double time, Interaction& inter,
                     siconos::algebra::blocks::size_type derivativeNumber = 0) override;
  /** to compute the input
   *
   *  \param time the current time
   *  \param inter the Interaction owning lambda
   *  \param level "derivative" order of lambda used to compute input
   */
  void computeInput(double time, Interaction& inter,
                    siconos::algebra::blocks::size_type level = 0) override;

  /** compute H and D */
  void computeJach(double time, Interaction& inter) override;

  /** \name Deprecated names, kept for one release */
  ///@{
  [[deprecated("use setComputeContactMatrixFunction")]] void setComputeJacobianhOver_qFunction(
      const siconos::modeling::func_prototypes::FunctionBVV_M& fct) {
    setComputeContactMatrixFunction(fct);
  }
  /** \deprecated former name of computeContactMatrix. Still called for classes that override
   * it. */
  [[deprecated("override/call computeContactMatrix")]] virtual void computeJacobianhOver_q(
      const siconos::algebra::BlockVector& q,
      const Eigen::Ref<const siconos::algebra::SiconosVector>& lambda);
  [[deprecated("use setComputeDFunction")]] void setComputeJacobianhOver_lambdaFunction(
      const siconos::modeling::func_prototypes::FunctionBVV_M& fct) {
    setComputeDFunction(fct);
  }
  /** \deprecated former name of computeD. Still called for classes that override it. */
  [[deprecated("override/call computeD")]] virtual void computeJacobianhOver_lambda(
      const siconos::algebra::BlockVector& q,
      const Eigen::Ref<const siconos::algebra::SiconosVector>& lambda);
  ///@}
};
}  // namespace siconos::modeling

#endif
