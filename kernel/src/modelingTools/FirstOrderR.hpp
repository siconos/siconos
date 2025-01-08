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

/*! \file FirstOrderR.hpp
  General interface for relations.
*/

#ifndef FirstOrderR_H
#define FirstOrderR_H

#include "Relation.hpp"

namespace siconos::modeling {
/**
   FirstOrder Relation

   This is an abstract class for all relation operating on first order systems.
   The following subclasses can be used:

   - FirstOrderNonlinearR: for fully nonlinear relations: \f$ y = h(X,t,\lambda) \f$ , \f$ R =
   g(X,t, lambda) \f$ .
   - FirstOrderType2R: specialization with \f$ y = h(X, \lambda) \f$ , \f$ R
   = g(\lambda) \f$ .
   - FirstOrderType1R: further specialization with \f$ y = h(X) \f$ , \f$ R
   = g(\lambda) \f$ .
   - FirstOrderLinearR: linear case: \f$ y = C(t)X + D(t)\lambda + e(t) \f$ , \f$ R =
   B(t)\lambda \f$ .
   - FirstOrderLinearTIR: time-invariant linear case: \f$ y = CX + D\lambda+ e \f$ , \f$ R =
   B\lambda \f$ .

   If the relation involves only one DynamicalSystem, then  \f$ R = r \f$ ,  \f$ X =
   x \f$. With two, then  \f$ R = [r_1, r_2] \f$,  \f$ X = [x_1 x_2]
   \f$.

   Remember that \f$ y \f$ and \f$ \lambda \f$ are relation from the Interaction, and have
   the same size.

*/
class FirstOrderR : public Relation {
 public:
  enum FirstOrderRDS { Xxx, Rrr, DSlinkSize };
  enum FirstOrderRVec { e, relationVectorsSize };

 protected:
  ACCEPT_SERIALIZATION(FirstOrderR);

  /** basic constructor
   *
   *  \param newType the type of the relation
   */
  FirstOrderR(RelationSubType newType) : Relation(RelationType::FirstOrder, newType) {};

  /* The following matrices are used if the relation is linear w.r.t to some
   * variables. If the matrices are allocated, the computation of the Jacobian
   * is not done.
   */

  /** The Jacobian of the constraints with respect to the generalized coordinates \f$ x\f$
   *  i.e. \f$ \nabla_x h(x,t,\lambda) \f$
   */
  std::shared_ptr<siconos::algebra::MapType> jacobianhOver_state_view_{nullptr};

  /** internal (optional) storage used for \f$ \nabla_x h(x,t,\lambda) \f$ */
  std::unique_ptr<std::vector<double>> jacobianhOver_state_internal_storage_{nullptr};

  /** True if \f$ \nabla_x h(x,t,\lambda) \f$ is a constant matrix */
  bool hasConstantJacobianhOver_state_{false};

  /** The Jacobian of the constraints with respect to the generalized coordinates \f$
   * \lambda\f$ i.e. \f$ \nabla_{\lambda} h(x,t,\lambda) \f$
   */
  std::shared_ptr<siconos::algebra::MapType> jacobianhOver_lambda_view_{nullptr};

  /** internal (optional) storage used for\f$ \nabla_{\lambda} h(x,t,\lambda) \f$*/
  std::unique_ptr<std::vector<double>> jacobianhOver_lambda_internal_storage_{nullptr};

  /** True if \f$ \nabla_{\lambda} h(x,t,\lambda) \f$ is a constant matrix */
  bool hasConstantJacobianhOver_lambda_{false};

  /** \f$ \nabla_x g(x,t,\lambda) \f$
   */
  std::shared_ptr<siconos::algebra::MapType> jacobiangOver_state_view_{nullptr};

  /** internal (optional) storage used for \f$ \nabla_x g(x,t,\lambda) \f$ */
  std::unique_ptr<std::vector<double>> jacobiangOver_state_internal_storage_{nullptr};

  /** True if \f$ \nabla_x g(x,t,\lambda) \f$ is a constant matrix */
  bool hasConstantJacobiangOver_state_{false};

  /** \f$ \nabla_{\lambda} g(x,t,\lambda) \f$ */
  std::shared_ptr<siconos::algebra::MapType> jacobiangOver_lambda_view_{nullptr};

  /** internal (optional) storage used for\f$ \nabla_{\lambda} g(x,t,\lambda) \f$*/
  std::unique_ptr<std::vector<double>> jacobiangOver_lambda_internal_storage_{nullptr};

  /** True if \f$ \nabla_{\lambda} g(x,t,\lambda) \f$ is a constant matrix */
  bool hasConstantJacobiangOver_lambda_{false};

 public:
  /** destructor */
  virtual ~FirstOrderR() noexcept = default;

  /** initialize the relation (check sizes, memory allocation ...)
   *
   *  \param inter the interaction using this relation
   */
  void initialize(Interaction &inter) override;

  /** \return True if \f$ \nabla_x h(x,t,\lambda) \f$ is taken into account */
  bool hasJacobianhOver_state() const { return jacobianhOver_state_view_ != nullptr; }

  /** \return True if \f$ \nabla_x h(x,t,\lambda) \f$ is taken into account */
  bool hasJacobianhOver_lambda() const override {
    return jacobianhOver_lambda_view_ != nullptr;
  }

  /** \return True if \f$ \nabla_x h(x,t,\lambda) \f$ is taken into account */
  bool hasJacobiangOver_state() const { return jacobiangOver_state_view_ != nullptr; }

  /** \return True if \f$ \nabla_x h(x,t,\lambda) \f$ is taken into account */
  bool hasJacobiangOver_lambda() const { return jacobiangOver_lambda_view_ != nullptr; }

  /** \return True if \f$ \nabla_x h(x,t,\lambda) \f$ is taken into account */
  bool hasConstantJacobiangOver_lambda() const { return hasConstantJacobiangOver_lambda_; }

  /** \return a read-only view on \f$ \nabla_x h(x,t,\lambda) \f$ matrix */
  inline const siconos::algebra::ConstMapType jacobianhOver_state() const {
    return siconos::algebra::ConstMapType(jacobianhOver_state_view_->data(),
                                          jacobianhOver_state_view_->rows(),
                                          jacobianhOver_state_view_->cols());
  }

  /** \return a read-only view on \f$ \nabla_{\lambda} h(x,t,\lambda) \f$ matrix */
  inline const siconos::algebra::ConstMapType jacobianhOver_lambda() const override {
    return siconos::algebra::ConstMapType(jacobianhOver_lambda_view_->data(),
                                          jacobianhOver_lambda_view_->rows(),
                                          jacobianhOver_lambda_view_->cols());
  }

  /** \return a read-only view on \f$ \nabla_x g(x,t,\lambda) \f$ matrix */
  inline const auto jacobiangOver_state() const {
    return siconos::algebra::ConstMapType(jacobiangOver_state_view_->data(),
                                          jacobiangOver_state_view_->rows(),
                                          jacobiangOver_state_view_->cols());
  }

  /** \return a read-only view on \f$ \nabla_{\lambda} g(x,t,\lambda) \f$ matrix */
  inline const auto jacobiangOver_lambda() const {
    return siconos::algebra::ConstMapType(jacobiangOver_lambda_view_->data(),
                                          jacobiangOver_lambda_view_->rows(),
                                          jacobiangOver_lambda_view_->cols());
  }

  virtual void display() const override;
  // Visitors stuff
  void accept(std::shared_ptr<siconos::internal::SiconosVisitor> tourist) const override;
};
}  // namespace siconos::modeling
#endif
