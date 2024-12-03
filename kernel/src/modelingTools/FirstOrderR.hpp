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

   - FirstOrderNonlinearR: for fully nonlinear relations: \f$ y = h(t, X,
   \lambda, Z) \f$ , \f$ R = g(t, X, \lambda, Z) \f$ .
   - FirstOrderType2R: specialization with \f$ y = h(t, X, \lambda, Z) \f$ , \f$ R
   = g(t, \lambda, Z) \f$ .
   - FirstOrderType1R: further specialization with \f$ y = h(t, X, Z) \f$ , \f$ R
   = g(t, \lambda, Z) \f$ .
   - FirstOrderLinearR: linear case: \f$ y = C(t)x + D(t)\lambda + F(t) z +
   e \f$ , \f$ R = B(t)\lambda \f$ .
   - FirstOrderLinearTIR: time-invariant linear case: \f$ y = Cx + D\lambda + F z
   + e \f$ , \f$ R = B\lambda \f$ .

   If the relation involves only one DynamicalSystem, then  \f$ R = r \f$ ,  \f$ X =
   x \f$ , and  \f$ Z = z \f$ . With two, then  \f$ R = [r_1, r_2] \f$,  \f$ X = [x_1 x_2]
   \f$, and  \f$ Z = [z_1 z_2] \f$ .

   Remember that \f$ y \f$ and \f$ \lambda \f$ are relation from the Interaction, and have
   the same size.

*/
class FirstOrderR : public Relation {
 public:
  enum FirstOrderRDS { x, z, r, DSlinkSize };
  enum FirstOrderRVec { e, relationVectorsSize };
  enum FirstOrderRMat { mat_C, mat_D, mat_F, mat_B, mat_K, relationMatricesSize };

 protected:
  ACCEPT_SERIALIZATION(FirstOrderR);

  /** basic constructor
   *
   *  \param newType the type of the relation
   */
  FirstOrderR(RelationSubType newType) : Relation(RelationType::FirstOrder, newType){};

  /* The following matrices are used if the relation is linear w.r.t to some
   * variables. If the matrices are allocated, the computation of the Jacobian
   * is not done.
   */

  /** A matrix to store the constant Jacobian of h(t, X, lambda, Z) w.r.t X */
  std::shared_ptr<siconos::algebra::MapType> _C{nullptr};

  /** A matrix to store the constant Jacobian of h(t, X, lambda, Z) w.r.t lambda
   */
  std::shared_ptr<siconos::algebra::SiconosMatrix> _D{nullptr};

  /** A matrix to store the constant Jacobian of h(t, X, lambda, Z) w.r.t Z */
  std::shared_ptr<siconos::algebra::SiconosMatrix> _F{nullptr};

  /** A matrix to store the constant Jacobian of g(t, X, lambda, Z) w.r.t lambda
   */
  std::shared_ptr<siconos::algebra::MapType> _B{nullptr};

  /** A matrix to store the constant Jacobian of g(t, X, lambda, Z) w.r.t X */
  std::shared_ptr<siconos::algebra::SiconosMatrix> _K{nullptr};

 public:
  /** destructor
   */
  virtual ~FirstOrderR() noexcept = default;

  /** initialize the relation (check sizes, memory allocation ...)
   *
   *  \param inter the interaction using this relation
   */
  void initialize(Interaction &inter) override;

  /** set C to pointer newC
   *
   *  \param newC the C matrix
   */
  inline void setCPtr(std::shared_ptr<siconos::algebra::SiconosMatrix> newC) {
    _C = std::make_shared<siconos::algebra::MapType>(newC->data(), newC->rows(), newC->cols());
  }

  /** set B to pointer newB
   *
   *  \param newB the B matrix
   */
  inline void setBPtr(std::shared_ptr<siconos::algebra::SiconosMatrix> newB) {
    _B = std::make_shared<siconos::algebra::MapType>(newB->data(), newB->rows(), newB->cols());
  }

  /** set D to pointer newPtr
   *
   *  \param newD the D matrix
   */
  inline void setDPtr(std::shared_ptr<siconos::algebra::SiconosMatrix> newD) { _D = newD; }

  /** set F to pointer newPtr
   *
   *  \param newF the F matrix
   */
  inline void setFPtr(std::shared_ptr<siconos::algebra::SiconosMatrix> newF) { _F = newF; }

  /** get C
   *
   *  \return C matrix
   */
  inline std::shared_ptr<siconos::algebra::MapType> C() const override { return _C; }
  // Note FP: final would be better than override but swig cannot handle it.

  /** get H
   *
   *  \return C matrix
   */
  inline std::shared_ptr<siconos::algebra::MapType> H() const override { return _C; }

  /** get D
   *
   *  \return D matrix
   */
  inline std::shared_ptr<siconos::algebra::SiconosMatrix> D() const { return _D; }

  /** get F
   *
   *  \return F matrix
   */
  inline std::shared_ptr<siconos::algebra::SiconosMatrix> F() const { return _F; }

  /** get B
   *
   *  \return B matrix
   */
  inline std::shared_ptr<siconos::algebra::MapType> B() const { return _B; }

  /** get K
   *
   *  \return K matrix
   */
  inline std::shared_ptr<siconos::algebra::SiconosMatrix> K() const { return _K; }
  // Visitors stuff
  void accept(std::shared_ptr<siconos::internal::SiconosVisitor> tourist) const override;
};
}  // namespace siconos::modeling
#endif
