// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file SphereLDSPlanR.hpp
  \brief SphereLDS relation with a plan - Inherits from LagrangianScleronomousR
*/

#ifndef SphereLDSPlanR_h
#define SphereLDSPlanR_h

#include "LagrangianScleronomousR.hpp"

namespace siconos::collision::native::bodies {

class SphereLDSPlanR : public siconos::modeling::LagrangianScleronomousR,
                       public std::enable_shared_from_this<SphereLDSPlanR> {
 private:
  ACCEPT_SERIALIZATION(SphereLDSPlanR);

  /* Ax + By + Cz + D = 0 */
  double r, A, B, C, D, nN, nU;
  /* u ^ v  = n */
  double u1, u2, u3, v1, v2, v3, n1, n2, n3, ru1, ru2, ru3, rv1, rv2, rv3;

 public:
  /** Constructor
   *
   *  \param r disk radius
   *  \param A
   *  \param B
   *  \param C
   *  \param D
   */
  SphereLDSPlanR(double r, double A, double B, double C, double D);

  ~SphereLDSPlanR() noexcept = default;

  double distance(double, double, double, double);

  /**
   to compute the output y = h(q) of the Relation

   \param q coordinates of the dynamical systems involved in the relation
   \param y the resulting vector
 */
  void computeh(const siconos::algebra::BlockVector& q,
                Eigen::Ref<siconos::algebra::SiconosVector> y) override;

  /** Computes \f$ \nabla^\top_q h(q) \f$
   * \param q coordinates of the dynamical systems involved in the relation
   */
  void computeContactMatrix(const siconos::algebra::BlockVector& q) override;

  bool equal(double _A, double _B, double _C, double _D, double _r) const {
    return (A == _A && B == _B && C == _C && D == _D && r == _r);
  }

  virtual void accept(modeling::relations::Visitor& tourist) const override {
    tourist.visit(*this);
  }
};
}  // namespace siconos::collision::native::bodies
#endif /* SphereLDSPlanR_h */
