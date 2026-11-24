/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#include "SphereLDSPlanR.hpp"

#include <cmath>

#include "BlockVector.hpp"
#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "op3x3.h"  // for orthoBaseFromVector

siconos::collision::native::bodies::SphereLDSPlanR::SphereLDSPlanR(double r, double A,
                                                                   double B, double C,
                                                                   double D)
    : siconos::modeling::LagrangianScleronomousR{}, r{r}, A{A}, B{B}, C{C}, D{D} {
  n1 = A;
  n2 = B;
  n3 = C;

  nN = sqrt(A * A + B * B + C * C);

  if (orthoBaseFromVector(&n1, &n2, &n3, &u1, &u2, &u3, &v1, &v2, &v3))  // siconos::numerics
    THROW_EXCEPTION(
        "siconos::collision::native::bodies::SphereLDSPlanR::SphereLDSPlanR. Problem in "
        "calling "
        "orthoBaseFromVector");
  // r*u & r *v

  ru1 = r * u1;
  ru2 = r * u2;
  ru3 = r * u3;

  rv1 = r * v1;
  rv2 = r * v2;
  rv3 = r * v3;
}

double siconos::collision::native::bodies::SphereLDSPlanR::distance(double x, double y,
                                                                    double z, double rad) {
  return (fabs(A * x + B * y + C * z + D) / nN - rad);
}

void siconos::collision::native::bodies::SphereLDSPlanR::computeh(
    const siconos::algebra::BlockVector& q, Eigen::Ref<siconos::algebra::SiconosVector> y) {
  double q_0 = q(0);
  double q_1 = q(1);
  double q_2 = q(2);

  y(0) = distance(q_0, q_1, q_2, r);
};

void normalize(std::shared_ptr<siconos::algebra::SiconosVector>, unsigned int);

void siconos::collision::native::bodies::SphereLDSPlanR::computeContactMatrix(
    const siconos::algebra::BlockVector& q) {
  double theta = q(3);
  double phi = q(4);

  double cthe = cos(theta);
  double sthe = sin(theta);
  double cphi = cos(phi);
  double sphi = sin(phi);

  contactMatrix_view_->setValue(0, 0, n1);
  contactMatrix_view_->setValue(1, 0, u1);
  contactMatrix_view_->setValue(2, 0, v1);
  contactMatrix_view_->setValue(0, 1, n2);
  contactMatrix_view_->setValue(1, 1, u2);
  contactMatrix_view_->setValue(2, 1, v2);
  contactMatrix_view_->setValue(0, 2, n3);
  contactMatrix_view_->setValue(1, 2, u3);
  contactMatrix_view_->setValue(2, 2, v3);
  contactMatrix_view_->setValue(0, 3, 0);
  contactMatrix_view_->setValue(1, 3, -rv1 * cphi - rv2 * sphi);
  contactMatrix_view_->setValue(2, 3, ru1 * cphi + ru2 * sphi);
  contactMatrix_view_->setValue(0, 4, 0);
  contactMatrix_view_->setValue(1, 4, -rv3);
  contactMatrix_view_->setValue(2, 4, ru3);
  contactMatrix_view_->setValue(0, 5, 0);
  contactMatrix_view_->setValue(1, 5, -rv3 * cthe + rv2 * cphi * sthe - rv1 * sphi * sthe);
  contactMatrix_view_->setValue(2, 5, ru3 * cthe + ru1 * sphi * sthe - ru2 * cphi * sthe);
}
