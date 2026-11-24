/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#include "ControlLinearAdditionalTermsTS.hpp"

void siconos::control::ControlLinearAdditionalTermsTS::addSmoothTerms(
    siconos::graphs::DynamicalSystemsGraph& DSG0,
    const siconos::graphs::DynamicalSystemsGraph::VDescriptor& dsgVD, const double h,
    siconos::algebra::SiconosVector& xfree) {
  // check whether we have a system with a control input
  if (DSG0.u.hasKey(dsgVD)) {
    assert(DSG0.B.hasKey(dsgVD));
    xfree.noalias() += h * *DSG0.B[dsgVD] * *DSG0.u[dsgVD];  // xfree += h*B*u
  }
  // check whether the DynamicalSystem is an Observer
  if (DSG0.eVector.hasKey(dsgVD)) {
    assert(DSG0.L.hasKey(dsgVD));
    xfree.noalias() += h * *DSG0.L[dsgVD] * *DSG0.eVector[dsgVD];  // xfree += -h*L*e
  }
}

void siconos::control::ControlLinearAdditionalTermsTS::addJacobianRhsContribution(
    siconos::graphs::DynamicalSystemsGraph& DSG0,
    const siconos::graphs::DynamicalSystemsGraph::VDescriptor& dsgVD, const double t,
    siconos::algebra::SiconosVector& jacRhs) {
  // nothing to be done here
}
