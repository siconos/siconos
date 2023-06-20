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

#include "ControlLinearAdditionalTermsTS.hpp"

#include "SiconosMatrixVectorOp.hpp"

void siconos::control::ControlLinearAdditionalTermsTS::addSmoothTerms(
    siconos::graphs::DynamicalSystemsGraph& DSG0,
    const siconos::graphs::DynamicalSystemsGraph::VDescriptor& dsgVD, const double h,
    siconos::algebra::SiconosVector& xfree) {
  // check whether we have a system with a control input
  if (DSG0.u.hasKey(dsgVD)) {
    assert(DSG0.B.hasKey(dsgVD));
    siconos::algebra::prod(h, *DSG0.B[dsgVD], *DSG0.u[dsgVD], xfree, false);  // xfree += h*B*u
  }
  // check whether the DynamicalSystem is an Observer
  if (DSG0.e.hasKey(dsgVD)) {
    assert(DSG0.L.hasKey(dsgVD));
    siconos::algebra::prod(h, *DSG0.L[dsgVD], *DSG0.e[dsgVD], xfree,
                           false);  // xfree += -h*L*e
  }
}

void siconos::control::ControlLinearAdditionalTermsTS::addJacobianRhsContribution(
    siconos::graphs::DynamicalSystemsGraph& DSG0,
    const siconos::graphs::DynamicalSystemsGraph::VDescriptor& dsgVD, const double t,
    siconos::algebra::SiconosMatrix& jacRhs) {
  // nothing to be done here
}
