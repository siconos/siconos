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

#include "ControlLinearAdditionalTermsED.hpp"

#include "DynamicalSystem.hpp"
#include "SiconosMatrixVectorOp.hpp"
#include "SiconosVector.hpp"
#include "SiconosMatrix.hpp"

namespace siconos::control {

typedef void (*AdditionalTermsEDfctU)(double, unsigned, double*, unsigned, double*, double*,
                                      unsigned, double*);
}

void siconos::control::ControlLinearAdditionalTermsED::init(
    siconos::graphs::DynamicalSystemsGraph& DSG0,
    const siconos::modeling::NonSmoothDynamicalSystem& nsds,
     std::shared_ptr<siconos::simulation::TimeDiscretisation> td) {
  siconos::graphs::DynamicalSystemsGraph::VIterator dsvi, dsvdend;
  for (std::tie(dsvi, dsvdend) = DSG0.vertices(); dsvi != dsvdend; ++dsvi) {
    auto& ds = *DSG0.bundle(*dsvi);
    if (DSG0.pluginU.hasKey(*dsvi)) {
      DSG0.tmpXdot[*dsvi] =
          std::make_shared<siconos::algebra::SiconosVector>(ds.getx().size());
    }
    if (DSG0.pluginJacgx.hasKey(*dsvi)) {
      DSG0.jacgx[*dsvi] =
          std::make_shared<siconos::algebra::SiconosMatrix>(ds.getx().size(), ds.getx().size());
    }
  }
}

void siconos::control::ControlLinearAdditionalTermsED::addSmoothTerms(
    siconos::graphs::DynamicalSystemsGraph& DSG0,
    const siconos::graphs::DynamicalSystemsGraph::VDescriptor& dsgVD, const double t,
    siconos::algebra::SiconosVector& xdot) {
  // check whether we have a system with a control input
  if (DSG0.u.hasKey(dsgVD)) {
    if (DSG0.B.hasKey(dsgVD)) {
      siconos::algebra::prod(DSG0.B.getRef(dsgVD), DSG0.u.getRef(dsgVD), xdot,
                             false);  // xdot += B*u
    } else if (DSG0.pluginU.hasKey(dsgVD)) {
      auto& ds = *DSG0.bundle(dsgVD);
      auto& u = DSG0.u.getRef(dsgVD);
      auto& tmpXdot = DSG0.tmpXdot.getRef(dsgVD);
      ((AdditionalTermsEDfctU)DSG0.pluginU.getRef(dsgVD).fPtr)(
          t, xdot.size(), const_cast<double*>(ds.getx().data()), u.size(), u.data(), tmpXdot.data(),
          ds.getz().size(), const_cast<double*>(ds.getz().data()));
      xdot += tmpXdot;  // xdot += g(x, u)
    } else {
      THROW_EXCEPTION(
          "siconos::control::ControlLinearAdditionalTermsED :: input u but no B nor pluginU");
    }
  }
  // check whether the DynamicalSystem is an Observer
  if (DSG0.e.hasKey(dsgVD)) {
    assert(DSG0.L.hasKey(dsgVD));
    siconos::algebra::prod(*DSG0.L[dsgVD], *DSG0.e[dsgVD], xdot, false);  // xdot += -L*e
  }
}

void siconos::control::ControlLinearAdditionalTermsED::addJacobianRhsContribution(
    siconos::graphs::DynamicalSystemsGraph& DSG0,
    const siconos::graphs::DynamicalSystemsGraph::VDescriptor& dsgVD, const double t,
    siconos::algebra::SiconosMatrix& jacRhs) {
  // check whether we have a system with a control input
  if (DSG0.pluginJacgx.hasKey(dsgVD)) {
    auto& ds = *DSG0.bundle(dsgVD);
    auto& u = DSG0.u.getRef(dsgVD);
    auto& tmpJacgx = DSG0.jacgx.getRef(dsgVD);
    ((AdditionalTermsEDfctU)DSG0.pluginJacgx.getRef(dsgVD).fPtr)(
        t, ds.getx().size(), const_cast<double*>(ds.getx().data()), u.size(), u.data(), tmpJacgx.data(),
        ds.getz().size(), const_cast<double*>(ds.getz().data()));
    jacRhs += tmpJacgx;  // JacRhs += \nabla_x g(x, u)
  } else {
    THROW_EXCEPTION(
        "siconos::control::ControlLinearAdditionalTermsED :: input u but no B nor pluginU");
  }
}
