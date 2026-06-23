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

#include <cstddef>

#include "BlockVector.hpp"
#include "DynamicalSystem.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

namespace siconos::control {

typedef void (*AdditionalTermsEDfctU)(double, size_t, double*, size_t, double*, double*,
                                      size_t, double*);
}

void siconos::control::ControlLinearAdditionalTermsED::init(
    siconos::graphs::DynamicalSystemsGraph& DSG0,
    const siconos::modeling::NonSmoothDynamicalSystem& nsds,
    std::shared_ptr<siconos::simulation::TimeDiscretisation> td) {
  siconos::graphs::DynamicalSystemsGraph::VIterator dsvi, dsvdend;
  for (std::tie(dsvi, dsvdend) = DSG0.vertices(); dsvi != dsvdend; ++dsvi) {
    auto& ds = *DSG0.bundle(*dsvi);
    if (DSG0.pluginU.hasKey(*dsvi)) {
      DSG0.tmpXdot[*dsvi] = std::make_shared<siconos::algebra::SiconosVector>(ds.dimension());
    }
    if (DSG0.pluginJacgx.hasKey(*dsvi)) {
      DSG0.jacgx[*dsvi] =
          std::make_shared<siconos::algebra::SiconosMatrix>(ds.dimension(), ds.dimension());
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
      xdot.noalias() += DSG0.B.getRef(dsgVD) * DSG0.u.getRef(dsgVD);
    } else if (DSG0.pluginU.hasKey(dsgVD)) {
      auto& ds = *DSG0.bundle(dsgVD);
      auto& u = DSG0.u.getRef(dsgVD);
      auto tmpXdot = DSG0.tmpXdot[dsgVD];
      siconos::algebra::BlockVector xb;
      xb.insertPtr(ds.x());
      siconos::algebra::BlockVector xdotb;
      xdotb.insertPtr(tmpXdot);
      (DSG0.pluginU[dsgVD])(xb, t, u, xdotb);
      xdot += *tmpXdot;  // xdot += g(x, u)
    } else {
      THROW_EXCEPTION(
          "siconos::control::ControlLinearAdditionalTermsED :: input u but no B nor pluginU");
    }
  }
  // check whether the DynamicalSystem is an Observer
  if (DSG0.eVector.hasKey(dsgVD)) {
    assert(DSG0.L.hasKey(dsgVD));
    xdot.noalias() += *DSG0.L[dsgVD] * *DSG0.eVector[dsgVD];  // xdot += -L*e
  }
}

void siconos::control::ControlLinearAdditionalTermsED::addJacobianRhsContribution(
    siconos::graphs::DynamicalSystemsGraph& DSG0,
    const siconos::graphs::DynamicalSystemsGraph::VDescriptor& dsgVD, const double t,
    siconos::algebra::SiconosVector& jacRhs) {
  // check whether we have a system with a control input
  if (DSG0.pluginJacgx.hasKey(dsgVD)) {
    auto& ds = *DSG0.bundle(dsgVD);
    auto& u = DSG0.u.getRef(dsgVD);
    auto& tmpJacgx = DSG0.jacgx.getRef(dsgVD);
    siconos::algebra::BlockVector xb;
    xb.insertPtr(ds.x());
    DSG0.pluginJacgx[dsgVD](xb, t, u, tmpJacgx);
    const double* result_data = tmpJacgx.data();
    for (auto i = 0; i < jacRhs.size(); ++i) {
      jacRhs(i) +=
          result_data[i];  // remind that the jacobian is saved as a flat vector (column-major)
    }
  } else {
    THROW_EXCEPTION(
        "siconos::control::ControlLinearAdditionalTermsED :: input u but no B nor "
        "pluginJacgx");
  }
}
