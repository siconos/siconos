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
#include "ProjectionOnConstraintsTools.hpp"

#include <algorithm>

#include "LagrangianR.hpp"
#include "NewtonEulerDS.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "SecondOrderDS.hpp"
#include "SiconosException.hpp"

namespace {
// Lagrangian systems (LagrangianDS, LagrangianSparseDS) are the SecondOrderDS that are not
// NewtonEulerDS. Their positions q and level-0 input p(0) live in the same space.
std::shared_ptr<siconos::modeling::SecondOrderDS> lagrangian(
    const std::shared_ptr<siconos::modeling::DynamicalSystem>& ds) {
  if (std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds)) return nullptr;
  return std::dynamic_pointer_cast<siconos::modeling::SecondOrderDS>(ds);
}
}  // namespace

void siconos::simulation::saveProjectionOrigin(
    const siconos::modeling::NonSmoothDynamicalSystem& nsds, ProjectionOrigin& a) {
  auto dsGraph = nsds.dynamicalSystems();
  for (auto vi : *dsGraph) {
    auto ds = dsGraph->bundle(vi);
    auto lds = lagrangian(ds);
    if (!lds)
      THROW_EXCEPTION(
          "siconos::simulation::saveProjectionOrigin: ProjectionUpdate::from_origin is only "
          "implemented for Lagrangian systems (LagrangianDS, LagrangianSparseDS).");
    a.insert_or_assign(ds.get(), siconos::algebra::SiconosVector(*lds->q()));
  }
}

void siconos::simulation::shiftProjectionRhs(siconos::graphs::InteractionsGraph& indexSet,
                                             const ProjectionOrigin& a) {
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  for (std::tie(ui, uiend) = indexSet.vertices(); ui != uiend; ++ui) {
    auto inter = indexSet.bundle(*ui);
    auto lr = std::dynamic_pointer_cast<siconos::modeling::LagrangianR>(inter->relation());
    if (!lr)
      THROW_EXCEPTION(
          "siconos::simulation::shiftProjectionRhs: ProjectionUpdate::from_origin is only "
          "implemented for Lagrangian relations.");

    // (a - q_j), concatenated over the dynamical systems of the interaction, in the order
    // of the columns of H (source, then target if different).
    const auto& prop = indexSet.properties(*ui);
    std::vector<std::shared_ptr<siconos::modeling::DynamicalSystem>> dss{prop.source};
    if (prop.target && prop.target != prop.source) dss.push_back(prop.target);

    auto H = lr->H();  // G(q_j), computed by computeCriteria()
    siconos::algebra::SiconosVector shift(H.cols());
    siconos::algebra::Index pos = 0;
    for (auto& ds : dss) {
      auto lds = lagrangian(ds);
      auto it = a.find(ds.get());
      if (!lds || it == a.end())
        THROW_EXCEPTION(
            "siconos::simulation::shiftProjectionRhs: no projection origin for a dynamical "
            "system of the interaction.");
      auto n = lds->q()->size();
      shift.segment(pos, n) = it->second - *lds->q();
      pos += n;
    }
    if (pos != H.cols())
      THROW_EXCEPTION(
          "siconos::simulation::shiftProjectionRhs: size of H inconsistent with the "
          "dynamical systems of the interaction.");

    // y0 = g(q_j) + G(q_j) (a - q_j): linearization at q_j of g evaluated at a, with
    // the Jacobian at q_j.
    inter->y(0)->head(H.rows()) += H * shift;
  }
}

double siconos::simulation::updatePositionsAfterProjection(
    const siconos::modeling::NonSmoothDynamicalSystem& nsds, const ProjectionOrigin& a,
    ProjectionUpdate mode) {
  double dqmax = 0.0;
  auto dsGraph = nsds.dynamicalSystems();
  for (auto vi : *dsGraph) {
    auto ds = dsGraph->bundle(vi);
    if (auto neds = std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds)) {
      if (mode != ProjectionUpdate::incremental)
        THROW_EXCEPTION(
            "siconos::simulation::updatePositionsAfterProjection: "
            "ProjectionUpdate::from_origin is not implemented for NewtonEulerDS.");
      if (neds->p(0)) {
        *neds->q() += neds->p_read(0);
        dqmax = std::max(dqmax, neds->p_read(0).lpNorm<Eigen::Infinity>());
      }
      neds->normalizeq();
      neds->computeT(neds->q_read());
    } else if (auto lds = lagrangian(ds)) {
      if (!lds->p(0)) continue;
      auto& q = *lds->q();
      if (mode == ProjectionUpdate::incremental) {
        // q_{j+1} = q_j + G_j^T lambda_j
        q += lds->p_read(0);
        dqmax = std::max(dqmax, lds->p_read(0).lpNorm<Eigen::Infinity>());
      } else {
        // q_{j+1} = a + G_j^T mu_j
        auto it = a.find(ds.get());
        if (it == a.end())
          THROW_EXCEPTION(
              "siconos::simulation::updatePositionsAfterProjection: no projection origin "
              "for a dynamical system.");
        siconos::algebra::SiconosVector qnew = it->second + lds->p_read(0);
        dqmax = std::max(dqmax, (qnew - q).lpNorm<Eigen::Infinity>());
        q = qnew;
      }
    } else
      THROW_EXCEPTION(
          "siconos::simulation::updatePositionsAfterProjection: Ds is not from NewtonEulerDS "
          "neither from LagrangianDS.");
  }
  return dqmax;
}
