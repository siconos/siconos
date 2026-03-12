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
#include "FEM-GlobalFrictionContact.hpp"

#include <memory>

#include "FEM-MoreauJeanGOSI.hpp"
#include "MoreauJeanGOSI.hpp"
#include "NewtonImpactFrictionNSL.hpp"
#include "Simulation.hpp"
#include "SimulationGraphs.hpp"
#include "SolidLinearTIDS.hpp"
#include "Tools.hpp"
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "SecondOrderDS.hpp"
#include "siconos_debug.h"

void siconos::mechanics::fem::nonsmooth_formulations::GlobalFrictionContact::compute_q() {
  if (_q->size() != _sizeGlobalOutput) _q->resize(_sizeGlobalOutput);

  siconos::algebra::Index offset = 0;
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  auto& DSG0 = *simulation()->nonSmoothDynamicalSystem()->dynamicalSystems();
  for (std::tie(dsi, dsend) = DSG0.vertices(); dsi != dsend; ++dsi) {
    auto ds = DSG0.bundle(*dsi);
    auto ds_size = ds->real_size();

    if (auto mjgosi =
            std::dynamic_pointer_cast<siconos::mechanics::fem::integrators::MoreauJeanGOSI>(
                DSG0.properties(DSG0.descriptor(ds)).osi)) {
      auto& ds_work_vectors = *DSG0.properties(DSG0.descriptor(ds)).workVectors;
      auto solid = std::dynamic_pointer_cast<siconos::mechanics::fem::SolidLinearTIDS>(ds);
      assert(solid);
      auto& vfree = *ds_work_vectors[tools::enum_to_index(
          siconos::mechanics::fem::integrators::MoreauJeanGOSI::wk_ds::vfree)];
      auto& sigmafree = *ds_work_vectors[tools::enum_to_index(
          siconos::mechanics::fem::integrators::MoreauJeanGOSI::wk_ds::residu_sigma_free)];

      // q_sigma = [v; sigma]
      siconos::algebra::SiconosVector qSigmafree{vfree.size() + sigmafree.size()};
      qSigmafree.head(solid->dimension()) = vfree;
      qSigmafree.tail(solid->stressDimension()) = sigmafree;
      _q->segment(offset, ds_size) = qSigmafree;

    } else if (auto mjgosi = std::dynamic_pointer_cast<siconos::integrators::MoreauJeanGOSI>(
                   DSG0.properties(DSG0.descriptor(ds)).osi)) {
      auto& ds_work_vectors = *DSG0.properties(DSG0.descriptor(ds)).workVectors;
      if ((std::dynamic_pointer_cast<siconos::modeling::SecondOrderDS>(ds))) {
        auto& vfree = *ds_work_vectors[tools::enum_to_index(
            siconos::integrators::MoreauJeanGOSI::wk_ds::vfree)];
        _q->segment(offset, ds_size) = vfree.head(ds_size);
      }
    } else
      THROW_EXCEPTION(
          "siconos::mechanics::fem::nonsmooth_formulations::GlobalFrictionContact::"
          "computeq: only implemented for a MoreauJeanGOSI integrator.");

    offset += ds_size;
  }
}

void siconos::mechanics::fem::nonsmooth_formulations::GlobalFrictionContact::
    compute_nslaw_contribution(siconos::graphs::InteractionsGraph& indexSet) {
  if (_b->size() != _sizeOutput) _b->resize(_sizeOutput);
  _b->setZero();
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  for (std::tie(ui, uiend) = indexSet.vertices(); ui != uiend; ++ui) {
    auto inter = indexSet.bundle(*ui);

    if (auto nslaw = std::dynamic_pointer_cast<siconos::modeling::NewtonImpactFrictionNSL>(
            inter->nonSmoothLaw()))
      _mu->push_back(nslaw->mu());

    // -- Dynamical systems of the current interaction
    auto ds1 = indexSet.properties(*ui).source;
    auto ds2 = indexSet.properties(*ui).target;

    // We need to be sure that the integrators of these DS are MoreauJeanGOSI
    auto& DSG0 = *simulation()->nonSmoothDynamicalSystem()->dynamicalSystems();
    auto osi1 = DSG0.properties(DSG0.descriptor(ds1)).osi;
    auto osi2 = DSG0.properties(DSG0.descriptor(ds2)).osi;
    [[maybe_unused]] bool is_ds2_mjgosi =
        std::dynamic_pointer_cast<siconos::integrators::MoreauJeanGOSI>(osi2) ||
        std::dynamic_pointer_cast<siconos::mechanics::fem::integrators::MoreauJeanGOSI>(osi2);

    if (auto mjgosi1 = std::dynamic_pointer_cast<siconos::integrators::MoreauJeanGOSI>(osi1)) {
      assert(is_ds2_mjgosi);
      mjgosi1->NonSmoothLawContributionToOutput(inter, *this);
    } else if (auto fem_mjgosi1 = std::dynamic_pointer_cast<
                   siconos::mechanics::fem::integrators::MoreauJeanGOSI>(osi1)) {
      assert(is_ds2_mjgosi);
      fem_mjgosi1->NonSmoothLawContributionToOutput(inter, *this);
    } else {
      THROW_EXCEPTION(
          "siconos::nonsmooth_formulations::GlobalFrictionContact::computeq. Not yet "
          "implemented for "
          "the given Integrator type ");
    }
    auto& osnsp_rhs = *(*indexSet.properties(*ui).workVectors)[tools::enum_to_index(
        siconos::integrators::MoreauJeanGOSI::wk_inter::osnsp_rhs)];
    auto pos = indexSet.properties(*ui).absolute_position;
    auto sizeY = inter->dimension();
    _b->segment(pos, sizeY) = osnsp_rhs.head(sizeY);
  }
  DEBUG_EXPR(siconos::algebra::print(*_b););
}

void siconos::mechanics::fem::nonsmooth_formulations::GlobalFrictionContact::
    update_dynamicalsystems_state() {
  auto& DSG0 = *simulation()->nonSmoothDynamicalSystem()->dynamicalSystems();

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = DSG0.vertices(); dsi != dsend; ++dsi) {
    auto ds = DSG0.bundle(*dsi);
    auto& ds_work_vectors = *DSG0.properties(*dsi).workVectors;
    auto pos = DSG0.properties(*dsi).absolute_position;
    if (auto solid = std::dynamic_pointer_cast<siconos::mechanics::fem::SolidLinearTIDS>(ds)) {
      auto& osi = static_cast<siconos::mechanics::fem::integrators::MoreauJeanGOSI&>(
          *DSG0.properties(*dsi).osi);
      auto& v_iter = osi.get_v_iter(ds_work_vectors);
      auto& stress_iter = osi.get_sigma_iter(ds_work_vectors);
      v_iter = _globalVelocities->segment(pos, solid->dimension());
      stress_iter =
          _globalVelocities->segment(pos + solid->dimension(), solid->stressDimension());
    } else if (auto sds = std::dynamic_pointer_cast<siconos::modeling::SecondOrderDS>(ds)) {
      auto& osi =
          static_cast<siconos::integrators::MoreauJeanGOSI&>(*DSG0.properties(*dsi).osi);
      auto& v_iter = osi.get_v_iter(ds_work_vectors);
      auto sizeDS = sds->dimension();
      v_iter = _globalVelocities->segment(pos, sizeDS);
    }
  }
}