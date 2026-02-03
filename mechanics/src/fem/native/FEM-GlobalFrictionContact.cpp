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

#include "FEM-MoreauJeanGOSI.hpp"
#include "OSNSMatrix.hpp"
#include "Simulation.hpp"
#include "SolidLinearTIDS.hpp"
#include "Tools.hpp"
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

void siconos::mechanics::fem::nonsmooth_formulations::GlobalFrictionContact::compute_q() {
  if (_q->size() != _sizeGlobalOutput) _q->resize(_sizeGlobalOutput);

  size_t offset = 0;
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  auto& DSG0 = *simulation()->nonSmoothDynamicalSystem()->dynamicalSystems();
  for (std::tie(dsi, dsend) = DSG0.vertices(); dsi != dsend; ++dsi) {
    auto ds = DSG0.bundle(*dsi);
    auto ds_size = ds->dimension();

    auto mjgosi =
        std::dynamic_pointer_cast<siconos::mechanics::fem::integrators::MoreauJeanGOSI>(
            DSG0.properties(DSG0.descriptor(ds)).osi);
    if (mjgosi) {
      auto& ds_work_vectors = *DSG0.properties(DSG0.descriptor(ds)).workVectors;
      auto solid = std::dynamic_pointer_cast<siconos::mechanics::fem::SolidLinearTIDS>(ds);
      assert(solid);
      auto dss = ds->real_size();
      auto& vfree = *ds_work_vectors[tools::enum_to_index(
          siconos::mechanics::fem::integrators::MoreauJeanGOSI::wk_ds::vfree)];
      auto& sigmafree = *ds_work_vectors[tools::enum_to_index(
          siconos::mechanics::fem::integrators::MoreauJeanGOSI::wk_ds::residu_sigma_free)];

      // q_sigma = [v; sigma]
      siconos::algebra::SiconosVector qSigmafree{vfree.size() + sigmafree.size()};
      qSigmafree.head(solid->dimension()) = vfree;
      qSigmafree.tail(solid->stressDimension()) = sigmafree;
      _q->segment(offset, dss) = qSigmafree;

    } else
      THROW_EXCEPTION(
          "siconos::mechanics::fem::nonsmooth_formulations::GlobalFrictionContact::"
          "computeq: only implemented for a MoreauJeanGOSI integrator.");

    offset += ds_size;
  }
}
void siconos::mechanics::fem::nonsmooth_formulations::GlobalFrictionContact::
    update_dynamicalsystems_state() {
  auto& DSG0 = *simulation()->nonSmoothDynamicalSystem()->dynamicalSystems();

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = DSG0.vertices(); dsi != dsend; ++dsi) {
    auto ds = DSG0.bundle(*dsi);
    auto solid = std::dynamic_pointer_cast<siconos::mechanics::fem::SolidLinearTIDS>(ds);
    assert(solid);
    auto velocity = solid->velocity();
    auto stress = solid->stress();
    DEBUG_PRINTF("ds.number() : %i \n", ds.number());
    DEBUG_EXPR(velocity->display(););
    DEBUG_EXPR(stress->display(););
    DEBUG_EXPR(_globalVelocities->display(););
    auto pos = DSG0.properties(*dsi).absolute_position;
    // std::cout << "_globalVelocities in postCompute:" << std::endl;
    // _globalVelocities->display();
    auto& osi = static_cast<siconos::mechanics::fem::integrators::MoreauJeanGOSI&>(
        *DSG0.properties(*dsi).osi);
    auto& ds_work_vectors = *DSG0.properties(*dsi).workVectors;
    auto& v_iter = osi.get_v_iter(ds_work_vectors);
    auto& stress_iter = osi.get_sigma_iter(ds_work_vectors);
    v_iter = _globalVelocities->segment(pos, solid->dimension());
    stress_iter =
        _globalVelocities->segment(pos + solid->dimension(), solid->stressDimension());
  }
}
