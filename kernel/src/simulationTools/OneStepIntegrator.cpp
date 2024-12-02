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

#include "OneStepIntegrator.hpp"

#include "BlockVector.hpp"
#include "DynamicalSystem.hpp"
#include "EventsManager.hpp"
#include "Interaction.hpp"
#include "Relation.hpp"
#include "SiconosException.hpp"
#include "SiconosVector.hpp"
#include "SiconosVisitor.hpp"
#include "Simulation.hpp"
#include "Tools.hpp"  // enum_to_string
#include "Topology.hpp"

// FP WIP
// struct siconos::integrators::OneStepIntegrator::IterationMatrixVisitor : public
// siconos::internal::SiconosVisitor {

//   std::shared_ptr<siconos::algebra::SimpleMatrix> visit(const MoreauJeanOSI& osi,
//   std::shared_ptr<siconos::modeling::DynamicalSystem> ds) const
//   {
//     return osi.W(ds);
//   }

// };

std::shared_ptr<std::vector<std::shared_ptr<siconos::algebra::SiconosVector>>>
siconos::integrators::OneStepIntegrator::_initializeDSWorkVectors(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds) {
  const siconos::graphs::DynamicalSystemsGraph::VDescriptor& dsv =
      _dynamicalSystemsGraph->descriptor(ds);

  // Create new work buffers, store in the graph
  auto wv = std::make_shared<std::vector<std::shared_ptr<siconos::algebra::SiconosVector>>>();
  auto wm = std::make_shared<std::vector<std::shared_ptr<siconos::algebra::SiconosMatrix>>>();
  _dynamicalSystemsGraph->properties(dsv).workVectors = wv;
  _dynamicalSystemsGraph->properties(dsv).workMatrices = wm;

  // Initialize memory buffers
  ds->initMemory(getSizeMem());

  // Force dynamical system to its initial state
  ds->resetToInitialState();

  return wv;
}

void siconos::integrators::OneStepIntegrator::initialize() {
  if (_extraAdditionalTerms) {
    _extraAdditionalTerms->init(*_simulation->nonSmoothDynamicalSystem()->topology()->dSG(0),
                                *_simulation->nonSmoothDynamicalSystem(),
                                _simulation->eventsManager()->timeDiscretisation());
  }

  initialize_nonsmooth_problems();
  _isInitialized = true;
}
void siconos::integrators::OneStepIntegrator::updateAndSwapAllOutput(double time) {
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  auto& indexSet0 = *_simulation->nonSmoothDynamicalSystem()->topology()->indexSet0();
  for (std::tie(ui, uiend) = indexSet0.vertices(); ui != uiend; ++ui) {
    if (!checkInteractionOSI(indexSet0, ui)) continue;
    auto& inter = *indexSet0.bundle(*ui);
    updateAndSwapAllOutput(inter, time);
  }
}

void siconos::integrators::OneStepIntegrator::OneStepIntegrator::updateAndSwapAllOutput(
    siconos::modeling::Interaction& inter, double time) {
  // - compute interaction output (y) for all levels
  // - swaps in memory

  // This function:
  // - is called during osi->initialize()
  // - may be called by simulation when an interaction is inserted on the fly
  //    for instance:
  //      - contact detection
  //      - inter = ew interaction + link with ds
  //      - simu->osi->initializeWorkVectorsForInteraction(inter)
  //      - simu->osi->updateAndSwapAllOutput()

  if (_steps > 1)  // Multi--step methods
  {
    // Compute the old Values of Output with stored values in Memory
    for (unsigned int k = 0; k < _steps - 1; k++) {
      /** ComputeOutput to fill the Memory
       * We assume the state x is stored in xMemory except for the  initial
       * condition which has not been swap yet.
       */
      //        relation()->LinkDataFromMemory(k);
      for (unsigned int i = 0; i < inter.upperLevelForOutput() + 1; ++i) {
        inter.computeOutput(time, i);
        //_yMemory[i]->swap(*_y[i]);
      }
    }
    inter.swapInMemory();
  }

  
  // Compute a first value for the output
  // VA 10/04/2024 What is the interest of the following line ?
  inter.computeOutput(time, 0);

  // prepare the gradients
  inter.relation()->computeJach(time, inter);
  for (unsigned int i = 0; i < inter.upperLevelForOutput() + 1; ++i) {
    inter.computeOutput(time, i);
  }
  inter.swapInMemory();
}

void siconos::integrators::OneStepIntegrator::_check_and_update_interaction_levels(
    siconos::modeling::Interaction& inter) {
  bool isInitializationNeeded = false;
  if (!(inter.lowerLevelForOutput() <= _levelMinForOutput &&
        inter.upperLevelForOutput() >= _levelMaxForOutput)) {
    inter.setLowerLevelForOutput(_levelMinForOutput);
    inter.setUpperLevelForOutput(_levelMaxForOutput);
    isInitializationNeeded = true;
  }

  if (!(inter.lowerLevelForInput() <= _levelMinForInput &&
        inter.upperLevelForInput() >= _levelMaxForInput)) {
    inter.setLowerLevelForInput(_levelMinForInput);
    inter.setUpperLevelForInput(_levelMaxForInput);
    isInitializationNeeded = true;
  }

  if (isInitializationNeeded) inter.reset();
}

void siconos::integrators::OneStepIntegrator::resetAllNonSmoothParts() {
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    _dynamicalSystemsGraph->bundle(*dsi)->resetAllNonSmoothParts();
  }
}

void siconos::integrators::OneStepIntegrator::resetNonSmoothPart(unsigned int level) {
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    _dynamicalSystemsGraph->bundle(*dsi)->resetNonSmoothPart(level);
  }
}

void siconos::integrators::OneStepIntegrator::updateOutput(double time, unsigned int level) {
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  auto& indexSet0 = *_simulation->nonSmoothDynamicalSystem()->topology()->indexSet0();
  for (std::tie(ui, uiend) = indexSet0.vertices(); ui != uiend; ++ui) {
    if (!checkInteractionOSI(indexSet0, ui)) continue;
    auto& inter = *indexSet0.bundle(*ui);
    assert(inter.lowerLevelForOutput() <= level);
    assert(inter.upperLevelForOutput() >= level);
    inter.computeOutput(time, level);
  }
}

void siconos::integrators::OneStepIntegrator::updateOutput(double time) {
  for (auto level = _levelMinForOutput; level <= _levelMaxForOutput; ++level)
    updateOutput(time, level);
}

void siconos::integrators::OneStepIntegrator::updateInput(double time) {
  for (auto level = _levelMinForInput; level < _levelMaxForInput + 1; level++)
    updateInput(time, level);
}

void siconos::integrators::OneStepIntegrator::updateInput(double time, unsigned int level) {
  // resetNonSmoothPart(level);
  // We compute input using lambda(level).
  //_simulation->nonSmoothDynamicalSystem()->updateInput(time,level);
  // resetNonSmoothPart(level);
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  auto& indexSet0 = *_simulation->nonSmoothDynamicalSystem()->topology()->indexSet0();
  for (std::tie(ui, uiend) = indexSet0.vertices(); ui != uiend; ++ui) {
    if (!checkInteractionOSI(indexSet0, ui)) continue;
    auto& inter = *indexSet0.bundle(*ui);
    assert(inter.lowerLevelForInput() <= level);
    assert(inter.upperLevelForInput() >= level);
    inter.computeInput(time, level);
  }
}

double siconos::integrators::OneStepIntegrator::computeResiduOutput(
    double time, std::shared_ptr<siconos::graphs::InteractionsGraph> indexSet) {
  double residu = 0.0;
  THROW_EXCEPTION(
      "siconos::integrators::OneStepIntegrator::computeResiduOutput not implemented for "
      "integrator of type " +
      std::to_string(
          static_cast<std::underlying_type<IntegratorType>::type>(_integratorType)));
  return residu;
}

double siconos::integrators::OneStepIntegrator::computeResidu() {
  // default : error
  THROW_EXCEPTION("OneStepIntegrator::computeResidu not implemented for integrator of type " +
                  std::to_string(static_cast<std::underlying_type<IntegratorType>::type>(
                      _integratorType)));
  return 0.0;
}

void siconos::integrators::OneStepIntegrator::computeFreeState() {
  // default : error
  THROW_EXCEPTION(
      "OneStepIntegrator::computeFreeState not implemented for integrator of type " +
      std::to_string(
          static_cast<std::underlying_type<IntegratorType>::type>(_integratorType)));
}

void siconos::integrators::OneStepIntegrator::computeFreeOutput(
    siconos::graphs::InteractionsGraph::VDescriptor& vertex_inter,
    siconos::nonsmooth_formulations::OneStepNSProblem* osnsp) {
  // default : error
  THROW_EXCEPTION(
      "OneStepIntegrator::computeFreeOutput not implemented for integrator of type " +
      std::to_string(
          static_cast<std::underlying_type<IntegratorType>::type>(_integratorType)));
}

double siconos::integrators::OneStepIntegrator::computeResiduInput(
    double time, std::shared_ptr<siconos::graphs::InteractionsGraph> indexSet) {
  double residu = 0.0;
  THROW_EXCEPTION(
      "siconos::integrators::OneStepIntegrator::computeResiduInput not implemented for "
      "integrator of type " +
      std::to_string(
          static_cast<std::underlying_type<IntegratorType>::type>(_integratorType)));
  return residu;
}

bool siconos::integrators::OneStepIntegrator::addInteractionInIndexSet(
    std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int i) {
  THROW_EXCEPTION(
      "OneStepIntegrator::addInteractionInIndexSet - Should be called at this level");
  return 0;
};

bool siconos::integrators::OneStepIntegrator::removeInteractionFromIndexSet(
    std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int i) {
  THROW_EXCEPTION(
      "OneStepIntegrator::removeInteractionFromIndexSet - Should not be called at this "
      "level");
  return 0;
};

void siconos::integrators::OneStepIntegrator::display() const {
  std::cout << "==== OneStepIntegrator display =====\n";
  std::cout << "| _integratorType : " << siconos::tools::enum_to_string(_integratorType)
            << "\n";
  std::cout << "| _sizeMem: " << _sizeMem << "\n";
  std::cout << "====================================\n";
}
