/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#include "OneStepIntegrator.hpp"

#include "DynamicalSystem.hpp"
#include "EventsManager.hpp"  // IWYU pragma: keep
#include "Interaction.hpp"
#include "Relation.hpp"  // IWYU pragma: keep
#include "SiconosException.hpp"
#include "SiconosVector.hpp"
#include "Simulation.hpp"
#include "Tools.hpp"  // enum_to_string

std::shared_ptr<siconos::algebra::blocks::SharedVector>
siconos::integrators::OneStepIntegrator::_initializeDSWorkVectors(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds) {
  const siconos::graphs::DynamicalSystemsGraph::VDescriptor& dsv =
      _dynamicalSystemsGraph->descriptor(ds);

  // Create new work buffers, store in the graph
  auto wv = std::make_shared<siconos::algebra::blocks::SharedVector>();
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
    auto sim = simulation();
    _extraAdditionalTerms->init(*sim->nonSmoothDynamicalSystem()->topology()->dSG(0),
                                *sim->nonSmoothDynamicalSystem(),
                                sim->eventsManager()->timeDiscretisation());
  }

  initialize_nonsmooth_problems();
  _isInitialized = true;
}

void siconos::integrators::OneStepIntegrator::updateAndSwapAllOutput(double time) {
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  auto& indexSet0 = *simulation()->nonSmoothDynamicalSystem()->topology()->indexSet0();
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

  // Multi-step methods (SchatzmanPaoliOSI): this function is called at the beginning of
  // each time step (TimeStepping::initializeNewtonSolve), after the memories have been
  // pushed at the end of the previous step. Once the memories hold _steps values, the
  // output must only refresh the latest entry: pushing it again would replace y_{k-1}
  // with y_k. (For one-step methods, the memory has one entry and both are the same.)
  const bool refreshOnly =
      _steps > 1 && inter.yMemory(inter.lowerLevelForOutput()).nbVectorsInMemory() >= _steps;

  if (_steps > 1 && !refreshOnly)  // Multi--step methods: fill the memory
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
  // FP: thought it was useless but for NewtonEuler5D for example,
  // we need Nc to be set before any call to computeJach (below)
  // so we need a first call to computeOutput (for h)
  inter.computeOutput(time, 0);

  // prepare the gradients
  // --> Note FP:
  // useless. This is done in computeOutput ?
  // No, it's not (why ?)
  // This means that when we call computeOutput, operators depending
  // on the state might not be uptodate. See V.A remark in NewtonEulerR computeOutput
  inter.relation()->computeJach(time, inter);
  for (unsigned int i = 0; i < inter.upperLevelForOutput() + 1; ++i) {
    inter.computeOutput(time, i);
  }
  if (refreshOnly)
    inter.updateLatestInMemory();
  else
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

void siconos::integrators::OneStepIntegrator::resetNonSmoothPart(
    siconos::algebra::blocks::size_type level) {
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    _dynamicalSystemsGraph->bundle(*dsi)->resetNonSmoothPart(level);
  }
}

std::shared_ptr<siconos::algebra::SiconosMatrix>
siconos::integrators::OneStepIntegrator::iterationMatrix(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds) {
  assert(ds && "siconos::integrators::OneStepIntegrator::iterationMatrix(ds): ds == nullptr.");
  return _dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds))
      .iterationMatrix;
}

std::shared_ptr<siconos::algebra::SiconosDenseLUMatrix>
siconos::integrators::OneStepIntegrator::LUiterationMatrix(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds) {
  assert(ds && "siconos::integrators::OneStepIntegrator::iterationMatrix(ds): ds == nullptr.");
  return _dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).LUW;
}

void siconos::integrators::OneStepIntegrator::updateOutput(double time, unsigned int level) {
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  auto& indexSet0 = *simulation()->nonSmoothDynamicalSystem()->topology()->indexSet0();
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
  /** VA. 16/02/2017 This should normally be done only for interaction managed by the osi
   */
  for (auto level = _levelMinForInput; level < _levelMaxForInput + 1; level++)
    updateInput(time, level);
}

void siconos::integrators::OneStepIntegrator::updateInput(double time, unsigned int level) {
  // resetNonSmoothPart(level);
  // We compute input using lambda(level).
  // simulation()->>nonSmoothDynamicalSystem()->updateInput(time,level);
  // resetNonSmoothPart(level);
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  auto& indexSet0 = *simulation()->nonSmoothDynamicalSystem()->topology()->indexSet0();
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
    std::shared_ptr<siconos::modeling::Interaction> inter,
    siconos::graphs::InteractionsGraph::size_type i) {
  THROW_EXCEPTION(
      "OneStepIntegrator::addInteractionInIndexSet - Should be called at this level");
  return 0;
};

bool siconos::integrators::OneStepIntegrator::removeInteractionFromIndexSet(
    std::shared_ptr<siconos::modeling::Interaction> inter,
    siconos::graphs::InteractionsGraph::size_type i) {
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
