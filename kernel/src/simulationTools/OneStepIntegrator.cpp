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

#include "OneStepIntegrator.hpp"
#include "Simulation.hpp"
#include "DynamicalSystem.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "ExtraAdditionalTerms.hpp"
#include "Relation.hpp"
#include "EventsManager.hpp"
#include <SiconosConfig.h>
#include <functional>
using namespace std::placeholders;

SP::VectorOfVectors
OneStepIntegrator::_initializeDSWorkVectors(SP::DynamicalSystem ds)
{
  const DynamicalSystemsGraph::VDescriptor& dsv =
    _dynamicalSystemsGraph->descriptor(ds);

  // Create new work buffers, store in the graph
  SP::VectorOfVectors wv = std::make_shared<VectorOfVectors>();
  SP::VectorOfMatrices wm = std::make_shared<VectorOfMatrices>();
  _dynamicalSystemsGraph->properties(dsv).workVectors = wv;
  _dynamicalSystemsGraph->properties(dsv).workMatrices = wm;

  // Initialize memory buffers
  ds->initMemory(getSizeMem());

  // Force dynamical system to its initial state
  ds->resetToInitialState();

  return wv;
}


void OneStepIntegrator::initialize()
{
  if(_extraAdditionalTerms)
  {
    _extraAdditionalTerms->init(*_simulation->nonSmoothDynamicalSystem()->topology()->dSG(0),
                                *_simulation->nonSmoothDynamicalSystem(),
                                _simulation->eventsManager()->timeDiscretisation());
  }

  initialize_nonsmooth_problems();
  _isInitialized=true;

}

void OneStepIntegrator::update_interaction_output(Interaction& inter, double time, InteractionProperties& interaction_properties)
{
  // - compute interaction output (y) for all levels
  // - swaps in memory

  // This function:
  // - is called during osi->initialize()
  // - may be called by simulation when an interaction is inserted on the fly
  //    for instance:
  //      - contact detection
  //      - inter = ew interaction + link with ds
  //      - simu->osi->initializeWorkVectorsForInteraction(inter)
  //      - simu->osi->update_interaction_output()

  if(_steps > 1)  // Multi--step methods
  {
    // Compute the old Values of Output with stored values in Memory
    for(unsigned int k = 0; k < _steps - 1; k++)
    {
      /** ComputeOutput to fill the Memory
       * We assume the state x is stored in xMemory except for the  initial
       * condition which has not been swap yet.
       */
      //        relation()->LinkDataFromMemory(k);
      for(unsigned int i = 0; i < inter.upperLevelForOutput() + 1; ++i)
      {
        inter.computeOutput(time, i);
        //_yMemory[i]->swap(*_y[i]);
      }
    }
    inter.swapInMemory();
  }
  // Compute a first value for the output
  inter.computeOutput(time,  0);

  // prepare the gradients
  inter.relation()->computeJach(time, inter);
  for(unsigned int i = 0; i < inter.upperLevelForOutput() + 1; ++i)
  {
    inter.computeOutput(time, i);
  }
  inter.swapInMemory();
}

void OneStepIntegrator::_check_and_update_interaction_levels(Interaction& inter)
{

  bool isInitializationNeeded = false;
  if(!(inter.lowerLevelForOutput() <= _levelMinForOutput && inter.upperLevelForOutput()  >= _levelMaxForOutput))
  {
    inter.setLowerLevelForOutput(_levelMinForOutput);
    inter.setUpperLevelForOutput(_levelMaxForOutput);
    isInitializationNeeded = true;
  }

  if(!(inter.lowerLevelForInput() <= _levelMinForInput && inter.upperLevelForInput() >= _levelMaxForInput))
  {
    inter.setLowerLevelForInput(_levelMinForInput);
    inter.setUpperLevelForInput(_levelMaxForInput);
    isInitializationNeeded = true;
  }

  if(isInitializationNeeded)
    inter.reset();
}

void OneStepIntegrator::resetAllNonSmoothParts()
{
  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    _dynamicalSystemsGraph->bundle(*dsi)->resetAllNonSmoothParts();
  }
}

void OneStepIntegrator::resetNonSmoothPart(unsigned int level)
{
  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    _dynamicalSystemsGraph->bundle(*dsi)->resetNonSmoothPart(level);
  }
}

void OneStepIntegrator::updateOutput(double time, unsigned int level)
{
  InteractionsGraph::VIterator ui, uiend;
  InteractionsGraph & indexSet0 = *_simulation->nonSmoothDynamicalSystem()->topology()->indexSet0();
  for(std::tie(ui, uiend) = indexSet0.vertices(); ui != uiend; ++ui)
  {
    if(!checkInteractionOSI(indexSet0, ui)) continue;
    Interaction & inter = *indexSet0.bundle(*ui);
    assert(inter.lowerLevelForOutput() <= level);
    assert(inter.upperLevelForOutput() >= level);
    inter.computeOutput(time, level);
  }
}


void OneStepIntegrator::updateOutput(double time)
{
  for(unsigned int level = _levelMinForOutput; level<=_levelMaxForOutput; ++level)
    updateOutput(time, level);
}


void OneStepIntegrator::updateInput(double time)
{
  for(unsigned int level = _levelMinForInput;
          level < _levelMaxForInput + 1;
       level++)
    updateInput(time, level);
}

void OneStepIntegrator::updateInput(double time, unsigned int level)
{
  // Warning: This reset may be prone to issue with multiple osis.
  // resetNonSmoothPart(level);
  // We compute input using lambda(level).
  InteractionsGraph::VIterator ui, uiend;
  SP::Interaction inter;
  InteractionsGraph & indexSet0 = *_simulation->nonSmoothDynamicalSystem()->topology()->indexSet0();
  for(std::tie(ui, uiend) = indexSet0.vertices(); ui != uiend; ++ui)
  {
    if(!checkInteractionOSI(indexSet0, ui)) continue;
    Interaction & inter = *indexSet0.bundle(*ui);
    assert(inter.lowerLevelForInput() <= level);
    assert(inter.upperLevelForInput() >= level);
    inter.computeInput(time, level);
  }

}


double OneStepIntegrator::computeResiduOutput(double time, SP::InteractionsGraph indexSet)
{
  double residu =0.0;
  THROW_EXCEPTION("OneStepIntegrator::computeResiduOutput not implemented for integrator of type " + std::to_string(_integratorType));
  return residu;
}

double OneStepIntegrator::computeResiduInput(double time, SP::InteractionsGraph indexSet)
{
  double residu =0.0;
  THROW_EXCEPTION("OneStepIntegrator::computeResiduInput not implemented for integrator of type " + std::to_string(_integratorType));
  return residu;
}





void OneStepIntegrator::display()
{
  std::cout << "==== OneStepIntegrator display =====" <<std::endl;
  std::cout << "| _integratorType : " << _integratorType <<std::endl;
  std::cout << "| _sizeMem: " << _sizeMem <<std::endl;
  std::cout << "====================================" <<std::endl;
}
