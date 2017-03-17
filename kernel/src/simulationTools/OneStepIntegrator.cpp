/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
#include "Model.hpp"
#include "DynamicalSystem.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "ExtraAdditionalTerms.hpp"
#include "Relation.hpp"
#include <SiconosConfig.h>
#if defined(SICONOS_STD_FUNCTIONAL) && !defined(SICONOS_USE_BOOST_FOR_CXX11)
#include <functional>
using namespace std::placeholders;
#else
#include <boost/bind.hpp>
#include <boost/weak_ptr.hpp>
#endif

void OneStepIntegrator::initialize( Model& m )
{
  if (_extraAdditionalTerms)
  {
    _extraAdditionalTerms->init(*m.nonSmoothDynamicalSystem()->topology()->dSG(0), m);
  }
  // a subgraph has to be implemented.
  _dynamicalSystemsGraph = _simulation->nonSmoothDynamicalSystem()->topology()->dSG(0);

  double t0 = _simulation->startingTime();
  
  // 1 - Loop over all dynamical systems
  //  For each ds, allocate/initialize the required buffers in the ds graph
  // (workDS and workMatrices properties, that depend both on osi and ds types)
  // Note FP : what if a DS is associated with more than one osi?
  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
    _dynamicalSystemsGraph->properties(*dsi).workVectors.reset(new VectorOfVectors());
    _dynamicalSystemsGraph->properties(*dsi).workMatrices.reset(new VectorOfMatrices());

    // Initialize memory buffers
    ds->initMemory(getSizeMem());
    // Force dynamical system to its initial state
    ds->resetToInitialState();
    // part of the work which depends on osi and/or DS types.
    initializeDynamicalSystem(m, t0, ds);
  }

  // 2 - Nonsmooth problems : set levels and initialize. Depends on OSI type.
  // Note FP : is it the right place for this initialization??
  initialize_nonsmooth_problems();

  // 3 - Loop over all interactions of index set 0.
  // For each interaction, allocate/initialize buffers in the interaction graph
  // (DSlink property) and connect/fill these buffers with DS buffers.
  // This strongly depends on the DS and OSI types.
  SP::InteractionsGraph indexSet0 = m.nonSmoothDynamicalSystem()->topology()->indexSet0();
  InteractionsGraph::VIterator ui, uiend;
  for (std11::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
  {
    Interaction& inter = *indexSet0->bundle(*ui);
    unsigned int nslawSize = inter.nonSmoothLaw()->size();
    InteractionProperties& interaction_properties = indexSet0->properties(*ui);
    // init block property. Note FP: this should probably be moved
    // to OSNSPb init?
    interaction_properties.block.reset(new SimpleMatrix(nslawSize, nslawSize));

    // Update DSlink : this depends on OSI and relation types.
    initializeInteraction(t0, inter, interaction_properties, *_dynamicalSystemsGraph);

    // Update interaction attributes (output)
    if (_steps > 1) // Multi--step methods
      {
	// Compute the old Values of Output with stored values in Memory
	for (unsigned int k = 0; k < _steps - 1; k++)
	  {
	    /** ComputeOutput to fill the Memory
	     * We assume the state x is stored in xMemory except for the  initial
	     * condition which has not been swap yet.
	     */
	    //        relation()->LinkDataFromMemory(k);
	    for (unsigned int i = 0; i < inter.upperLevelForOutput() + 1; ++i)
	      {
		inter.computeOutput(t0, interaction_properties, i);
		//_yMemory[i]->swap(*_y[i]);
	      }
	  }
	inter.swapInMemory();
	
      }
    // Compute a first value for the output
    inter.computeOutput(t0, interaction_properties, 0);
    
    // prepare the gradients
    inter.relation()->computeJach(t0, inter, interaction_properties);
    for (unsigned int i = 0; i < inter.upperLevelForOutput() + 1; ++i)
      {
      inter.computeOutput(t0, interaction_properties, i);
    }
    inter.swapInMemory();
  }
}


void OneStepIntegrator::_check_and_update_interaction_levels(Interaction& inter)
{
  
  bool isInitializationNeeded = false;
  if (!(inter.lowerLevelForOutput() <= _levelMinForOutput && inter.upperLevelForOutput()  >= _levelMaxForOutput ))
  {
    inter.setLowerLevelForOutput(_levelMinForOutput);
    inter.setUpperLevelForOutput(_levelMaxForOutput);
    isInitializationNeeded = true;
  }

  if (!(inter.lowerLevelForInput() <= _levelMinForInput && inter.upperLevelForInput() >= _levelMaxForInput ))
  {
    inter.setLowerLevelForInput(_levelMinForInput);
    inter.setUpperLevelForInput(_levelMaxForInput);
    isInitializationNeeded = true;
  }

  if (isInitializationNeeded)
    inter.reset();
}

// void OneStepIntegrator::initializeDynamicalSystem(Model& m, double t, SP::DynamicalSystem ds)
// {
//   RuntimeException::selfThrow("OneStepIntegrator::initializeDynamicalSystem not implemented for integrator of type " + _integratorType);
// }


void OneStepIntegrator::resetNonSmoothPart()
{
 DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if (!checkOSI(dsi)) continue;
    _dynamicalSystemsGraph->bundle(*dsi)->resetAllNonSmoothPart();
  }
}

void OneStepIntegrator::resetNonSmoothPart(unsigned int level)
{
 DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if (!checkOSI(dsi)) continue;
    _dynamicalSystemsGraph->bundle(*dsi)->resetNonSmoothPart(level);
  }
}

void OneStepIntegrator::updateOutput(double time)
{
  /** VA. 16/02/2017 This should normally be done only for interaction managed by the osi */
  for (unsigned int level = _levelMinForOutput;
       level < _levelMaxForOutput + 1;
       level++)
    _simulation->nonSmoothDynamicalSystem()->updateOutput(time,level);
}

void OneStepIntegrator::updateInput(double time)
{
  /** VA. 16/02/2017 This should normally be done only for interaction managed by the osi */
  for (unsigned int level = _levelMinForInput;
       level < _levelMaxForInput + 1;
       level++)
    _simulation->nonSmoothDynamicalSystem()->updateInput(time,level);
}

void OneStepIntegrator::updateOutput(double time, unsigned int level)
{
  /** VA. 16/02/2017 This should normally be done only for interaction managed by the osi */
  _simulation->nonSmoothDynamicalSystem()->updateOutput(time,level);
}

void OneStepIntegrator::updateInput(double time, unsigned int level)
{
  /** VA. 16/02/2017 This should normally be done only for interaction managed by the osi */
  _simulation->nonSmoothDynamicalSystem()->updateInput(time,level);
}




void OneStepIntegrator::display()
{
  std::cout << "==== OneStepIntegrator display =====" <<std::endl;
  std::cout << "| _integratorType : " << _integratorType <<std::endl;
  std::cout << "| _sizeMem: " << _sizeMem <<std::endl;
  std::cout << "====================================" <<std::endl;
}
