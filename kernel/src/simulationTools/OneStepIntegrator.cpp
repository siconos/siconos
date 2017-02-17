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

#include <SiconosConfig.h>
#if defined(SICONOS_STD_FUNCTIONAL) && !defined(SICONOS_USE_BOOST_FOR_CXX11)
#include <functional>
using namespace std::placeholders;
#else
#include <boost/bind.hpp>
#include <boost/weak_ptr.hpp>
#endif



OneStepIntegrator::OneStepIntegrator(const OSI::TYPES& id):
  _integratorType(id), _sizeMem(1)
{
}

void OneStepIntegrator::initialize( Model& m )
{
  if (_extraAdditionalTerms)
  {
    _extraAdditionalTerms->init(*m.nonSmoothDynamicalSystem()->topology()->dSG(0), m);
  }
  // a subgraph has to be implemented.
  _dynamicalSystemsGraph = _simulation->nonSmoothDynamicalSystem()->topology()->dSG(0);
}
void OneStepIntegrator::initializeDynamicalSystem(Model& m, double t, SP::DynamicalSystem ds)
{
  RuntimeException::selfThrow("OneStepIntegrator::initializeDynamicalSystem not implemented for integrator of type " + _integratorType);
}


void OneStepIntegrator::computeInitialNewtonState()
{
  // Default behavior :  do nothing and used the current state as starting state of the Newton iteration
}

double OneStepIntegrator::computeResidu()
{
  RuntimeException::selfThrow("OneStepIntegrator::computeResidu not implemented for integrator of type " + _integratorType);
  return 0.0;
}

void OneStepIntegrator::computeFreeState()
{
  RuntimeException::selfThrow("OneStepIntegrator::computeFreeState not implemented for integrator of type " + _integratorType);
}

void OneStepIntegrator::computeFreeOutput(InteractionsGraph::VDescriptor& vertex_inter, OneStepNSProblem* osnsp)
{
  RuntimeException::selfThrow("OneStepIntegrator::computeFreeOutput not implemented for integrator of type " + _integratorType);
}

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
