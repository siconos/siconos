/* Siconos-Kernel, Copyright INRIA 2005-2012.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
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

void OneStepIntegrator::initialize()
{
  if (_extraAdditionalTerms)
  {
    Model& m = *_simulation->model();
    _extraAdditionalTerms->init(*m.nonSmoothDynamicalSystem()->topology()->dSG(0), m);
  }

  // a subgraph has to be implemented.
  _dynamicalSystemsGraph = _simulation->nonSmoothDynamicalSystem()->topology()->dSG(0);
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

void OneStepIntegrator::display()
{
  std::cout << "==== OneStepIntegrator display =====" <<std::endl;
  std::cout << "| _integratorType : " << _integratorType <<std::endl;
  std::cout << "| _sizeMem: " << _sizeMem <<std::endl;
  std::cout << "====================================" <<std::endl;
}
