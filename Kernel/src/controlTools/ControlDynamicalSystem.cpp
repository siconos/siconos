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
#include "ControlDynamicalSystem.hpp"

using namespace std;

ControlDynamicalSystem::ControlDynamicalSystem(double t0, double T, double h):
  _t0(t0), _T(T), _h(h), _theta(0.5)
{
}

void ControlDynamicalSystem::setTheta(unsigned int newTheta)
{
  _theta = newTheta;
}

void ControlDynamicalSystem::initialize()
{
  _model.reset(new Model(_t0, _T));
  _model->nonSmoothDynamicalSystem()->insertDynamicalSystem(_processDS);
  _processTD.reset(new TimeDiscretisation(_t0, _h));
  _processSimulation.reset(new TimeStepping(_processTD, 0));
  _processSimulation->setName("plant simulation");
  _processIntegrator.reset(new Moreau(_processDS, _theta));
  _processSimulation->insertIntegrator(_processIntegrator);
  _model->initialize(_processSimulation);
}
