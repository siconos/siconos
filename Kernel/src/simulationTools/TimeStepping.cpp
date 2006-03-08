/* Siconos-Kernel version 1.1.2, Copyright INRIA 2005-2006.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/

#include "TimeStepping.h"
using namespace std;

// --- Default constructor ---
TimeStepping::TimeStepping(Model * newModel): Strategy(newModel)
{
  strategyType = "TimeStepping";
}

// --- From Model ---
TimeStepping::TimeStepping(Model& newModel): Strategy(newModel)
{
  strategyType = "TimeStepping";
}

// --- XML constructor ---
TimeStepping::TimeStepping(StrategyXML* strxml, Model *newModel): Strategy(strxml, newModel)
{
  strategyType = "TimeStepping";
}

// --- Destructor ---
TimeStepping::~TimeStepping()
{}

TimeStepping* TimeStepping::convert(Strategy *str)
{
  TimeStepping* ts = dynamic_cast<TimeStepping*>(str);
  return ts;
}

void TimeStepping::initialize()
{
  // initialization of the OneStepIntegrators
  for (unsigned int i = 0; i < integratorVector.size(); i++)
    integratorVector[i]->initialize();
  // initialization of  OneStepNonSmoothProblem
  if (nsProblem != NULL)
    nsProblem->initialize();
}

void TimeStepping::run()
{
  // Current Step
  unsigned int k = timeDiscretisation->getK();
  // Number of time steps
  unsigned int nSteps = timeDiscretisation->getNSteps();
  while (k < nSteps)
  {
    // transfer of state i+1 into state i and time incrementation
    nextStep();
    // update current time step
    k = timeDiscretisation->getK();

    computeOneStep();

  }
}

// compute simulation between ti and ti+1, ti+1 being currentTime (?)
// Initial DS/interaction state is given by memory vectors (check that?)
// and final state is the one saved in DS/Interaction at the end of this function
void TimeStepping::computeOneStep()
{
  // solve ...
  computeFreeState();
  computeOneStepNSProblem();
  // update
  update();
}

