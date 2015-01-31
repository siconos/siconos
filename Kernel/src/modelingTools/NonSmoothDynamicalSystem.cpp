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
#include "NonSmoothDynamicalSystem.hpp"
#include "Topology.hpp"
#include "Interaction.hpp"
#include "LagrangianLinearTIDS.hpp"
#include "FirstOrderLinearTIDS.hpp"
#include "Relation.hpp"

using namespace RELATION;

// --- CONSTRUCTORS/DESTRUCTOR ---

// Default constructor
NonSmoothDynamicalSystem::NonSmoothDynamicalSystem(): _BVP(false), _mIsLinear(true)
{
  // === Builds an empty topology ===
  _topology.reset(new Topology());
};


NonSmoothDynamicalSystem::~NonSmoothDynamicalSystem()
{
  clear();
}

// === DynamicalSystems management ===

void NonSmoothDynamicalSystem::display() const
{
  std::cout << " ===== Non Smooth Dynamical System display ===== " <<std::endl;
  std::cout << "---> isBVP = " << _BVP <<std::endl;
  dynamicalSystems()->begin();
  _topology->indexSet0()->display();
  std::cout << "===================================================" <<std::endl;
}

#include <limits>
double NonSmoothDynamicalSystem::nsdsConvergenceIndicator()
{
  // calculate the max value of all DS convergence indicators
  double convergenceIndicator = -std::numeric_limits<double>::infinity();
  double dsIndic ;
  DynamicalSystemsGraph::VIterator vi;
  for (vi = dynamicalSystems()->begin(); vi != dynamicalSystems()->end(); ++vi)
  {
    dsIndic = dynamicalSystems()->bundle(*vi)->dsConvergenceIndicator();
    if (dsIndic > convergenceIndicator) convergenceIndicator = dsIndic;
  }
  return(convergenceIndicator);
}

void NonSmoothDynamicalSystem::link(SP::Interaction inter, SP::DynamicalSystem ds1, SP::DynamicalSystem ds2)
{
  _mIsLinear = ((inter)->relation()->isLinear() && _mIsLinear);
  _topology->link(inter, ds1, ds2);
};


void NonSmoothDynamicalSystem::clear()
{
  _topology->clear();
}

void NonSmoothDynamicalSystem::setSymmetric(bool val)
{
  _topology->setSymmetric(val);
}

