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

/*! \file InteractionManager.hpp
  \brief Definition of a class that manages dynamic interactions.
*/

#ifndef InteractionManager_h
#define InteractionManager_h

#include "Interaction.hpp"
#include "Model.hpp"
#include "DynamicalSystem.hpp"

#include "SiconosVisitor.hpp"

class InteractionManager
{
public:
  virtual ~InteractionManager() {}

  virtual void updateInteractions(SP::Simulation simulation) {}

  virtual SP::SiconosVisitor getDynamicalSystemsVisitor(SP::Simulation simulation)
    { return SP::SiconosVisitor(); }

  virtual void insertNonSmoothLaw(SP::NonSmoothLaw, int group1, int group2) = 0;
  virtual SP::NonSmoothLaw nonSmoothLaw(int group1, int group2) = 0;

protected:
  // May only be called from updateInteractions()
  void link(SP::Interaction inter,
            SP::DynamicalSystem ds1,
            SP::DynamicalSystem ds2 = SP::DynamicalSystem());

  // May only be called from updateInteractions()
  void unlink(SP::Interaction inter);

  SP::Simulation _simulation;

  // Will be called by simulation before updateInteractions() for
  // link() and unlink() to work.
  void setSimulation(SP::Simulation sim) { _simulation = sim; }

  friend class Simulation;
  friend class TimeStepping;
  friend class EventDriven;
};

#endif /* InteractionManager_h */
