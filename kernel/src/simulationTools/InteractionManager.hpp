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
#include "NSLawMatrix.hpp"

class InteractionManager
{
public:
  InteractionManager() : _nslaws(1) {}
  virtual ~InteractionManager() {}

  /** Called by Simulation prior to advancing to the next Event */
  virtual void updateInteractions(SP::Simulation simulation) {}

  /** Specify a non-smooth law to use for a given combination of
   *  interaction groups.
   * \param nsl a SP::NonSmoothLaw
   * \param group1 first group
   * \param group2 second group */
  virtual void insertNonSmoothLaw(SP::NonSmoothLaw nslaw,
                                  long unsigned int group1,
                                  long unsigned int group2);

  /** Retrieve a non-smooth law to use for a given combination of
   *  interaction groups.
   * \return nsl a SP::NonSmoothLaw
   * \param group1 first group
   * \param group2 second group */
  virtual SP::NonSmoothLaw nonSmoothLaw(long unsigned int group1,
                                        long unsigned int group2)
    { return _nslaws(group1, group2); }

protected:

  /** nslaws */
  NSLawMatrix _nslaws;

  friend class Simulation;
  friend class TimeStepping;
  friend class EventDriven;

  ACCEPT_SERIALIZATION(InteractionManager);
};

#endif /* InteractionManager_h */
