/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
/*! \file
*/

#ifndef SICONOS_DATAC_H
#define SICONOS_DATAC_H



#include "SiconosKernel.hpp"

enum DATAC_STATUS {DATAC_NULL, DATAC_MODEL, DATAC_INIT, DATAC_FULL };
typedef  std::vector <NonSmoothLaw*> NonSmoothLawSet;
typedef  std::vector <Relation*> RelationsSet;
typedef  std::vector <TimeDiscretisation*> TimesSet;


/** Class to manage a global structure for siconos simulation
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 2.1.0.
 *  \date (Creation) Apr 26, 2004
 *
 *  This class is dedicated to be declared in a global scope and to be use
 *  by the global API C to store informations about siconos simulation.
 *
 *  It allows to describe an API C interface with integer id and not object
 *   pointer.
 *
 */

class DataC
{

protected:

  /** The object is initialized */
  int _Status;

  DynamicalSystemsSet *_DSSet;
  NonSmoothLawSet *_NSLawSet;
  RelationsSet *_RelationsSet;
  InteractionsSet     *_InteractionsSet;
  NonSmoothDynamicalSystem *_NSDS;
  Model *_Model;
  TimesSet *_TimesSet;
  SP::Simulation _Simulation;
  EventsManager  *_EventsManager;

public:

  /** default constructor - no parameters
  *  \param
  *  \param
  */

  DataC();

  /** destructor
  */
  ~DataC();

  // GETTERS/SETTERS

  /** get initialization state
  *  \return int : DATAC_STATUS
  */
  int getStatus();
  void setStatus(int status);

  Model *model();
  void  setModelPtr(Model *ptrModel);

  SP::Simulation simulation();
  void setSimulationPtr(SP::Simulation ptrSimulation);

  NonSmoothDynamicalSystem * nonSmoothDynamicalSystem();
  void setNonSmoothDynamicalSystemPtr(NonSmoothDynamicalSystem * ptr);

  DynamicalSystemsSet * dynamicalSystemsSet();

  InteractionsSet* interactionsSet();

  RelationsSet* relationsSet();

  TimesSet* timesSet();

  NonSmoothLawSet* nonSmoothLawSet();

  EventsManager* eventsManager();
  void setEventsManagerPtr(EventsManager* ptr);

};

#endif // SICONOS_DATAC_H
