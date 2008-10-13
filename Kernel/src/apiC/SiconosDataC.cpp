/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
* Siconos is a program dedicated to the modeling, the simulation and the control
* of non smooth dynamical systems
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

#include "SiconosDataC.h"


DataC::DataC()
{
  _Status = DATAC_NULL;


  _DSSet = new DynamicalSystemsSet();
  _NSLawSet = new NonSmoothLawSet();
  _RelationsSet = new RelationsSet();
  _InteractionsSet = new InteractionsSet();
  _TimesSet = new TimesSet();
  _NSDS = NULL;
  _Model = NULL;
  _Simulation = NULL;

  _EventsManager = NULL;
};

DataC::~DataC()
{
  //  delete _DSSet;
  // delete _InteractionSet;
}

int DataC::getStatus()
{
  return  _Status;
};
void DataC::setStatus(int status)
{
  _Status = status;
};

Model * DataC::getModelPtr()
{
  if (!_Model)
    RuntimeException::selfThrow("ApiC::Model is NULL");

  return _Model;
};

void  DataC::setModelPtr(Model *ptrModel)
{
  _Model = ptrModel;
};

Simulation * DataC::getSimulationPtr()
{
  if (!_Simulation)
    RuntimeException::selfThrow("ApiC::Simulation is NULL");

  return _Simulation;
};
void DataC::setSimulationPtr(Simulation * ptrSimulation)
{
  _Simulation = ptrSimulation;
};

NonSmoothDynamicalSystem * DataC::getNonSmoothDynamicalSystemPtr()
{
  if (!_NSDS)
    RuntimeException::selfThrow("ApiC::NSDS is NULL");

  return _NSDS;
};

void DataC::setNonSmoothDynamicalSystemPtr(NonSmoothDynamicalSystem * ptr)
{
  _NSDS = ptr;
};

TimesSet * DataC::getTimesSetPtr()
{
  if (!_TimesSet)
    RuntimeException::selfThrow("ApiC::TimesSet is NULL");

  return _TimesSet;
};


DynamicalSystemsSet * DataC::getDynamicalSystemsSetPtr()
{
  if (!_DSSet)
    RuntimeException::selfThrow("ApiC::DSSet is NULL");

  return _DSSet;
};

InteractionsSet* DataC::getInteractionsSetPtr()
{
  if (!_InteractionsSet)
    RuntimeException::selfThrow("ApiC::InteractionsSet is NULL");

  return _InteractionsSet;
};


RelationsSet* DataC::getRelationsSetPtr()
{
  if (!_RelationsSet)
    RuntimeException::selfThrow("ApiC::RelationsSet is NULL");

  return _RelationsSet;
};


NonSmoothLawSet* DataC::getNonSmoothLawSetPtr()
{
  if (!_NSLawSet)
    RuntimeException::selfThrow("ApiC::NonSmoothLawSet is NULL");

  return _NSLawSet;
};

EventsManager* DataC::getEventsManagerPtr()
{
  if (!_EventsManager)
    RuntimeException::selfThrow("ApiC::EventsManager is NULL");

  return _EventsManager;
};
void DataC::setEventsManagerPtr(EventsManager* ptr)
{
  _EventsManager = ptr;
}
