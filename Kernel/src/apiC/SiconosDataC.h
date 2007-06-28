/* Siconos-Kernel version 2.1.0, Copyright INRIA 2005-2006.
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


#include <math.h>
#include <stdio.h>

#include <iostream>
#include "Model.h"
#include "DynamicalSystem.h"
#include "DynamicalSystemsSet.h"
#include "InteractionsSet.h"

using namespace std;

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
  int _initOk;

  Model *_Model;
  Simulation * _Simulation;
  NonSmoothDynamicalSystem *_NSDS;
  TimeDiscretisation *_Time;
  DynamicalSystemsSet *_SetDS;
  CheckInsertDS       *_CheckDS;
  InteractionsSet     *_InteractionsSet;

public:

  /** default constructor
  *  \param
  *  \param
  */
  DataC() {};


  /** destructor
  */
  ~DataC() {};

  // GETTERS/SETTERS

  /** get initialization state
  *  \return int : 0 the object is not initialized, 1 the object is not initialized
  */
  inline int getInitialized()
  {
    return _initOk;
  };

  inline Model * getModelPtr()
  {
    return _Model;
  };
  inline void  setModelPtr(Model *ptrModel)
  {
    _Model = ptrModel;
  };

  inline Simulation * getSimulationPtr()
  {
    return _Simulation;
  };
  inline void setSimulationPtr(Simulation * ptrSimulation)
  {
    _Simulation = ptrSimulation;
  };

  inline NonSmoothDynamicalSystem * getNonSmoothDynamicalSystemPtr()
  {
    return _NSDS;
  };
  inline void setNonSmoothDynamicalSystemPtr(NonSmoothDynamicalSystem * ptr)
  {
    _NSDS = ptr;
  };

  inline TimeDiscretisation * getTimeDiscretisationPtr()
  {
    return _Time;
  };
  inline void setTimeDiscretisationPtr(TimeDiscretisation *ptr)
  {
    _Time = ptr;
  };

  inline DynamicalSystemsSet * getDynamicalSystemsSetPtr()
  {
    return _SetDS;
  };
  inline  void setDynamicalSystemsSetPtr(DynamicalSystemsSet *ptr)
  {
    _SetDS = ptr;
  };

  inline CheckInsertDS* getCheckInsertDStPtr()
  {
    return _CheckDS;
  };
  inline void setCheckInsertDStPtr(CheckInsertDS* ptr)
  {
    _CheckDS = ptr;
  };

  inline InteractionsSet* getInteractionsSetPtr()
  {
    return _InteractionsSet;
  };
  inline void setInteractionsSetPtr(InteractionsSet* ptr)
  {
    _InteractionsSet = ptr;
  };

};

#endif // SICONOS_DATAC_H
