/* Siconos version 1.0, Copyright INRIA 2005.
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
#ifndef ONESTEPINTEGRATOR_H
#define ONESTEPINTEGRATOR_H

#include "TimeDiscretisation.h"
#include "DynamicalSystem.h"
#include "SiconosMatrix.h"
#include "SiconosVector.h"
#include "OneStepIntegratorXML.h"
#include "Strategy.h"
#include "SiconosConst.h"
#include "SiconosNumerics.h"
#include "check.h"
#include <iostream>
#include <vector>

const std::string MOREAU_INTEGRATOR = "Moreau";
const std::string ADAMS_INTEGRATOR = "Adams";
const std::string LSODAR_INTEGRATOR = "LSODAR";

class DynamicalSystem;
class TimeDiscretisation;
class OneStepIntegratorXML;

/** \class OneStepIntegrator
 *  \brief It's the generic object which can integre a DynamicalSystem
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) Apr 26, 2004
 *
 * !!! This is a virtual class, interface for some specific integrators !!!
 *
 * At the time, available integrators are: Moreau, LSodar and Adams.
 *
 */
class OneStepIntegrator
{
protected:

  // -- Members --
  /** type of the Integrator */
  std::string integratorType;

  /** The DynamicalSystem to integrate */
  DynamicalSystem *ds;

  /** size of the memory for the integrator */
  int sizeMem;

  /** The time discretisation */
  TimeDiscretisation *timeDiscretisation;

  /** the corresponding XML object */
  OneStepIntegratorXML *integratorXml;

  // -- Default constructor --
  /** \fn OneStepIntegrator()
   *  \brief default constructor
   */
  OneStepIntegrator();

public:

  /** \fn OneStepIntegrator(OneStepIntegratorXML*)
   *  \brief constructor from xml file
   *  \param OneStepIntegratorXML* : the corresponding XML object
   */
  OneStepIntegrator(OneStepIntegratorXML*);

  /** \fn OneStepIntegrator(TimeDiscretisation*, DynamicalSystem* )
   *  \brief constructor from a minimum set of data
   *  \param TimeDiscretisation* : the TimeDiscretisation of the OneStepIntegrator
   *  \param DynamicalSystem* : the DynamicalSystem to be integrated
   */
  OneStepIntegrator(TimeDiscretisation*, DynamicalSystem*);


  /** \fn ~OneStepIntegrator()
   *  \brief destructor
   */
  virtual ~OneStepIntegrator();

  // --- GETTERS/SETTERS ---

  /** \fn inline const string getType() const
   *  \brief get the type of the OneStepIntegrator
   *  \return string : the type of the OneStepIntegrator
   */
  inline const std::string getType() const
  {
    return integratorType;
  }

  /** \fn inline void setType(const string&)
   *  \brief set the type of the OneStepIntegrator
   *  \return string : the type of the OneStepIntegrator
   */
  inline void setType(const std::string& newType)
  {
    integratorType = newType;
  }

  /** \fn DynamicalSystem* getDynamicalSystemPtr()
    *  \brief get the DynamicalSystem associated with the Integrator
    *  \return a DS
    */
  inline DynamicalSystem* getDynamicalSystemPtr() const
  {
    return ds;
  };

  /** \fn void setDynamicalSystemPtr(DynamicalSystem*)
   *  \brief set the DynamicalSystem of this Integrator
   *  \param a pointer to a DynamicalSystem
   */
  inline void setDynamicalSystemPtr(DynamicalSystem* newDs)
  {
    ds = newDs;
  };

  /** \fn const unsigned int getSizeMem() const
   *  \brief get sizeMem value
   *  \return an unsigned int
   */
  inline const unsigned int getSizeMem() const
  {
    return sizeMem;
  };

  /** \fn void setSizeMem(const unsigned int&)
   *  \brief set sizeMem
   *  \param an unsigned int
   */
  inline void setSizeMem(const unsigned int& newValue)
  {
    sizeMem = newValue;
  };

  /** \fn TimeDiscretisation* getTimeDiscretisationPtr()
   *  \brief get the TimeDiscretisation pointer
   *  \return a pointer on the TimeDiscretisation
   */
  inline TimeDiscretisation* getTimeDiscretisationPtr() const
  {
    return timeDiscretisation;
  };

  /** \fn void setTimeDiscretisationPtr(TimeDiscretisation*)
   *  \brief set the TimeDiscretisation pointer
   *  \param a pointer on a TimeDiscretisation
   */
  inline void setTimeDiscretisationPtr(TimeDiscretisation* td)
  {
    timeDiscretisation = td;
  };

  /** \fn inline OneStepIntegratorXML* getOneStepIntegratorXMLPtr()
   *  \brief get the OneStepIntegratorXML of the OneStepIntegrator
   *  \return a pointer on the OneStepIntegratorXML of the OneStepIntegrator
   */
  inline OneStepIntegratorXML* getOneStepIntegratorXMLPtr() const
  {
    return integratorXml;
  }

  /** \fn inline void setOneStepIntegratorXMLPtr(OneStepIntegratorXML* integratorXml)
   *  \brief set the OneStepIntegratorXML of the OneStepIntegrator
   *  \param OneStepIntegratorXML* : the pointer to set the OneStepIntegratorXML
   */
  inline void setOneStepIntegratorXMLPtr(OneStepIntegratorXML* newIntegratorXml)
  {
    integratorXml = newIntegratorXml;
  }


  // --- OTHERS ... ---

  /** \fn void initialize()
   *  \brief initialise the integrator
   */
  virtual void initialize() = 0;

  /** \fn void nextStep()
   *  \brief prepares the DynamicalSystem for the next time step, push x and xDot in Memory
   *  \param to be defined
   *  \exception to be defined
   *  \return void
   */
  virtual void nextStep();

  /** \fn void computeFreeState()
   *  \brief set r = 0 and integrates the DS
   *  \exception RuntimeException if the integration method for this type of DS is not implemented
   *  \return void
   */
  virtual void computeFreeState() = 0;

  /** \fn void integrate()
   *  \brief not implemented
   */
  virtual void integrate() = 0;

  /** \fn void updateState()
   *  \brief update the state of the DynamicalSystem attached to this Integrator
   */
  virtual void updateState() = 0;

  /** \fn void display()
   *  \brief print the data to the screen
   */
  virtual  void display() const;

  /** \fn void saveIntegratorToXML()
   *  \brief copy the data of the OneStepNSProblem to the XML tree
   *  \exception RuntimeException
   */
  virtual void saveIntegratorToXML();

};

#endif // ONESTEPINTEGRATOR_H
