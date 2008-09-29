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

/*! \file Actuator.h
  \brief General interface to define an actuator
*/

#ifndef Actuator_H
#define Actuator_H

#include "SiconosPointers.h"
#include"DynamicalSystemsSet.h"
#include"EventsManager.h"
#include<string>

class TimeDiscretisation;
class Model;
class Event;
class DynamicalSystem;
class Sensor;

/** A set of Sensors */
typedef std::set<Sensor*> Sensors;

/** An iterator through a set of Sensors */
typedef Sensors::iterator SensorsIterator;

/** Return-type for Actuators insertion. */
typedef std::pair<SensorsIterator, bool> SensorsCheckInsert;

TYPEDEF_SPTR(Sensors);

/** Actuators Base Class

   \author SICONOS Development Team - copyright INRIA
   \version 3.0.0.
   \date (Creation) February 09, 2007

   Abstract class, interface to user-defined actuators.

   An Actuator is dedicated to act on parameters of the Model (especially z param. in DynamicalSystem) according to some specific values recorded thanks to sensors. It gives an interface for User who can implement its own Actuator.
   clearly define which data he needs to save.

   An Actuator handles a TimeDiscretisation, which defines the set of all instants where the Actuator must operate \n
   (i.e. each times where actuate() function will be called). An Event, inserted into the EventsManager of the Simulation, is linked to this TimeDiscretisation.

   Moreover, an Actuator is identified thanks to an id and a type (a number associated to the derived class type indeed).

   \section BActuator Construction
   To build an Actuator it is necessary to use the factory. Inputs are a number which identify the derived class type and a TimeDiscretisation:
   \code
   // Get the registry
   ActuatorFactory::Registry& regActuator(ActuatorFactory::Registry::get()) ;
   // Build an Actuator of type "myType" with t as a TimeDiscretisation.
   regActuator.instantiate(myType, t);
   \endcode

   The best way is to use the controlManager:
   \code
   // cm a ControlManager
   cm->addActuator(myType,t);
   // or if cm has already been initialized:
   cm->addAndRecordActuator(myType,t)
   \endcode
*/
class Actuator
{
protected:

  /** type of the Actuator */
  int type;

  /** id of the Actuator */
  std::string id;

  /** Sensors linked to this actuator */
  SP::Sensors  allSensors;

  /** Dynamical Systems list: all the systems on which this actuator may act. */
  SP::DynamicalSystemsSet allDS;

  /** The model linked to this actuator */
  SP::Model  model;

  /** A time discretisation scheme */
  SP::TimeDiscretisation timeDiscretisation;

  /** The event which will linked this actuator to the eventsManager of the simulation */
  SP::Event eActuator;

  /** default constructor
   */
  Actuator();

  /** copy constructor
   * Private => no copy nor pass-by value allowed.
   */
  Actuator(const Actuator&);

public:

  /** Constructor with a TimeDiscretisation.
   * \param a string, the type of the Actuator, which corresponds to the class type
   * \param a TimeDiscretisation*, (linked to a model).
   */
  Actuator(int, SP::TimeDiscretisation);

  /** Constructor with a TimeDiscretisation.
   * \param a string, the type of the Actuator, which corresponds to the class type
   * \param a TimeDiscretisation*, (linked to a model).
   * \param a Sensor linked to this Actuator.
   */
  Actuator(int, SP::TimeDiscretisation, const Sensors&);

  /** destructor
   */
  virtual ~Actuator();

  /** set id of the Actuator
   *  \param a string
   */
  inline void setId(const std::string& newId)
  {
    id = newId;
  };

  /** get id of the Actuator
   *  \return a string
   */
  inline const std::string getId() const
  {
    return id;
  };

  /** get the type of the Actuator (ie class name)
   *  \return a string
   */
  inline const int getType() const
  {
    return type;
  };

  /** get all the Sensors linked to this actuator.
   *  \return a Sensors object (list of Sensor)
   */
  inline const SP::Sensors getSensors() const
  {
    return allSensors;
  };

  /** Ass a set of Sensors into this actuator.
   *  \param a Sensors object (list of Sensor)
   */
  void addSensors(const Sensors&);

  /** add a Sensor in the actuator.
   *  \param a pointer to Sensor
   */
  void addSensorPtr(SP::Sensor);

  /** get all the Dynamical Systems linked to this actuator.
   *  \return a DynamicalSystemsSet.
   */
  inline const SP::DynamicalSystemsSet getDynamicalSystems() const
  {
    return allDS;
  };

  /** Add a set of DynamicalSystem into this actuator.
   *  \param a DynamicalSystemsSet.
   */
  void addDynamicalSystems(const DynamicalSystemsSet&);

  /** add a DynamicalSystem into the actuator.
   *  \param a pointer to DynamicalSystem
   */
  void addDynamicalSystemPtr(SP::DynamicalSystem);

  /** get the Model linked to this Actuator
   *  \return a pointer to Model
   */
  inline SP::Model getModelPtr() const
  {
    return model;
  };

  /** get the TimeDiscretisation linked to this Actuator
  *  \return a pointer to TimeDiscretisation.
  */
  inline SP::TimeDiscretisation getTimeDiscretisationPtr() const
  {
    return timeDiscretisation;
  };

  /** get the Event associated with this actuator
   *  \return an Event*
   */
  inline SP::Event getEventPtr() const
  {
    return eActuator;
  };

  /** initialize actuator data.
   */
  virtual void initialize();

  /** Add the actuator into the simulation EventsManager.
   */
  void recordInSimulation();

  /** capture data when the ActuatorEvent is processed => set data[ActuatorEvent]=...
   */
  virtual void actuate() = 0;

  /** display the data of the Actuator on the standard output
   */
  void display() const;

};

#endif
