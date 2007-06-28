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

/*! \file Actuator.h
  General interface to define an actuator ...
*/

#ifndef Actuator_H
#define Actuator_H

#include"DynamicalSystemsSet.h"
#include"EventsManager.h"
#include<string>

class SiconosVector;
class TimeDiscretisation;
class Model;
class Event;
class DynamicalSystemsSet;
class DynamicalSystem;
class Sensor;

/** A set of Sensors */
typedef std::set<Sensor*> Sensors;

/** An iterator through a set of Sensors */
typedef Sensors::iterator SensorsIterator;

/** Return-type for Actuators insertion. */
typedef std::pair<SensorsIterator, bool> SensorsCheckInsert;

/** Actuator Base Class
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 2.1.0.
 *  \date (Creation) February 09, 2007
 *
 * See \ref UMSiconosControl for details on how to define its own Actuator.
 *
 *
 */
class Actuator
{
protected:

  /** type of the Actuator */
  std::string type;

  /** id of the Actuator */
  std::string id;

  /** data list - Each vector of data is identified with a string. */
  Sensors allSensors;

  /** Dynamical Systems list: all the systems on which this actuator may act. */
  DynamicalSystemsSet allDS;

  /** The model linked to this sensor */
  Model * model;

  /** A time discretisation scheme */
  TimeDiscretisation *timeDiscretisation;

  /** The set of Events (ActuatorEvents) linked to this sensor */
  EventsContainer eventsSet;

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
  Actuator(const std::string&, TimeDiscretisation*);

  /** Constructor with a TimeDiscretisation.
   * \param a string, the type of the Actuator, which corresponds to the class type
   * \param a TimeDiscretisation*, (linked to a model).
   * \param a list of Sensors linked to this Actuator.
   */
  Actuator(const std::string&, TimeDiscretisation*, const Sensors&);

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
  inline const std::string getType() const
  {
    return type;
  };

  /** get all the Sensors linked to this actuator.
   *  \return a Sensors object (list of Sensor)
   */
  inline const Sensors getSensors() const
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
  void addSensorPtr(Sensor*);

  /** get all the Dynamical Systems linked to this actuator.
   *  \return a DynamicalSystemsSet.
   */
  inline const DynamicalSystemsSet getDynamicalSystems() const
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
  void addDynamicalSystemPtr(DynamicalSystem*);

  /** get the Model linked to this Actuator
   *  \return a pointer to Model
   */
  inline Model* getModelPtr() const
  {
    return model;
  };

  /** get the TimeDiscretisation linked to this Actuator
  *  \return a pointer to TimeDiscretisation.
  */
  inline TimeDiscretisation* getTimeDiscretisationPtr() const
  {
    return timeDiscretisation;
  };

  /** get the list of Events associated to this sensor
  *  \return an EventsContainer.
  */
  inline const EventsContainer getEvents() const
  {
    return eventsSet ;
  };

  /** initialize sensor data.
   */
  virtual void initialize();

  /** capture data when the ActuatorEvent is processed => set data[ActuatorEvent]=...
   */
  virtual void actuate() = 0;

  /** display the data of the Actuator on the standard output
   */
  void display() const;

};

#endif
