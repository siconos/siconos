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

/*! \file ControlManager.h
  \brief Tools to control a Model: Sensors and Actuators.
*/

#ifndef ControlManager_H
#define ControlManager_H

#include "Actuator.h"

class Actuator;
class Sensor;
class Model;
class TimeDiscretisation;

/** A set of Actuators */
typedef std::set<Actuator*> Actuators;

/** An iterator through a set of Actuators */
typedef Actuators::iterator ActuatorsIterator;

/** Return-type for Actuators insertion. */
typedef std::pair<ActuatorsIterator, bool> ActuatorsCheckInsert;

/** ControlManager Class: tools to control a Model (Sensors, Actuators ...)

    \author SICONOS Development Team - copyright INRIA
    \version 3.0.0.
    \date (Creation) February 08, 2007

    This class is used to handle all the sensors and actuators declared by the user and to
    schedule them into the simulation.

    A ControlManager has:
    - a list of sensors
    - a list of actuators
    - a link to an existing model

    The usual way to define control over a system is as follows:
    - declare a ControlManager and associate it with a Model
    - add some sensors and actuators into the ControlManager
    - initialize the ControlManager (which will result in the recording of all actuators and sensors into the
    list of events processed during the simulation)
    - optionally add some new sensor/actuator at any time but with a specific function: addAndRecord(...).
    A call to this function results in the creation of a Sensor/Actuator and in the insertion of the corresponding event
    into the simulation eventsManager.


*/
class ControlManager
{
protected:

  /** A list of Sensors */
  Sensors allSensors;

  /** A list of Actuators */
  Actuators allActuators;

  /** The model linked to this ControlManager */
  Model * model;

  /** default constructor
   */
  ControlManager();

  /** copy constructor
   * Private => no copy nor pass-by value allowed.
   */
  ControlManager(const ControlManager&);

public:

  /** Constructor with a Model, to which control will be applied.
   * \param a Model*
   */
  ControlManager(Model*);

  /** destructor
   */
  virtual ~ControlManager();

  /** get the Model linked to this ControlManager
   *  \return a pointer to Model
   */
  inline Model* getModelPtr() const
  {
    return model;
  };

  /** get the list of Sensors associated to this manager.
   *  \return a Sensors object.
   */
  inline const Sensors getSensors() const
  {
    return allSensors ;
  };

  /** get the list of Actuators associated to this manager.
   *  \return a Actuators object.
   */
  inline const Actuators getActuators() const
  {
    return allActuators ;
  };

  /** To build and add a new Sensor in the Manager
   * \param the type (int) of the Sensor
   * \param the TimeDiscretisation of the Sensor
   * \return a pointer to the added Sensor
   */
  Sensor* addSensor(int, TimeDiscretisation*);

  /** To build, add, initialize a new Sensor in the Manager and record it in the simulation
   * This function is only useful to add a new Sensor after the initialization of the manager
   * else call addSensor()
   * \param the type (int) of the Sensor
   * \param the TimeDiscretisation of the Sensor
   * \return a pointer to the added Sensor
   */
  Sensor* addAndRecordSensor(int, TimeDiscretisation*);

  /** To build and add a new Actuator in the Manager
   * \param the type (int) of the Actuator
   * \param the TimeDiscretisation of the Actuator
   * \return a pointer to the added Actuator
   */
  Actuator* addActuator(int, TimeDiscretisation*);

  /** To build, add, initialize a new Actuator in the Manager and record it in the simulation
   * This function is only useful to add a new Actuator after the initialization of the manager
   * else call addActuator()
   * \param the type (int) of the Actuator
   * \param the TimeDiscretisation of the Actuator
   * \return a pointer to the added Actuator
   */
  Actuator* addAndRecordActuator(int, TimeDiscretisation*);



  /** initialize sensor data.
   */
  void initialize();

  /** display the data of the ControlManager on the standard output
   */
  void display() const;

};

#endif
