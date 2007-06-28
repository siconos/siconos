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

/*! \file ControlManager.h
  Tools to control a Model: Sensors and Actuators.
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
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 2.1.0.
 *  \date (Creation) February 08, 2007
 *
 *
 *
 */
class ControlManager
{
protected:

  /** A list of Sensors */
  Sensors allSensors;

  /** A list of Sensors */
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
   * \param the type (class name) of the Sensor
   * \param the TimeDiscretisation of the Sensor
   * \return a pointer to the added Sensor
   */
  Sensor* addSensor(const std::string&, TimeDiscretisation*);

  /** To build and add a new Actuator in the Manager
   * \param the type (class name) of the Actuator
   * \param the TimeDiscretisation of the Actuator
   * \return a pointer to the added Actuator
   */
  Actuator* addActuator(const std::string&, TimeDiscretisation*);

  /** initialize sensor data.
   */
  void initialize();

  /** display the data of the ControlManager on the standard output
   */
  void display() const;

};

#endif
