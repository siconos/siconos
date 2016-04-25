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

/*! \file ControlManager.hpp
  \brief Tools to provide control in a Simulation: Sensors, Observer and Actuators.
*/

#ifndef ControlManager_H
#define ControlManager_H

#include "Actuator.hpp"


#include <set>

/** A set of Sensors */
typedef std::set<SP::Sensor> Sensors;

/** A set of Actuators */
typedef std::set<SP::Actuator> Actuators;

/** A set of Observers */
typedef std::set<SP::Observer> Observers;

/** An iterator through a set of Sensors */
typedef Sensors::iterator SensorsIterator;

/** An iterator through a set of Actuators */
typedef Actuators::iterator ActuatorsIterator;

/** An iterator through a set of Observers */
typedef Observers::iterator ObserversIterator;


/** ControlManager Class: tools to provide control in a Simulation (Sensors, Actuators, Observers)

    \author SICONOS Development Team - copyright INRIA
    \version 3.0.0.
    \date (Creation) February 08, 2007

    This class is used to handle all the sensors and actuators declared by the user and to
    schedule them into the simulation.

    A ControlManager has:
    - a list of Sensor
    - a list of Actuator
    - a link to an existing Simulation

    The usual way to define control over a system is as follows:
    - declare a ControlManager and associate it with a Simulation
    - add some sensors and actuators into the ControlManager
    - initialize the ControlManager (which will result in the recording of all actuators and sensors into the
    list of events processed during the simulation)
    - optionally add some new sensor/actuator at any time but with a specific function: addAndRecord(...).
    A call to this function results in the creation of a Sensor/Actuator and in the insertion of the corresponding event
    into the simulation eventsManager.
*/

class ControlManager
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(ControlManager);

protected:
  /** A list of Sensors */
  Sensors _allSensors;

  /** A list of Actuators */
  Actuators _allActuators;

  /** A list of Observers */
  Observers _allObservers;

  /** The simulation linked to this ControlManager */
  SP::Simulation _sim;

  /** default constructor
   */
  ControlManager(): _sim(SP::Simulation()) {};

  /** Create associated Event and give the opportunity to get the TimeDiscretisation
   * \param s a Sensor
   * \param td a TimeDiscretisation asociated with this Sensor
   */
  void linkSensorSimulation(SP::Sensor s, SP::TimeDiscretisation td);

  /** Create associated Event and give the opportunity to get the TimeDiscretisation
   * \param act a Sensor
   * \param td a TimeDiscretisation asociated with this Sensor
   */
  void linkActuatorSimulation(SP::Actuator act, SP::TimeDiscretisation td);

  /** Create associated Event and give the opportunity to get the TimeDiscretisation
   * \param obs a Sensor
   * \param td a TimeDiscretisation asociated with this Sensor
   */
  void linkObserverSimulation(SP::Observer obs, SP::TimeDiscretisation td);

public:

  /** Constructor with a Simulation, to which control will be applied.
   * \param sim the Simulation
   */
  ControlManager(SP::Simulation sim);

  /** destructor
   */
  virtual ~ControlManager();

  /** get the Simulation linked to this ControlManager
   *  \return a SP::Simulation
   */
  inline SP::Simulation simulation() const
  {
    return _sim;
  };

  /** get the list of Sensors associated to this manager.
   *  \return a Sensors object.
   */
  inline const Sensors getSensors() const
  {
    return _allSensors ;
  };

  /** get the list of Actuators associated to this manager.
   *  \return a Actuators object.
   */
  inline const Actuators getActuators() const
  {
    return _allActuators ;
  };

  /** get the list of Observers associated to this manager.
   *  \return a Observers object.
   */
  inline const Observers getObservers() const
  {
    return _allObservers ;
  };

  /** To build and add a new Sensor in the Manager
   * \param name the type of the Sensor
   * \param td the SP::TimeDiscretisation of the Sensor
   * \param ds the DynamicalSystem used in the Sensor
   * \return a SP::Sensor to the added Sensor
   */
  SP::Sensor addSensor(int name, SP::TimeDiscretisation td, SP::DynamicalSystem ds);

  /** To build, add, initialize a new Sensor in the Manager and record
   * it in the simulation This function is only useful to add a new
   * Sensor after the initialization of the manager else call
   * addSensor()
   * \param name the type (int) of the Sensor
   * \param td the SP::TimeDiscretisation of the Sensor
   * \param ds the DynamicalSystem used in the Sensor
   * \param m the model
   * \return a SP::Sensor to the added Sensor
   */
  SP::Sensor addAndRecordSensor(int name, SP::TimeDiscretisation td, SP::DynamicalSystem ds, const Model& m);

  /** Add an existing Sensor to the Manager
   * \param s a SP::Sensor to the Sensor we want to add
   * \param td the TimeDiscretisation used for the associated Event
   */
  void addSensorPtr(SP::Sensor s, SP::TimeDiscretisation td);

  /** To add, initialize an existing Sensor in the manager and record
   * it in the simulation This function is only useful to add a new
   * Sensor after the initialization of the manager else call
   * addSensor()
   * \param s a SP::Sensor to the Sensor we want to add
   * \param td the TimeDiscretisation used for the associated Event
   */
  void addAndRecordSensorPtr(SP::Sensor s, SP::TimeDiscretisation td, const Model& m);

  /** To build and add a new Actuator in the Manager
   * \param name the type of the Actuator
   * \param td the SP::TimeDiscretisation of the Actuator
   * \param sensor the ControlSensor used to feed the Actuator
   * \return a SP::ACtuator to the added Actuator
   */
  SP::Actuator addActuator(int name, SP::TimeDiscretisation td, SP::ControlSensor sensor);

  /** To build, add, initialize a new Actuator in the manager and
   * record it in the simulation This function is only useful to add a
   * new Actuator after the initialization of the manager else call
   * addActuator()
   * \param name the type of the Actuator
   * \param t the SP::TimeDiscretisation of the Actuator
   * \param sensor the ControlSensor used to feed the Actuator
   * \param m the model
   * \return a SP::Actuator to the added Actuator
   */
  SP::Actuator addAndRecordActuator(int name, SP::TimeDiscretisation t, SP::ControlSensor sensor, const Model& m);

  /** Add an existing Actuator to the manager
   * \param act a SP::Actuator to the Actuator we want to add
   * \param td the TimeDiscretisation used for the associated Event
   */
  void addActuatorPtr(SP::Actuator act, SP::TimeDiscretisation td);

  /** To add, initialize an existing Actuator in the manager and record
   * it in the simulation This function is only useful to add a new
   * Actuator after the initialization of the manager otherwise call
   * addActuator()
   * \param act a SP::Actuator to the Actuator we want to add
   * \param td the TimeDiscretisation used for the associated Event
   */
  void addAndRecordActuatorPtr(SP::Actuator act, SP::TimeDiscretisation td, const Model& m);

  /** To build and add a new Observer in the Manager
   * \param name the type of the Observer
   * \param td the SP::TimeDiscretisation of the Observer
   * \param sensor the ControlSensor feeding the Observer
   * \param xHat0 the initial guess for the state
   * \return a SP::ACtuator to the added Observer
   */
  SP::Observer addObserver(int name, SP::TimeDiscretisation td, SP::ControlSensor sensor, const SiconosVector& xHat0);

  /** To build, add, initialize a new Observer in the manager and
   * record it in the simulation This function is only useful to add a
   * new Observer after the initialization of the manager else call
   * addObserver()
   * \param name the type of the Observer
   * \param td the SP::TimeDiscretisation of the Observer
   * \param sensor the ControlSensor feeding the Observer
   * \param xHat0 the initial guess for the state
   * \return a SP::Observer to the added Observer
   */
  SP::Observer addAndRecordObserver(int name, SP::TimeDiscretisation td, SP::ControlSensor sensor, const SiconosVector& xHat0, const Model& m);

  /** Add an existing Observer to the manager
   * \param obs a SP::Observer to the Observer we want to add
   * \param td the TimeDiscretisation used for the associated Event
   */
  void addObserverPtr(SP::Observer obs, SP::TimeDiscretisation td);

  /** To add, initialize an existing Observer in the manager and record
   * it in the simulation This function is only useful to add a new
   * Observer after the initialization of the manager otherwise call
   * addObserver()
   * \param obs a SP::Observer to the Observer we want to add
   * \param td the TimeDiscretisation used for the associated Event
   */
  void addAndRecordObserverPtr(SP::Observer obs, SP::TimeDiscretisation td, const Model& m);


  /** initialize all Sensors, Observers and Actuators.
   * \param m the Model
   */
  void initialize(const Model& m);

  /** display the data of the ControlManager on the standard output
   */
  void display() const;


};
#endif
