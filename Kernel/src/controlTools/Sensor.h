/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2006.
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

/*! \file Sensor.h
  General interface to define a sensor ...
*/

#ifndef Sensor_H
#define Sensor_H

#include "Tools.h"
#include "SiconosAlgebra.h"
#include "EventsManager.h"

class SiconosVector;
class TimeDiscretisation;
class Model;
class Event;

/** The object used to store data in the Sensor. To each Event corresponds a Data */
typedef std::map<Event*, VectorMap>  DataSet;

/** Sensor Base Class
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 2.1.1.
 *  \date (Creation) February 01, 2007
 *
 *  This class is dedicated to data capture. It gives an interface for User who can implement its own Sensor to
 *  clearly define which data he needs to save.
 *  A Sensor is linked to a set of Events, and each time an Event occurs, capture() function is called.
 *  Then to define a Sensor, you need:
 *  - a link to a Model, where from data will be caught,
 *  - a TimeDiscretisation, to set instants of capture,
 *
 *  Moreover, a Sensor is identified thanks to an id and a type (the derived class type indeed).
 *
 *  The TimeDiscretisation object is used to build a set of SensorEvent, which determines when data will be captured.
 *  More precisely, for each tk in the TimeDiscretisation, an Event is added into the EventsManager of the Simulation.
 *  And a simulation->processEvents called leads to SensorEvent->process and thus Sensor->capture() call.
 *
 * The data are saved in a DataSet object named data, a map which associate to each Event another map.
 * This second map links a string, used to identify the data, and a SiconosVector.
 * As an example consider the case where you need to save the state vector x of a DynamicalSystem, then you can define
 * a Data object, with "myDS_X" as an id and yourDS->getXPtr() as the SiconosVector*. For myEvent being an Event where you
 * need to save data, you get:
 *  (data[myEvent])["myDS_X] = model->getNonSmoothDynamicalSystemPtr()->getDynamicalSystemPtr()->getXPtr()
 *
 * See \ref UMSiconosControl for details on how to define its own Sensor.
 *
 *
 */
class Sensor
{
protected:

  /** type of the Sensor */
  std::string type;

  /** id of the Sensor */
  std::string id;

  /** data list - Each vector of data is identified with a string. */
  DataSet data;

  /** The model linked to this sensor */
  Model * model;

  /** A time discretisation scheme */
  TimeDiscretisation *timeDiscretisation;

  /** The set of Events (SensorEvents) linked to this sensor */
  EventsContainer eventsSet;

  /** default constructor
   */
  Sensor();

  /** copy constructor
   * Private => no copy nor pass-by value allowed.
   */
  Sensor(const Sensor&);

public:

  /** Constructor with a TimeDiscretisation.
   * \param a string, the type of the Sensor, which corresponds to the class type
   * \param a TimeDiscretisation*, (linked to a model).
   */
  Sensor(const std::string&, TimeDiscretisation*);

  /** destructor
   */
  virtual ~Sensor();

  /** set id of the Sensor
   *  \param a string
   */
  inline void setId(const std::string& newId)
  {
    id = newId;
  };

  /** get id of the Sensor
   *  \return a string
   */
  inline const std::string getId() const
  {
    return id;
  };

  /** get the type of the Sensor
   *  \return a string
   */
  inline const std::string getType() const
  {
    return type;
  };

  /** get the Model linked to this Sensor
   *  \return a pointer to Model
   */
  inline Model* getModelPtr() const
  {
    return model;
  };

  /** get the TimeDiscretisation linked to this Sensor
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

  /** get all the data saved for this sensor
  *  \return a DataSet
  */
  inline const DataSet getData() const
  {
    return data;
  };

  /** initialize sensor data.
   */
  virtual void initialize();

  /** capture data when the SensorEvent is processed => set data[SensorEvent]=...
   */
  virtual void capture() = 0;

  /** display the data of the Sensor on the standard output
   */
  void display() const;

};

#endif
