/* Siconos-Kernel, Copyright INRIA 2005-2011.
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

/*! \file Sensor.hpp
  \brief General interface to define a sensor
*/

#ifndef Sensor_H
#define Sensor_H

#include "Tools.hpp"
#include "SiconosAlgebra.hpp"
#include "EventsManager.hpp"

class SiconosVector;
#include "TimeDiscretisation.hpp"
class Model;
class Event;

/** A map that links a string to a pointer to SiconosVector. */
typedef std::map<std::string, SP::SiconosVector> VectorMap;

/** An iterator through a map that links a string to a pointer to SiconosVector. */
typedef VectorMap::iterator VectorMapIterator;

/** A const iterator through a map that links a string to a pointer to SiconosVector. */
typedef VectorMap::const_iterator VectorMapConstIterator;

/** The object used to store data in the Sensor. To each Event corresponds a Data */
typedef std::map<SP::Event, VectorMap>  DataSet;

/** Sensor Base Class

   \author SICONOS Development Team - copyright INRIA
   \version 3.0.0.
   \date (Creation) February 01, 2007

   Abstract class, interface to user-defined sensors.

   A Sensor is dedicated to data capture. It gives an interface for
   User who can implement its own Sensor to clearly define which data
   he needs to save.

   A Sensor handles a TimeDiscretisation, which defines the set of all
   instants where the sensor must operate \n (i.e. each times where
   capture() function will be called). An Event, inserted into the
   EventsManager of the Simulation, is linked to this
   TimeDiscretisation.

   Moreover, a Sensor is identified thanks to an id and a type (a
   number associated to the derived class type indeed).

   \section BSensor Construction

   To build a Sensor it is necessary to use the factory. Inputs are a
   number which identify the derived class type and a
   TimeDiscretisation:

   \code

   // Get the registry
   SensorFactory::Registry& regSensor(SensorFactory::Registry::get()) ;
   // Build a Sensor of type "myType" with t as a TimeDiscretisation.
   regSensor.instantiate(myType, t);

   \endcode

   The best way is to use the controlManager:
   \code
   // cm a ControlManager
   cm->addSensor(myType,t);
   // or if cm has already been initialized:
   cm->addAndRecordSensor(myType,t)
   \endcode


   The data are saved in a DataSet object named data, a map which
   associate to each Event another map.  This second map links a
   string, used to identify the data, and a SiconosVector.  As an
   example consider the case where you need to save the state vector x
   of a DynamicalSystem, then you can define a Data object, with
   "myDS_X" as an id and yourDS->x() as the SiconosVector. For
   myEvent being an Event where you need to save data, you get:
   (data[myEvent])["myDS_X"] =
   model->nonSmoothDynamicalSystem()->dynamicalSystem()->x()

   See \ref UMSiconosControl for details on how to define its own Sensor.

 */
class Sensor : public boost::enable_shared_from_this<Sensor>
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(Sensor);


  /** type of the Sensor */
  int _type;

  /** id of the Sensor */
  std::string _id;

  /** data list - Each vector of data is identified with a string. */
  DataSet _data;

  /** The model linked to this sensor */
  SP::Model _model;

  /** A time discretisation scheme */
  SP::TimeDiscretisation _timeDiscretisation;

  /** The event which will linked this sensor to the eventsManager of the simulation */
  SP::Event _eSensor;

  /** default constructor
   */
  Sensor();

  /** copy constructor
   * Private => no copy nor pass-by value allowed.
   */
  Sensor(const Sensor&);

public:

  /** Constructor with a TimeDiscretisation.
   * \param int, the type of the Sensor, which corresponds to the class type
   * \param a TimeDiscretisation*, (linked to a model).
   */
  Sensor(int, SP::TimeDiscretisation);

  /** destructor
   */
  virtual ~Sensor();

  /** set id of the Sensor
   *  \param a string
   */
  inline void setId(const std::string& newId)
  {
    _id = newId;
  };

  /** get id of the Sensor
   *  \return a string
   */
  inline const std::string getId() const
  {
    return _id;
  };

  /** get the type of the Sensor
   *  \return an int
   */
  inline int getType() const
  {
    return _type;
  };

  /** get the Model linked to this Sensor
   *  \return a pointer to Model
   */
  inline SP::Model model() const
  {
    return _model;
  };

  /** get the TimeDiscretisation linked to this Sensor
  *  \return a pointer to TimeDiscretisation.
  */
  inline SP::TimeDiscretisation timeDiscretisation() const
  {
    return _timeDiscretisation;
  };

  /** get the Event associated with this sensor
   *  \return an Event*
   */
  inline SP::Event event() const
  {
    return _eSensor;
  };

  /** get all the data saved for this sensor
   *  \return a DataSet
   */
  inline const DataSet getData() const
  {
    return _data;
  };

  /** initialize sensor data.
   */
  virtual void initialize();

  /** Add the sensor into the simulation EventsManager.
   */
  void recordInSimulation();

  /** capture data when the SensorEvent is processed => set data[SensorEvent]=...
   */
  virtual void capture() = 0;

  /** display the data of the Sensor on the standard output
   */
  void display() const;

};
DEFINE_SPTR(Sensor);
#endif
