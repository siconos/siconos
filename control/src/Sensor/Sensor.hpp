/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

/*! \file Sensor.hpp
  \brief General interface to define a sensor
*/

#ifndef Sensor_H
#define Sensor_H

#include "RuntimeException.hpp"
#include "SiconosPointers.hpp"

#include "SiconosFwd.hpp"

#include "SiconosAlgebraTypeDef.hpp"

#include "ControlTypeDef.hpp"
#include "SiconosControlFwd.hpp"
/** A map that links a std::string to a pointer to SiconosVector. */
//typedef std::map<std::string, SP::SiconosVector> VectorMap;

/** An iterator through a map that links a std::string to a pointer to SiconosVector. */
//typedef VectorMap::iterator VectorMapIterator;

/** A const iterator through a map that links a std::string to a pointer to SiconosVector. */
//typedef VectorMap::const_iterator VectorMapConstIterator;

/** The object used to store data in the Sensor. To each Event corresponds a Data */
//typedef std::map<SP::Event, VectorMap>  DataSet;

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
   std::string, used to identify the data, and a SiconosVector.  As an
   example consider the case where you need to save the state vector x
   of a DynamicalSystem, then you can define a Data object, with
   "myDS_X" as an id and yourDS->x() as the SiconosVector. For
   myEvent being an Event where you need to save data, you get:
   (data[myEvent])["myDS_X"] =
   model->nonSmoothDynamicalSystem()->dynamicalSystem()->x()

   See users' manual for details on how to define its own Sensor.

 */
class Sensor
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(Sensor);

  /** type of the Sensor */
  unsigned int _type;

  /** id of the Sensor */
  std::string _id;

  /** pointer to the DynamicalSystem we are measuring */
  SP::DynamicalSystem _DS;

  /** pointer to the state of the DynamicalSystem */
  SP::SiconosVector _DSx;

  /** default constructor
   */
  Sensor();

  /** copy constructor
   * Private => no copy nor pass-by value allowed.
   */
  Sensor(const Sensor&);

public:

  /** Constructor with a TimeDiscretisation.
   * \param type the type of the Sensor, which corresponds to the class type.
   * \param ds the SP::DynamicalSystem we observe.
   */
  Sensor(unsigned int type, SP::DynamicalSystem ds);

  /** destructor
   */
  virtual ~Sensor();

  /** set id of the Sensor
   *  \param newId the id of the Sensor
   */
  inline void setId(const std::string& newId)
  {
    _id = newId;
  };

  /** get id of the Sensor
   *  \return a std::string
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

  /** get the DynamicalSystem linked to this Sensor
   *  \return SP::DynamicalSystem
   */
  inline SP::DynamicalSystem getDS() const
  {
    return _DS;
  };

  /** This is derived in child classes if they need to copy the TimeDiscretisation
   * associated with this Sensor
  *  \param td the TimeDiscretisation for this Sensor
  */
  virtual void setTimeDiscretisation(const TimeDiscretisation& td) {};

  /** get all the data saved for this sensor
   *  \return a DataSet
   */
  //  inline const DataSet getData() const
  //  {return _data;};

  /** initialize sensor data.
   * \param nsds the Model
   */
  virtual void initialize(const NonSmoothDynamicalSystem& nsds) {};

  /** capture data when the SensorEvent is processed => set data[SensorEvent]=...
   */
  virtual void capture() = 0;

  /** display the data of the Sensor on the standard output
   */
  void display() const;

};
#endif
