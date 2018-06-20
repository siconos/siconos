/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

/*! \file ControlSensor.hpp
 * A generic control sensor
*/

#ifndef ControlSensor_H
#define ControlSensor_H

#include "Sensor.hpp"
#include <boost/circular_buffer.hpp>

/** \class ControlSensor
 *  \brief Generic control Sensor to get the output of the system
 *
 * A generic control sensor
 *
 */
class ControlSensor : public Sensor
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(ControlSensor);

protected:
  /** A vector for the current value of the output */
  SP::SiconosVector _storedY;

  /** delay between the measurement on the DynamicalSystem and the avaibility of the value */
  double _delay;

  /** A buffer to store the value of \f$y_k\f$ if there is a delay */
  boost::circular_buffer<SP::SiconosVector> _bufferY;

  /** Default constructor
   */
  ControlSensor() {};

  /** Simple Constructor
   * \param type the type of the Sensor
   * \param ds the SP::DynamicalSystem it observes
   * \param delay the delay between the measurement and the avaibility of the data
   */
  ControlSensor(unsigned int type, SP::DynamicalSystem ds, double delay = 0):
    Sensor(type, ds), _delay(delay) {}

public:

  virtual void initialize(const NonSmoothDynamicalSystem& nsds);

  /** Get the dimension of the output
   * \return an unsigned int
   */
  unsigned int getYDim() const ;

  /** Get a pointer to the output
   * \return SP::SiconosVector to the output
   */
  inline const SiconosVector& y() const
  {
    if (_delay == 0)
      return *_storedY;
    else
      return *_bufferY.front();
  };

  inline SP::SiconosVector yTk() const
  {
    return _storedY;
  }

  /** capture data when the SensorEvent is processed => set data[SensorEvent]=...
   */
  // Note: This is redundant to abstract definition in Sensor.hpp,
  // just to help Doxygen's imperfect detection of abstract classes.
  virtual void capture() = 0;
};
#endif
