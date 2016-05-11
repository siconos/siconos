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

/*! \file ControlSensor.hpp
 * A generic control sensor
*/

#ifndef ControlSensor_H
#define ControlSensor_H

#include "Sensor.hpp"
#include <boost/circular_buffer.hpp>

/** \class ControlSensor
 *  \brief Generic control Sensor to get the output of the system
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.4.0.
 *  \date (Creation) november 09, 2011
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

  virtual void initialize(const Model& m);

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
