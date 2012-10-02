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

/*! \file SampledPIDActuator.hpp
  \brief General interface to define an actuator
*/

#ifndef SampledPIDActuator_H
#define SampledPIDActuator_H

#include "SiconosKernel.hpp"
#include <boost/circular_buffer.hpp>

class SampledPIDActuator : public Actuator
{
private:
  /** default constructor */
  SampledPIDActuator() {};

  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(SampledPIDActuator);

  /** error vector */
  std11::shared_ptr<boost::circular_buffer<double> > _err;

  /** reference we are tracking */
  double _ref;

  /** control variable */
  SP::SiconosVector _u;

  /** vector of gains */
  SP::SiconosVector _K;

  /** the sensor that feed the controller */
  SP::ControlSensor _sensor;

  /** boolean to determined if the controller has been correctly initialized */
  bool _initDone;

  /** current \f$ \Delta t\f$ (or timeStep) */
  double _curDeltaT;

public:

  /** Constructor with a TimeDiscretisation.
   * \param t the SP::TimeDiscretisation (/!\ it should not be used elsewhere !).
   * \param ds the SP::DynamicalSystem we are controlling
   */
  SampledPIDActuator(SP::TimeDiscretisation t, SP::DynamicalSystem ds);

  /** Constructor with a TimeDiscretisation.
   * \param t a SP::TimeDiscretisation (/!\ it should not be used elsewhere !).
   * \param ds the SP::DynamicalSystem we are controlling
   * \param sensorList set of Sensor linked to this SampledPIDActuator.
   */
  SampledPIDActuator(SP::TimeDiscretisation t, SP::DynamicalSystem ds, const Sensors& sensorList);

  /** destructor
   */
  virtual ~SampledPIDActuator();

  /** initialize actuator data.
   * \param m a SP::Model
   */
  void initialize(SP::Model m);

  /** Compute the new control law at each event
   * Here we are using the following formula:
   * \f$ u_k = u_{k-1} + c_1 e_k + c_2 e_{k-1} + c_3 e_{k-2} \f$ , where
   * \f{array} c_1 &= K_P - \frac{K_D}{\Delta t} + K_I \Delta t \\
   * c_2 &= -1 - \frac{2K_D}{\Delta t} \\
   * c_3 &= \frac{K_D}{\Delta t} \\
   * \f{array}
   */
  void actuate();

  /** Set the value of _K to newValue
   * * \param newValue SiconosVector \f$ [K_P, K_I, K_D] \f$
   */
  void setK(const SiconosVector& newValue);

  /** Set _K to pointer newPtr
   * \param newPtr SP::SiconosVector f$ [K_P, K_I, K_D] \f$
   */
  void setKPtr(SP::SiconosVector newPtr);

  /** Set the value of _ref to newValue
   * \param newValue
   */
  void inline setRef(const double newValue)
  {
    _ref = newValue;
  }
};
DEFINE_SPTR(SampledPIDActuator)
#endif
