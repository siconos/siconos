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

/*! \file sampledPIDActuator.hpp
  \brief General interface to define an actuator
*/

#ifndef sampledPIDActuator_H
#define sampledPIDActuator_H

#include "SiconosKernel.hpp"
#include <boost/circular_buffer.hpp>

class sampledPIDActuator : public Actuator
{
private:
  /** default constructor */
  sampledPIDActuator();

  /** error vector */
  boost::shared_ptr<boost::circular_buffer<double> > _err;

  /** dimension of the state space */
  unsigned int _nDim;

  /** reference we are tracking */
  double _ref;

  /** control variable */
  SP::SimpleVector _u;

  /** vector of gains */
  SP::SimpleVector _K;

  /** the sensor that feed the controller */
  SP::controlSensor _sensor;

  /** the dynamical system we are controlling */
  SP::FirstOrderLinearDS _DS;

  /** boolean to determined if the controller has been correctly initialized */
  bool _initDone;

  /** current \f$\Deltat\f$ (or timeStep) */
  double _curDeltaT;

public:

  /** Constructor with a TimeDiscretisation.
   * \param a string, the type of the Actuator, which corresponds to the class type
   * \param a SP::TimeDiscretisation (/!\ it should not be used elsewhere !)
   * \param a SP::Model
   */
  sampledPIDActuator(int, SP::TimeDiscretisation, SP::Model);

  /** Constructor with a TimeDiscretisation.
   * \param a string, the type of the Actuator, which corresponds to the class type
   * \param a SP::TimeDiscretisation (/!\ it should not be used elsewhere !)
   * \param a SP::Model
   * \param a set of Sensor linked to this Actuator.
   */
  sampledPIDActuator(int, SP::TimeDiscretisation, SP::Model, const Sensors&);

  /** destructor
   */
  virtual ~sampledPIDActuator();

  /** initialize actuator data.
   */
  virtual void initialize();

  /** Compute the new control law at each event
   * Here we are using the following formula:
   * \f$u_k = u_{k-1} + c_1e_k + c_2e_{k-1} + c_3e_{k-2}\f$, where
   * \f{align*}c_1 &= K_P - \frac{K_D}/{\Delta t} + K_I \Delta t \\
   * c_2 &= -1 - \frac{2K_D}/{\Delta t} \\
   * c_3 &= \frac{K_D}/{\Delta t} \f{align*}
   */
  void actuate();

  /** Set the value of _K to newValue
   * * \param a SimpleVector \f$[K_P, K_I, K_D]\f$
   */
  void setK(const SimpleVector&);

  /** Set _K to pointer newPtr
   * \param a SP::SimpleVector
   */
  void setKPtr(SP::SimpleVector);

  /** Set the value of _ref to newValue
   * \param a double
   */
  void inline setRef(const double newValue)
  {
    _ref = newValue;
  }

  /** Encapsulates an operation of dynamic casting. Needed by Python interface.
  * \param Actuator*
  * \return a pointer on the derived type
  */
  static sampledPIDActuator* convert(Actuator* s);
};
DEFINE_SPTR(sampledPIDActuator)
#endif
