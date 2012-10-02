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

/*! \file LinearChatteringSMC.hpp
  \brief General interface to define an actuator
  */

#ifndef LinearChatteringSMC_H
#define LinearChatteringSMC_H

#include "SiconosKernel.hpp"
#include <boost/circular_buffer.hpp>

class LinearChatteringSMC : public CommonSMC
{
private:
  /** default constructor */
  LinearChatteringSMC() {};

  /** serialization hooks */
  ACCEPT_SERIALIZATION(LinearChatteringSMC);

  /** sign of Cx*/
  SP::SiconosVector _s;

public:

  /** Constructor with a TimeDiscretisation and a Model.
   * \param t a SP::TimeDiscretisation (/!\ it should not be used elsewhere !)
   * \param ds the SP::DynamicalSystem we are controlling
   */
  LinearChatteringSMC(SP::TimeDiscretisation t, SP::DynamicalSystem ds);

  /** Constructor with a TimeDiscretisation, a Model and a set of Sensor.
   * \param t a SP::TimeDiscretisation (/!\ it should not be used elsewhere !)
   * \param ds the SP::DynamicalSystem we are controlling
   * \param sensorList a set of Sensor linked to this Actuator.
   */
  LinearChatteringSMC(SP::TimeDiscretisation t, SP::DynamicalSystem ds, const Sensors& sensorList);

  /** destructor
  */
  virtual ~LinearChatteringSMC();

  /** initialize actuator data.
   * \param m a SP::Model
  */
  void initialize(SP::Model m);

  /** Compute the new control law at each event
   * Here we are using the following formula:
   * TODO
   */
  void actuate();

};
DEFINE_SPTR(LinearChatteringSMC)
#endif
