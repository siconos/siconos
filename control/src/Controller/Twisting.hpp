/* Siconos-Kernel, Copyright INRIA 2005-2015
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

/*! \file Twisting.hpp
  \brief twisting algorithm
*/

#ifndef Twisting_H
#define Twisting_H

#include "CommonSMC.hpp"


class Twisting : public CommonSMC
{
private:
  /** serialization hooks */
  ACCEPT_SERIALIZATION(Twisting);


protected:
  /** default constructor */
  Twisting() {};

public:
  /** Constructor for the ActuatorFactory
   * \param sensor the ControlSensor feeding the Actuator
   */
  Twisting(SP::ControlSensor sensor): CommonSMC(TWISTING, sensor) {};

  /** Constructor for a nonlinear system.
   * \param sensor the ControlSensor feeding the Actuator
   * \param hControl sampling period
   */
  Twisting(SP::ControlSensor sensor, double hControl);

  /** Constructor for the linear case
   * \param sensor the ControlSensor feeding the Actuator
   * \param gain control magnitude
   * \param beta twisting parameter
   * \param hControl sampling period
   */
  Twisting(SP::ControlSensor sensor, double gain, double beta, double hControl);

  /** destructor
   */
  virtual ~Twisting();

  /** Compute the new control law at each event
   * Here we are using the following formula:
   */
  virtual void actuate();

  /** set nonsmooth data: NormalConeNSL and AVI osnsp
   * \param hControl sampling period
   */
  void setNSdata(double hControl);

  virtual void initialize(const Model& m);
};
#endif
