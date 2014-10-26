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

/*! \file ExplicitLinearSMC.hpp
  \brief General interface to define an actuator
  */

#ifndef ExplicitLinearSMC_H
#define ExplicitLinearSMC_H

#include "CommonSMC.hpp"

class ExplicitLinearSMC : public CommonSMC
{
private:
  /** default constructor */
  ExplicitLinearSMC() {};

  /** serialization hooks */
  ACCEPT_SERIALIZATION(ExplicitLinearSMC);

  /** \f$\sigma = Cx\f$ */
  SP::SiconosVector _sigma;

public:

  /** Constructor.
   * \param sensor the ControlSensor feeding the Actuator
   */
  ExplicitLinearSMC(SP::ControlSensor sensor);

  /** destructor
  */
  virtual ~ExplicitLinearSMC();

  /** Initializer
   * \param m the Model of the Simulation
   */
  virtual void initialize(const Model& m);

  /** Compute the new control law at each event
   * Here we are using the following formula:
   * TODO
   */
  void actuate();

};
#endif
