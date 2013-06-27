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

/*! \file PID.hpp
  \brief General interface to define an actuator
*/

#ifndef PID_H
#define PID_H

#include "Actuator.hpp"
#include "SiconosAlgebraTypeDef.hpp"
#include <boost/circular_buffer.hpp>

class PID : public Actuator
{
private:
  /** default constructor */
  PID() {};

  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(PID);

  /** error vector */
  std11::shared_ptr<boost::circular_buffer<double> > _err;

  /** reference we are tracking */
  double _ref;

  double _curDeltaT;

  /** vector of gains */
  SP::SiconosVector _K;

public:

  /** Constructor.
   * \param sensor the ControlSensor feeding the Actuator
   */
  PID(SP::ControlSensor sensor);

  /** destructor
   */
  virtual ~PID();

  /** initialize actuator data.
   * \param m a SP::Model
   */
  virtual void initialize(const Model& m);

  /** Compute the new control law at each event
   * Here we are using the following formula:
   * \f$ u_k = u_{k-1} + c_1 e_k + c_2 e_{k-1} + c_3 e_{k-2} \f$ , where
   * \f{array} c_1 &= K_P - \frac{K_D}{\Delta t} + K_I \Delta t \\
   * c_2 &= -1 - \frac{2K_D}{\Delta t} \\
   * c_3 &= \frac{K_D}{\Delta t} \\
   * \f}
   */
  void actuate();

  /** Set the value of _K to newValue
   * * \param newValue SiconosVector \f$ [K_P, K_I, K_D] \f$
   */
  void setK(const SiconosVector& newValue);

  /** Set _K to pointer newPtr
   * \param newPtr SP::SiconosVector \f$ [K_P, K_I, K_D] \f$
   */
  void setKPtr(SP::SiconosVector newPtr);

  /** Set the value of _ref to newValue
   * \param newValue
   */
  void inline setRef(const double newValue)
  {
    _ref = newValue;
  }

  /** Get the timestep from the TimeDiscretisation associated with this PID controller
  *  \param td the TimeDiscretisation for this Actuator
  */
  virtual void setTimeDiscretisation(const TimeDiscretisation& td);

/** display the data of the Actuator on the standard output
   */
  virtual void display() const;

};
DEFINE_SPTR(PID)
#endif
