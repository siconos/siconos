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

/*! \file LinearSMC.hpp
  \brief General interface to define an actuator
*/

#ifndef LinearSMC_H
#define LinearSMC_H

#include "CommonSMC.hpp"


class LinearSMC : public CommonSMC
{
private:
  /** serialization hooks */
  ACCEPT_SERIALIZATION(LinearSMC);


protected:
  /** default constructor */
  LinearSMC() {};


public:

  /** Constructor with a TimeDiscretisation and a DynamicalSystem.
   * \param t a SP::TimeDiscretisation (/!\ it should not be used elsewhere !)
   * \param ds the SP::DynamicalSystem we are controlling
   * \param name do not set this yourself ! this is used in derived classes
   */
  LinearSMC(SP::TimeDiscretisation t, SP::DynamicalSystem ds, int name = LINEAR_SMC);

  /** Constructor with a TimeDiscretisation and a DynamicalSystem.
   * \param t a SP::TimeDiscretisation (/!\ it should not be used elsewhere !)
   * \param ds the SP::DynamicalSystem we are controlling
   * \param B the B matrix in the FirstOrderLinearR
   * \param D the D matrix in the FirstOrderLinearR
   * \param name do not set this yourself ! this is used in derived classes
   */
  LinearSMC(SP::TimeDiscretisation t, SP::DynamicalSystem ds, SP::SiconosMatrix B, SP::SiconosMatrix D, int name = LINEAR_SMC);

  /** Constructor with a TimeDiscretisation, a DynamicalSystem and a set of Sensor.
   * \param t a SP::TimeDiscretisation (/!\ it should not be used elsewhere !)
   * \param ds the SP::DynamicalSystem we are controlling
   * \param sensorList a set of Sensor linked to this Actuator.
   * \param name do not set this yourself ! this is used in derived classes
   */
  LinearSMC(SP::TimeDiscretisation t, SP::DynamicalSystem ds, const Sensors& sensorList, int name = LINEAR_SMC);

  /** destructor
   */
  virtual ~LinearSMC();

  /** Compute the new control law at each event
   * Here we are using the following formula:
   * TODO
   */
  virtual void actuate();


  /** Set the D matrix
   * \param D the new D matrix
  */
  void setD(const SiconosMatrix & D);

    /** Set the D matrix
   * \param D the new D matrix
  */
  inline void setDPtr(SP::SiconosMatrix D)
  {
    _D = D;
  };


};
DEFINE_SPTR(LinearSMC)
#endif
