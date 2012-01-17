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

/*! \file commonSMC.hpp
  \brief General interface to define an actuator
  */

#ifndef commonSMC_H
#define commonSMC_H

#include "SiconosKernel.hpp"
#include <boost/circular_buffer.hpp>

class commonSMC : public Actuator
{
private:
  /** default constructor */
  // commonSMC();
  ACCEPT_SERIALIZATION(commonSMC);
protected:
  /** dimension of the output */
  unsigned int _nDim;

  /** index for saving data */
  unsigned int _indx;

  /** control variable */
  SP::SimpleVector _u;

  /** the vector defining the surface (\f$ s = Cx\f$) */
  SP::SimpleVector _Csurface;

  /** the sensor that feed the controller */
  SP::controlSensor _sensor;

  /** the dynamical system we are controlling */
  SP::FirstOrderLinearDS _DS;

  /** boolean to determined if the controller has been correctly initialized */
  bool _initDone;

  /** current \f$\Deltat\f$ (or timeStep) */
  double _curDeltaT;

  /** matrix describing the relation between the control value and sgn(s) */
  SP::SimpleMatrix _B;

  /** matrix describing the influence of \f$lambda\f$ on s */
  SP::SimpleMatrix _D;

  /** the Relation for the Controller */
  SP::FirstOrderLinearR _relationSMC;

  /** the NonSmoothLaw for the controller */
  SP::NonSmoothLaw _sign;

  /** */
  SP::Interaction _interactionSMC;

  /** easy access to lambda */
  SP::SiconosVector _lambda;

  /** easy access to the state */
  SP::SiconosVector _xController;

public:

  /** Constructor with a TimeDiscretisation and a Model.
   * \param an integer, the type of the Actuator, which corresponds to the class type
   * \param a SP::TimeDiscretisation (/!\ it should not be used elsewhere !)
   * \param a SP::Model
   */
  commonSMC(int i, SP::TimeDiscretisation t, SP::Model m): Actuator(i, t, m) {}

  /** Constructor with a TimeDiscretisation, a Model and a set of Sensor.
   * \param a string, the type of the Actuator, which corresponds to the class type
   * \param a SP::TimeDiscretisation (/!\ it should not be used elsewhere !)
   * \param a SP::Model
   * \param a set of Sensor linked to this Actuator.
   */
  commonSMC(int i, SP::TimeDiscretisation t, SP::Model m, const Sensors& s): Actuator(i, t, m, s) {}

  /** Compute the new control law at each event
   */
  virtual void actuate() = 0;

  /** Encapsulates an operation of dynamic casting. Needed by Python interface.
   * \param Actuator*
   * \return a pointer on the derived type
   */
  static commonSMC* convert(Actuator* s);
};
DEFINE_SPTR(commonSMC)
#endif
