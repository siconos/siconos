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

/*! \file LinearSMCOT2.hpp
  \brief General interface to define a sliding mode controller with
  disturbance compensation. Reference: Su, W.C.; Drakunov, S.; Özgüner, Ü.
  An O(T2) boundary layer in sliding mode for sampled-data systems
  */

#ifndef LinearSMCOT2_H
#define LinearSMCOT2_H

#include "SiconosKernel.hpp"

class LinearSMCOT2 : public CommonSMC
{
private:
  /** default constructor */
  LinearSMCOT2() {};
  /** Current value of the state (\f$ x_k\f$)*/
  SP::SimpleVector _X;
  /** Predicted current value of the state (\f$ \hat{x}_k = \Phi x_{k-1} + \Gamma u_{k-1}\f$)*/
  SP::SiconosVector _Xhat;
  /** Next value of the state only with the influence of the dynamic \f$ \XPhi = \Phi x_k\f$*/
  SP::SiconosVector _XPhi;
  /** Model for the computation of _XPhi*/
  SP::Model _modelPhi;
  /** DynamicalSystem for the computation of _XPhi*/
  SP::FirstOrderLinearDS _DSPhi;
  /** TimeDiscretisation for the computation of _XPhi*/
  SP::TimeDiscretisation _timeDPhi;
  /** OneSteoIntegrator for the computation of _XPhi*/
  SP::Lsodar _PhiOSI;
  /** Simulation for the computation of _XPhi*/
  SP::EventDriven _simulPhi;
  /** Model for the computation of Xhat*/
  SP::Model _modelPred;
  /** TimeDiscretisation for the computation of Xhat*/
  SP::TimeDiscretisation _timeDPred;
  /** OneSteoIntegrator for the computation of Xhat*/
  SP::Lsodar _PredOSI;
  /** Simulation for the computation of Xhat*/
  SP::EventDriven _simulPred;
  /** DynamicalSystem for the computation of _Xhat*/
  SP::FirstOrderLinearDS _DSPred;
  /** Coefficient*/
  double _coeff;

public:

  /** Constructor with a TimeDiscretisation and a Model.
   * \param t the SP::TimeDiscretisation (/!\ it should not be used elsewhere !)
   * \param ds the SP::DynamicalSystem we are controlling
   */
  LinearSMCOT2(SP::TimeDiscretisation t, SP::DynamicalSystem ds);

  /** Constructor with a TimeDiscretisation, a Model and a set of Sensor.
   * \param t SP::TimeDiscretisation (/!\ it should not be used elsewhere !)
   * \param ds the SP::DynamicalSystem we are controlling
   * \param sensorList a set of Sensor linked to this Actuator.
   */
  LinearSMCOT2(SP::TimeDiscretisation t, SP::DynamicalSystem ds, const Sensors& sensorList);

  /** destructor
  */
  virtual ~LinearSMCOT2();

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
DEFINE_SPTR(LinearSMCOT2)
#endif
