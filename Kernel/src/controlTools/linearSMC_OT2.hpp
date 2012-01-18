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

/*! \file linearSMC_OT2.hpp
  \brief General interface to define a sliding mode controller with
  disturbance compensation. Reference: Su, W.C.; Drakunov, S.; Özgüner, Ü.
  An O(T2) boundary layer in sliding mode for sampled-data systems
  */

#ifndef linearSMC_OT2_H
#define linearSMC_OT2_H

#include "SiconosKernel.hpp"

//XXX: to be removed
#include "commonSMC.hpp"
class linearSMC_OT2 : public commonSMC
{
private:
  /** default constructor */
  linearSMC_OT2();
  /** Current value of the state (\f$ x_k\f$)*/
  SP::SimpleVector _X;
  /** Predicted current value of the state (\f$ \hat{x}_k = \Phi x_{k-1} + \Gamma u_{k-1}\f$)*/
  SP::SimpleVector _Xhat;
  /** Next value of the state only with the influence of the dynamic \f$ \XPhi = \Phi x_k\f$*/
  SP::SimpleVector _XPhi;
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
   * \param a string, the type of the Actuator, which corresponds to the class type
   * \param a SP::TimeDiscretisation (/!\ it should not be used elsewhere !)
   * \param a SP::Model
   */
  linearSMC_OT2(int, SP::TimeDiscretisation, SP::Model);

  /** Constructor with a TimeDiscretisation, a Model and a set of Sensor.
   * \param a string, the type of the Actuator, which corresponds to the class type
   * \param a SP::TimeDiscretisation (/!\ it should not be used elsewhere !)
   * \param a SP::Model
   * \param a set of Sensor linked to this Actuator.
   */
  linearSMC_OT2(int, SP::TimeDiscretisation, SP::Model, const Sensors&);

  /** destructor
  */
  virtual ~linearSMC_OT2();

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

  /** Set the value of _Csurface to newValue
   * * \param a SimpleVector newValue
   */
  void setCsurface(const SimpleVector&);

  /** Set _Csurface to pointer newPtr
   * \param a SP::SimpleVector
   */
  void setCsurfacePtr(SP::SimpleVector);
};
DEFINE_SPTR(linearSMC_OT2)
#endif
