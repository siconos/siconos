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

/*! \file LinearSMCimproved.hpp
  \brief General interface to define an actuator
*/

#ifndef LinearSMCimproved_H
#define LinearSMCimproved_H

#include "LinearSMC.hpp"

#include <boost/circular_buffer_fwd.hpp>

typedef std11::shared_ptr<boost::circular_buffer<SP::SiconosVector> > BufferOfVectors;

class LinearSMCimproved : public LinearSMC
{
private:
  /** serialization hooks */
  ACCEPT_SERIALIZATION(LinearSMCimproved);

  /** default constructor */
  LinearSMCimproved() {};

protected:

  /** try to predict the perturbation */
  bool _predictionPerturbation;

  /** boolean to determine if we are in the discrete-time sliding phase */
  bool _inDisceteTimeSlidingPhase;

  /** Vector to store previous values of the perturbation */
  BufferOfVectors _measuredPert;

  /** Vector of predicted values for the perturbation */
  BufferOfVectors _predictedPert;

  /** Upper bound on the norm2 of the perturbation */
  double _ubPerturbation;

  /** Control input to counteract the effect of the perturbation */
  SP::SiconosVector _up;

  /** Predict the effect of the perturnation during the next timestep
   * \param xTk available state at the current time instant
   * \param CBstart matrix \f$CB^{*}\f$
   */
  void predictionPerturbation(const SiconosVector& xTk, SimpleMatrix& CBstar);


public:

  /** Constructor
   * \param sensor the ControlSensor feeding the Actuator
   */
  LinearSMCimproved(SP::ControlSensor sensor);

  /** Constructor with all the data
   * \param sensor the ControlSensor feeding the Actuator
   * \param B the B matrix in the FirstOrderLinearR
   * \param D the D matrix in the FirstOrderLinearR
   */
  LinearSMCimproved(SP::ControlSensor sensor, SP::SimpleMatrix B, SP::SimpleMatrix D);

  /** destructor
   */
  virtual ~LinearSMCimproved();

  /** Initialize Controller
   * \param m the Model
   */
  virtual void initialize(const Model& m);

  /** Compute the new control law at each event
   * Here we are using the following formula:
   * TODO
   */
  virtual void actuate();

  /** Enable perturbation prediction */
  void setPerturbationPrediction(double ub = std::numeric_limits<double>::quiet_NaN())
  {
    _ubPerturbation = ub;
    _predictionPerturbation = true;
  }

  /** Get the control input _up, acting against matched perturbations
   * \return a reference to _up
   */
  inline const SiconosVector& up() const { return *_up; };

  /** Set the order of the prediction
   * - 0 -> the predicted value is the same as the one we measured
   * - 1 -> \f$\widetilde{Cp_{k+1}} = 2Cp_k - Cp_{k-1}\f$
   * - 2 -> \f$\widetilde{Cp_{k+1}} = 3Cp_k - 3Cp_{k-1} + Cp_{k-2}\f$
   * \param order the order of the prediction
   */
  void setPredictionOrder(unsigned int order);
};
#endif
