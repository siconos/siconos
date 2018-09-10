/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
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
   * \param CBstar matrix \f$CB^{*}\f$
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
  LinearSMCimproved(SP::ControlSensor sensor, SP::SimpleMatrix B, SP::SimpleMatrix D = std11::shared_ptr<SimpleMatrix>());

  /** destructor
   */
  virtual ~LinearSMCimproved();

  /** Initialize Controller
   * \param nsds current nonsmooth dynamical system
   * \param s current simulation setup
   */
  virtual void initialize(const NonSmoothDynamicalSystem& nsds, const Simulation & s);

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
