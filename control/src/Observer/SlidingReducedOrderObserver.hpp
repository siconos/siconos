/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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

/*! \file SlidingReducedOrderObserver.hpp
  \brief Discrete-time Sliding Observer
  */

#ifndef SlidingReducedOrderObserver_H
#define SlidingReducedOrderObserver_H

#include "Observer.hpp"
#include "SiconosAlgebraTypeDef.hpp"

class SlidingReducedOrderObserver : public Observer
{
private:
  /** serialization hooks */
  ACCEPT_SERIALIZATION(SlidingReducedOrderObserver);

  /** default constructor */
  SlidingReducedOrderObserver() {};

protected:

  /** the vector defining the measurements (\f$ y = Cx \f$) */
  SP::SimpleMatrix _C;

  /** matrix multiplying the innovation term */
  SP::SimpleMatrix _L;

  double _theta;

  // clumsy hack to do nothing the first time this Observer is called
  bool _pass;

public:

  /** Constructor with the standard interface
   * \param sensor the SP::ControlSensor that feed us with measurements
   * \param xHat0 the initial guess for the state
   */
  SlidingReducedOrderObserver(SP::ControlSensor sensor, const SiconosVector& xHat0):
    Observer(SLIDING_REDUCED_ORDER, sensor, xHat0), _pass(false) {}

  /** Constructor with all the data
   * \param sensor the sensor that feeds the Observer
   * \param xHat0 the initial guess for the state
   * \param C observation matrix
   * \param L gain matrix
   */
  SlidingReducedOrderObserver(SP::ControlSensor sensor, const SiconosVector& xHat0, SP::SimpleMatrix C, SP::SimpleMatrix L):
    Observer(SLIDING_REDUCED_ORDER, sensor, xHat0), _C(C), _L(L), _pass(false) {}

  /** Update the estimate at each event
  */
  virtual void process();

  /** Initialization
   * \param nsds current nonsmooth dynamical system
   * \param s current simulation setup
   */
  virtual void initialize(const NonSmoothDynamicalSystem& nsds, const Simulation& s);

  /** Set the L matrix
   * \param L the new L matrix
   */
  inline void setLPtr(SP::SimpleMatrix L)
  {
    _L = L;
  };

  /** Set the C matrix
   * \param C the new C matrix
   */
  inline void setCPtr(SP::SimpleMatrix C)
  {
    _C = C;
  };

};
#endif
