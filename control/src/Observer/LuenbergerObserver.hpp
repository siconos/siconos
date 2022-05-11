/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

/*! \file LuenbergerObserver.hpp
  \brief Classical discrete-time Luenberger Observer
  */

#ifndef LuenbergerObserver_H
#define LuenbergerObserver_H

#include "Observer.hpp"
#include "SiconosAlgebraTypeDef.hpp"

class LuenbergerObserver : public Observer
{
private:
  /** serialization hooks */
  ACCEPT_SERIALIZATION(LuenbergerObserver);

  /** default constructor */
  LuenbergerObserver() {};

protected:

  /** the vector defining the measurements (\f$ y = Cx \f$) */
  SP::SiconosMatrix _C;

  /** matrix describing the relation between the control value and sgn(s) */
  SP::SiconosMatrix _L;

  double _theta;

  // clumsy hack to do nothing the first time this Observer is called
  bool _pass;

public:

  /** Constructor with a TimeDiscretisation, a ControlSensor and an initial estimate of the state.
   * \param sensor the SP::ControlSensor that feed us with measurements
   * \param xHat0 the initial guess for the state
   */
  LuenbergerObserver(SP::ControlSensor sensor, const SiconosVector& xHat0):
    Observer(LUENBERGER, sensor, xHat0), _pass(false) {}

  /** Constructor with all the data
   * \param sensor the ControlSensor that feeds the Observer
   * \param xHat0 the initial guess for the state
   * \param C the observation matrix
   * \param L the gain matrix
   */
  LuenbergerObserver(SP::ControlSensor sensor, const SiconosVector& xHat0, SP::SiconosMatrix C, SP::SiconosMatrix L):
    Observer(LUENBERGER, sensor, xHat0), _C(C), _L(L), _pass(false) {}

  /** Compute the new control law at each event
  */
  virtual void process();

  /** Initialization
   * \param nsds current nonsmooth dynamical system
   * \param s current simulation setup
   */
  virtual void initialize(const NonSmoothDynamicalSystem & nsds, const Simulation& s);

  /** Set the L matrix
   * \param L the new L matrix
   */
  inline void setLPtr(SP::SiconosMatrix L)
  {
    _L = L;
  };

  /** Set the C matrix
   * \param C the new C matrix
   */
  inline void setCPtr(SP::SiconosMatrix C)
  {
    _C = C;
  };

};
#endif
