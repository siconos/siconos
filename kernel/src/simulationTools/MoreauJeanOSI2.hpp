/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
/*! \file
  MoreauJeanOSI Time-Integrator for Dynamical Systems
*/

#ifndef MOREAU2_H
#define MOREAU2_H

#include "MoreauJeanOSI.hpp"
#include "SimpleMatrix.hpp"
#include "SiconosVector.hpp"
#include "FirstOrderLinearDS.hpp"


/**  MoreauJeanOSI Time-Integrator for Dynamical Systems
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Apr 26, 2004
 *
 * See User's guide for details.
 *
 * MoreauJeanOSI class is used to define some time-integrators methods for a list of dynamical systems.
 * Each DynamicalSystem is associated to a SiconosMatrix, named "W", and a double, "theta", through two
 * STL maps:
 * - thetaMap, thetaMap[ds] = a double
 * ds being a SP::DynamicalSystem
 *
 * W matrices are initialized and computed in initW and computeW. Depending on the DS type, they
 * may depend on time and DS state (x).
 *
 * Main functions:
 *
 * - computeFreeState(): computes Ffree (or vfree), dynamical systems state without taking non-smooth part into account \n
 * - updateState():
 *
 */
class MoreauJeanOSI2 : public MoreauJeanOSI
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(MoreauJeanOSI2);

public:

  /** constructor from a minimum set of data:  theta
   *  \param theta value for the theta parameter
   */
  MoreauJeanOSI2(double theta = 0.5);

  ~MoreauJeanOSI2();

  //  SP::SiconosVector  getFfree(FirstOrderLinearDS *d);

  /** integrates the Dynamical System linked to this integrator without boring the constraints
   */
  void computeFreeState();


  /** updates the state of the Dynamical Systems
   *  \param level level of interest for the dynamics: not used at the time
   */
  void updateState(const unsigned int level);

};

#endif // MOREAU2_H
