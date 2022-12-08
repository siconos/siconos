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

/*!\file
 * TimeSteppingD1Minus simulation
 */

#ifndef TIMESTEPPINGD1MINUS_H
#define TIMESTEPPINGD1MINUS_H

#include "Simulation.hpp"

/** TimeSteppingD1Minus Timestepping Strategy
 *
 *  see Schindler/Acary : Timestepping Schemes for Nonsmooth Dynamics Based
 *  on Discontinuous Galerkin Methods: Definition and Outlook
 */
class TimeSteppingD1Minus : public Simulation {
private:

  ACCEPT_SERIALIZATION(TimeSteppingD1Minus);

  /** default constructor */
  TimeSteppingD1Minus() {}

protected:
  /** initialisation specific to TimeSteppingD1Minus for OneStepNSProblem */
  void initOSNS() override;

public:
  /** constructor with the time-discretisation
   * \param nsds the current nonsmooth dynamical system
   * \param td pointer to a TimeDiscretisation
   * \param nb number of non smooth problem
   */
  TimeSteppingD1Minus(SP::NonSmoothDynamicalSystem nsds,
                      SP::TimeDiscretisation td, int nb);

  /** destructor */
  ~TimeSteppingD1Minus();

  /** updateIndexSet using current y and lambda values of interactions
   *  \param i the  number of the set to be updated
   *  0 : ALL interactions (NEVER)
   *  1 : ACTIVE interactions for IMPACTS
   *  2 : ACTIVE interactions for CONTACTS
   */
  void updateIndexSet(unsigned int i) override;

  /** run the simulation, from t0 to T */
  void run() override;

  /** step from current event to next event of EventsManager */
  void advanceToEvent() override;

  /** compute residu */
  void computeResidu();

  /** integrate DynamicalSystems taking not into account non-smooth part */
  void computeFreeState();

  ACCEPT_STD_VISITORS();
};

DEFINE_SPTR(TimeSteppingD1Minus)

#endif // TIMESTEPPINGD1MINUS_H
