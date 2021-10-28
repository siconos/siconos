/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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

/*! \file SimulationTypeDef.hpp
 * \brief Typedef for simulation-related objects
 */

#ifndef SimulationTypedef_H
#define SimulationTypedef_H

#include <vector>
#include <set>
#include "SiconosPointers.hpp" // for TYPEDEF_SPTR
#include "SiconosFwd.hpp" // for SP::OneStepIntegrator, ...

// ================== Objects to handle DS ==================

/** list of indices */
typedef std::vector<unsigned int> IndexInt;
TYPEDEF_SPTR(IndexInt)

// ================== Objects to handle OSI ==================

/** Vector of OneStepIntegrator */
typedef std::set<SP::OneStepIntegrator> OSISet;
TYPEDEF_SPTR(OSISet)

/** Iterator through vector of OSI*/
typedef OSISet::iterator OSIIterator;

// ================== Objects to handle OSNS ==================

/** Map of OSNS */
typedef std::vector<SP::OneStepNSProblem> OneStepNSProblems;

/** Iterator through OneStepNSProblems */
typedef OneStepNSProblems::iterator OSNSIterator;
TYPEDEF_SPTR(OneStepNSProblems)

// ================== Misc ==================

/** default tolerance value, used to update index sets */
#define DEFAULT_TOLERANCE 10 * MACHINE_PREC

enum SICONOS_OSNSP
{
  SICONOS_OSNSP_DEFAULT = 0
};
enum SICONOS_OSNSP_ED
{
  SICONOS_OSNSP_ED_SMOOTH_ACC,
  SICONOS_OSNSP_ED_IMPACT,
  SICONOS_OSNSP_ED_SMOOTH_POS
};
enum SICONOS_OSNSP_TS
{
  SICONOS_OSNSP_TS_VELOCITY = 0,
  SICONOS_OSNSP_TS_POS = 1
};
const int SICONOS_NB_OSNSP_TS = 1;
const int SICONOS_NB_OSNSP_TSP = 2;

/** Event constants */
#define TD_EVENT 1
#define NS_EVENT 2


#endif
