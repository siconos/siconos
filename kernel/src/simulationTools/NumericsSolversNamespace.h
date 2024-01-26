/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

/*! \file NumericsSolversNamespace.h"
  \brief numerics component solvers interface in a dedicated namespace
 */

#ifndef NumericsSolversNS
#define NumericsSolversNS

// WARNING : do not include this file into another header.
// Cpp include only !

// namespace siconos::numerics {

#include "AffineVariationalInequalities.h"
#include "FrictionContactProblem.h"
#include "GenericMechanicalProblem.h"
#include "GenericMechanical_Solvers.h"
#include "GlobalFrictionContactProblem.h"
#include "GlobalRollingFrictionContactProblem.h"
#include "LCP_Solvers.h"
#include "LinearComplementarityProblem.h"
#include "MLCP_Solvers.h"
#include "MixedLinearComplementarityProblem.h"
#include "NonSmoothDrivers.h"
#include "RelayProblem.h"
#include "RollingFrictionContactProblem.h"
#include "SolverOptions.h"
//}  // namespace siconos::numerics

#endif
