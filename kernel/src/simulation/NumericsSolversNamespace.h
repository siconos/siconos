/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
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
