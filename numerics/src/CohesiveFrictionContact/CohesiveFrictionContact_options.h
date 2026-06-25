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
#ifndef COHESIVEFRICTIONCONTACT_OPTIONS_H
#define COHESIVEFRICTIONCONTACT_OPTIONS_H

/*!\file CohesiveFrictionContact_options.h
  Options and solver IDs for cohesive friction-contact problems.

  This file defines the solver identifiers and option parameters for
  solving cohesive friction-contact problems using various algorithms.
*/

/**
 * Solver IDs for 3D Cohesive Friction Contact problems
 */
enum COHESIVE_FRICTION_SOLVER {
  /** Non-Smooth Gauss-Seidel (NSGS) solver with projection on cone */
  SICONOS_COHESIVE_FRICTION_3D_NSGS = 8000,

  /** Projection solver (on projection) */
  SICONOS_COHESIVE_FRICTION_3D_PROJECTION,

  /** NSGS with projection and local iteration */
  SICONOS_COHESIVE_FRICTION_3D_NSGS_PROJECTION_WITH_LOCAL_ITERATION,
  
};

/**
 * Index for double parameters in SolverOptions->dparam
 */
enum COHESIVE_FRICTION_DPARAM {
  
  /** Local tolerance for internal solvers */
  SICONOS_COHESIVE_FRICTION_DPARAM_LOCAL_TOLERANCE,
  
};

/**
 * Index for integer parameters in SolverOptions->iparam
 */
enum COHESIVE_FRICTION_IPARAM {

  /** Current contact number (for internal use) */
  SICONOS_COHESIVE_FRICTION_IPARAM_CURRENT_CONTACT_NUMBER,
  
};

#endif  // COHESIVEFRICTIONCONTACT_OPTIONS_H
