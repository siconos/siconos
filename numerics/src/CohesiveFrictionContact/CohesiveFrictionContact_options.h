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

  /** Anderson acceleration: Tikhonov regularization of the least-squares problem,
   *  relative to the largest diagonal entry of dF^T dF (default 1e-10) */
  SICONOS_COHESIVE_FRICTION_DPARAM_ANDERSON_REGULARIZATION = 20,

  /** Anderson acceleration: the history is cleared when the fixed-point residual
   *  grows by more than this factor between two sweeps (default 1.0) */
  SICONOS_COHESIVE_FRICTION_DPARAM_ANDERSON_SAFEGUARD = 21,

  /** Compliance coefficient epsilon (epsilon * I added to the diagonal blocks of M,
   *  see SICONOS_COHESIVE_FRICTION_IPARAM_COMPLIANCE) */
  SICONOS_COHESIVE_FRICTION_DPARAM_COMPLIANCE = 22,

};

/**
 * Index for integer parameters in SolverOptions->iparam
 */
enum COHESIVE_FRICTION_IPARAM {

  /** Current contact number (for internal use) */
  SICONOS_COHESIVE_FRICTION_IPARAM_CURRENT_CONTACT_NUMBER = 5,

  /** Anderson acceleration of the NSGS sweep: memory depth m (0 = off, default) */
  SICONOS_COHESIVE_FRICTION_IPARAM_ANDERSON_MEMORY = 20,

  /** Compliance strategy (see SICONOS_COHESIVE_FRICTION_COMPLIANCE) */
  SICONOS_COHESIVE_FRICTION_IPARAM_COMPLIANCE = 22,

};

/**
 * Compliance values for SICONOS_COHESIVE_FRICTION_IPARAM_COMPLIANCE.
 *
 * The problem is solved with u = (M + epsilon D) r + q, where D is the identity on the
 * selected blocks and epsilon = dparam[SICONOS_COHESIVE_FRICTION_DPARAM_COMPLIANCE].
 * This removes the hyperstaticity (self-stress modes) of the rigid problem, so that the
 * reactions are unique and bounded by O(incompatibility / epsilon).  The returned
 * velocity and the residual are those of the regularized problem.
 */
enum SICONOS_COHESIVE_FRICTION_COMPLIANCE {
  /** No compliance (rigid problem) */
  SICONOS_COHESIVE_FRICTION_COMPLIANCE_NO = 0,
  /** Add epsilon * I to the diagonal blocks of the cohesive points only */
  SICONOS_COHESIVE_FRICTION_COMPLIANCE_COHESIVE = 1,
  /** Add epsilon * I to the diagonal blocks of the contacts and the cohesive points */
  SICONOS_COHESIVE_FRICTION_COMPLIANCE_ALL = 2,
};

#endif  // COHESIVEFRICTIONCONTACT_OPTIONS_H
