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
#ifndef FRICTIONCONTACT3D_onecontact_nonsmooth_Newton_solvers_H
#define FRICTIONCONTACT3D_onecontact_nonsmooth_Newton_solvers_H

/*!\file fc3d_onecontact_nonsmooth_Newton_solvers.h
  \brief Typedef and functions declarations related to Newton solver for 3 dimension frictional contact problems.

  Each solver must have 4 functions in its interface:
  - initialize: link local static variables to the global ones (M,q,...)
  - update: link/fill the local variables corresponding to sub-blocks of the full problem, for a specific contact
  - solve: solve the local problem
  - free

*/
#include "NumericsMatrix.h"
#include "SolverOptions.h"
#include "FrictionContactProblem.h"
#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

enum SICONOS_FRICTION_3D_ONECONTACT_NSN_FORMULATION_ENUM
{
  SICONOS_FRICTION_3D_ONECONTACT_NSN_FORMULATION_ALARTCURNIER_STD =0,
  SICONOS_FRICTION_3D_ONECONTACT_NSN_FORMULATION_JEANMOREAU_STD =1,
  SICONOS_FRICTION_3D_ONECONTACT_NSN_FORMULATION_ALARTCURNIER_GENERATED =2,
  SICONOS_FRICTION_3D_ONECONTACT_NSN_FORMULATION_JEANMOREAU_GENERATED =3,
  SICONOS_FRICTION_3D_ONECONTACT_NSN_FORMULATION_NULL = 4 ,
  /** index in iparam to store the formulation */
  SICONOS_FRICTION_3D_ONECONTACT_NSN_FORMULATION = 10,
};
enum SICONOS_FRICTION_3D_ONECONTACT_NSN_LINESEARCH_ENUM
{
  SICONOS_FRICTION_3D_ONECONTACT_NSN_LINESEARCH_GOLDSTEINPRICE = 0 ,
  SICONOS_FRICTION_3D_ONECONTACT_NSN_LINESEARCH_ARMIJO = 1,
  /** index in iparam to store the line-search */
  SICONOS_FRICTION_3D_ONECONTACT_NSN_LINESEARCH = 11,
  /** index in iparam to store the maximum number of iterations */
  SICONOS_FRICTION_3D_ONECONTACT_NSN_LINESEARCH_MAXITER = 12
};

enum SICONOS_FRICTION_3D_ONECONTACT_NSN_HYBRID
{
  /** Maximum number of loop for the hybrid solver */
  SICONOS_FRICTION_3D_ONECONTACT_NSN_HYBRID_MAX_LOOP = 14,
  /** Maximum number of iterations for the projection solver */
  SICONOS_FRICTION_3D_ONECONTACT_NSN_HYBRID_MAX_ITER = 15
};


typedef void (*computeNonsmoothFunction)(double *, double * , double , double * , double *, double *, double *);

  /** initialize friction-contact 3D Newton solver
   * \param problem to solve
   * \param localproblem to solve
   * \param options of the solver
   */
  void fc3d_onecontact_nonsmooth_Newton_solvers_initialize(FrictionContactProblem* problem, FrictionContactProblem* localproblem, SolverOptions * options);

  /** solve friction-contact 3D problem with Newton
   * \param localproblem to solve
   * \param options of the solver
   * \return 0 iff successful.
   */
  int fc3d_onecontact_nonsmooth_Newton_solvers_solve(FrictionContactProblem* localproblem, double*, SolverOptions * options);

  /** free memory for friction contact 3D Newton solver
      \param problem the global problem to solve
      \param localproblem for freeing matrix0
      \param localsolver_options options of the solver
   */
  void fc3d_onecontact_nonsmooth_Newton_solvers_free(FrictionContactProblem * problem, FrictionContactProblem * localproblem, SolverOptions* localsolver_options);

  /** compute error for friction-contact 3D problem with Newton
   *  \param dimension of the global problem
   *   \param[in,out] velocity vector (\warning in-out parameter )
   *   \param reaction global reaction vector
   *   \param output_error
   */
  void fc3d_onecontact_nonsmooth_Newton_solvers_computeError(int dimension, double* velocity, double*reaction, double * output_error);

  /** Update friction-contact 3D problem: formalize local problem for one contact
      \param problem the global problem to solve
      \param localproblem the local problem to solve
      \param number (position in global matrix) of the considered contact
      \param reaction global reaction (only the block corresponding to the
      current contact will be modified
      \param options of the solver

      the rest is used to formalize the local problem)
  */
  void fc3d_onecontact_nonsmooth_Newton_AC_update(int number, FrictionContactProblem* problem, FrictionContactProblem* localproblem ,
                                   double * reaction, SolverOptions* options);

  int fc3d_onecontact_nonsmooth_Newton_solvers_solve_direct(FrictionContactProblem* localproblem,
                                 double * R, int *iparam, double *dparam);

  int fc3d_onecontact_nonsmooth_Newton_solvers_solve_damped(FrictionContactProblem* localproblem,
                                       double * R, int *iparam, double *dparam);


  /* Set the default solver options for the ONECONTACT_NSN_AC Solver
   * Some default values:
   * options.iparam[0] = 200 is the maximum number of iterations.
   * options.iparam[3] = 100000 is the nzmax parameter for sparse matrices.
   * options.iparam[10] = 0 : stands for STD Alart & Curnier function
   *  (other values may be 1 for JeanMoreau, 2 for STD generated, 3 for JeanMoreau generated)
   * options.dparam[0] = 1e-3 precision.
   * \param options  the solver options
   */
  int fc3d_onecontact_nonsmooth_Newtow_setDefaultSolverOptions(SolverOptions* options);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
