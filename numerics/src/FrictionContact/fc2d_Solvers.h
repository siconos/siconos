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
#ifndef fc2dSolvers_H
#define fc2dSolvers_H

/*!\file fc2d_Solvers.h
  \brief Subroutines for the resolution of contact problems with friction (2-dimensional case).
*/

#include "FrictionContactProblem.h"
#include "SolverOptions.h"
#include "Friction_cst.h"
#include "LinearComplementarityProblem.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /**  cpg (conjugated projected gradient) solver for global contact problems with friction (2D)
       \param[in]  problem the friction-contact problem
       \param[out] reaction vector
       \param[out] velocity vector
       \param[in,out] info termination value
       \param[in,out] options structure for options
  */
  void fc2d_cpg(FrictionContactProblem* problem , double *reaction , double *velocity , int *info, SolverOptions* options);

  /**  Non Linear Gauss Seidel solver (dense) for global contact problem with friction in 2D case.
       \param[in] problem the friction-contact problem
       \param[out] reaction vector
       \param[out] velocity vector
       \param[in,out] info termination value
       \param[in,out] options structure
  */
  void fc2d_nsgs_dense(FrictionContactProblem* problem , double *reaction , double *velocity , int *info, SolverOptions* options);

  /**  Non Linear Gauss Seidel solver (sbm) for global contact problem with friction in 2D case.
       \param[in] problem the friction-contact problem
       \param[out] reaction vector
       \param[out] velocity vector
       \param[in,out] info termination value
       \param[in,out] options structure
  */
  void fc2d_nsgs_sbm(FrictionContactProblem* problem, double *z, double *w, int *info, SolverOptions* options);

  
  /** fc2d_projc is a specific projection operator related to CPG (conjugated projected gradient) algorithm for global contact problem with friction.
   *
   *
   * \param[in] xi  the intermediate iterate which goes to be projected (projc1).
   * \param[in] n   the dimension of the system.
   * \param[in] statusi  a vector which contains the initial status.
   * \param[in] p       a vector which contains the components of the descent direction.
   * \param[in] fric a vector which contains the friction coefficient.
   * \param[out] reaction the corrected iterate.
   * \param[out] status  the new status.
   *
   */
  void fc2d_projc(double xi[], int *n, int statusi[], double p[], double fric[], double *reaction, int *status);

  /** fc2d_projf is a specific projection operator related to CPG (conjugated projected gradient) algorithm
   *              for global contact problem with friction.
   *
   *
   * \param[in] etat  parameter which represents the status vector.
   * \param[in] n      parameter which represents the dimension of the system.
   * \param[in] y    parameter which contains the components of the residue or descent direction vector.
   * \param[in] fric   parameter which contains the friction coefficient.
   * \param[out] projf1 parameter which contains the projected residue or descent direction.
   *
   */
  void fc2d_projf(int etat[], int *n, double y[], double fric[], double projf1[]);



  /** fc2d_lexicolemke is a Lemke solver for  frictionContact2D problems.
     * \param[in] problem structure that represents the fc2d (M, q...)
     * \param[in,out] reaction a n-vector of doubles which contains the initial solution and returns the solution of the problem.
     * \param[in,out] velocity a n-vector of doubles which returns the solution of the problem.
     * \param options
     * \param[out] info an integer which returns the termination value:
     0 = convergence,
     1 = no convergence,
     2 = Null diagonal term
    */
  void fc2d_lexicolemke(FrictionContactProblem* problem, double *reaction, double *velocity, int *info, SolverOptions* options);



  /** This function transform a FrictionContactProblem (2D) into a LinearComplementarityProblem
   * \param[in] problem A pointer to a FrictionContactProblem to transform
   * \param[out] lcp_problem A pointer to a LinearComplementarity_problem resulting from the reformulation
   */

  int fc2d_tolcp(FrictionContactProblem* problem, LinearComplementarityProblem * lcp_problem);

  /** fc2d_enum solver for  frictionContact2D problems.
     * \param[in] problem structure that represents the fc2d (M, q...)
     * \param[in,out] reaction a n-vector of doubles which contains the initial solution and returns the solution of the problem.
     * \param[in,out] velocity a n-vector of doubles which returns the solution of the problem.
     * \param options
     * \param[out] info an integer which returns the termination value:
     0 = convergence,
     1 = no convergence,
     2 = Null diagonal term
    */
  void fc2d_enum(FrictionContactProblem* problem, double *reaction, double *velocity, int *info, SolverOptions* options);

  /** @addtogroup SetSolverOptions
      @{
  */

  void fc2d_nsgs_set_default(SolverOptions* options);

  /** @} */



#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
