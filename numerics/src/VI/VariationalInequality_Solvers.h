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
#ifndef VISOLVERS_H
#define VISOLVERS_H

/*!\file VariationalInequality_Solvers.h
  \brief Subroutines for the resolution of Variational Inequalites (VI) problems
*/

#include "VariationalInequality.h"
#include "SolverOptions.h"
#include "VI_cst.h"
#include "SiconosCompat.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif
  /** set the default solver parameters and perform memory allocation for VariationalInequality
      \param options the pointer to the options to set
      \param solverId the identifier of the solver
  */
  int variationalInequality_setDefaultSolverOptions(SolverOptions* options, int solverId);

  /**Extra Gradient solver forvariational inequality problem based on the De Saxce Formulation
      \param problem the variational inequality problem to solve
      \param x global vector (n), in-out parameter
      \param w global vector (n), in-out parameters
      \param info return 0 if the solution is found
      \param options the solver options :
      iparam[0] : Maximum iteration number
      iparam[1] : Choice of the line search
           0 : Armijo rule with Khotbotov ratio (default)
           1 : Armijo rule with Solodov.Tseng ratio
      iparam[2] : bool activate the update in the loop (0:false default choice)
      iparam[3] : use rho_k * tau * min(1.0,a2/(rho_k*a1)) to decrease rho; commented in the code

      dparam[3] : rho  parameter.
         If rho >0, then self-adaptive (Armijo like) procedure.
         If rho <0, then constant rho parameter  (rho <-- -rho)
      Adaptive step-size parameters:
      Adaptive step-size parameters:
      dparam[4] = 2/3.0;  tau
      dparam[5] = 3.0/2.0;  tauinv
      dparam[6] = 0.9;   L
      dparam[7] = 0.3;   Lmin
  */
  void variationalInequality_ExtraGradient(VariationalInequality* problem, double *x, double *w, int* info, SolverOptions* options);

  /** set the default solver parameters and perform memory allocation for EG
    \param options the pointer to the array of options to set
  */

  int variationalInequality_ExtraGradient_setDefaultSolverOptions(SolverOptions* options);



  /** Fixed Point Projection solver for variational inequality problem based on the De Saxce Formulation
      \param problem the variational inequality problem to solve
      \param x global vector (n), in-out parameter
      \param w global vector (n), in-out parameters
      \param info return 0 if the solution is found
      \param options the solver options :
      iparam[0] : Maximum iteration number
      iparam[1] : Choice of the line search
           0 : Armijo rule with Khotbotov ratio (default)
           1 : Armijo rule with Solodov.Tseng ratio
           2 : Armijo rule with Han.Sun ratio
      iparam[2] : bool activate the update in the loop (0:false default choice)
      iparam[3] : use rho_k * tau * min(1.0,a2/(rho_k*a1)) to decrease rho; commented in the code
      dparam[3] : rho parameter.
         If rho >0, then self-adaptive (Armijo like) procedure.
         If rho <0, then constant rho parameter  (rho <-- -rho)
      Adaptive step-size parameters:
      dparam[4] = 2/3.0;  tau
      dparam[5] = 3.0/2.0;  tauinv
      dparam[6] = 0.9;   L
      dparam[7] = 0.3;   Lmin

  */
  void variationalInequality_FixedPointProjection(VariationalInequality* problem, double *x, double *w, int* info, SolverOptions* options);


  /** set the default solver parameters and perform memory allocation for EG
    \param options the pointer to the array of options to set
  */
  int variationalInequality_FixedPointProjection_setDefaultSolverOptions(SolverOptions* options);


  /** Hyperplane Projection solver for variational inequality problem based on the De Saxce Formulation
      \param problem the variational inequality problem to solve
      \param x global vector (n), in-out parameter
      \param w global vector (n), in-out parameters
      \param info return 0 if the solution is found
      \param options the solver options :
      iparam[0] : Maximum iteration number
      dparam[3] : rho >0
  */
  void variationalInequality_HyperplaneProjection(VariationalInequality* problem, double *x, double *w, int* info, SolverOptions* options);


  /** set the default solver parameters and perform memory allocation for EG
    \param options the pointer to the array of options to set
  */
  int variationalInequality_HyperplaneProjection_setDefaultSolverOptions(SolverOptions* options);


  /** VI Solver based on a merit function minimization with a line-search type algorithm
   * \param problem the variational inequality problem to solve
   * \param[in,out] x as input, the initial guess; as output the solution if
   * the algorithm is successful
   * \param[in,out] F value of the function
   * \param info 0 if a solution is found
   * \param options the solver options
   */
  void variationalInequality_box_newton_QiLSA(VariationalInequality* problem, double *x, double *F, int* info, SolverOptions* options);

  /** VI Solver based on the Newton-Josephy method globalized with a line-search type algorithm using the Qi merit function.
   * \param problem the variational inequality problem to solve
   * \param[in,out] x as input, the initial guess; as output the solution if
   * the algorithm is successful
   * \param[in,out] F value of the function
   * \param info 0 if a solution is found
   * \param options the solver options
   */
  void vi_box_AVI_LSA(VariationalInequality* problem, double *x, double *F, int* info, SolverOptions* options);

  /** Extra SolverOptions settings
   * \param options the struct to modify
   */
  void vi_box_AVI_extra_SolverOptions(SolverOptions* options);

  /** free the specific solverData
   * \param options the struct to modify
   */
  void vi_box_AVI_free_solverData(SolverOptions* options);

  /** Get the set from the VariationalInequality
   * \param problem the VI
   */
  void* vi_get_set(void* problem);

  /** Solver for box constrainted VI using PATH.
   * \param problem the variational inequality problem to solve
   * \param[in,out] z as input, the initial guess; as output the solution if
   * the algorithm is successful
   * \param[in,out] F value of the function
   * \param info 0 if a solution is found
   * \param options the solver options
   */
  void vi_box_path(VariationalInequality* problem, double *z, double* F, int *info , SolverOptions* options);

  /**  set the default solver parameters and perform memory allocation for a VI
   * solver
   * \param options the SolverOptions to set
   * \param solverId the id of the solver
   */
  int variationalInequality_common_setDefaultSolverOptions(SolverOptions* options, int solverId);



  /** Check for trivial solution in the variational inequality problem
      \param problem VariationalInequality*  the problem
      \param x global vector (n), in-out parameter
      \param fx global vector (n), in-out parameters
      \param options the pointer to the array of options to set
      \return info  =0 if a trivial solution has been found, else = -1
  */
  int checkTrivialCase_vi(VariationalInequality* problem , double* x, double* fx, SolverOptions* options);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
