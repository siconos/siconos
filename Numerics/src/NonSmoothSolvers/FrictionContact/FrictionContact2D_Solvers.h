/* Siconos-Numerics, Copyright INRIA 2005-2011.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
#ifndef FrictionContact2DSolvers_H
#define FrictionContact2DSolvers_H

/*!\file FrictionContact2D_Solvers.h
  \brief Subroutines for the resolution of contact problems with friction (2-dimensional case).\n
  \author Nineb Sheherazade and Dubois Frederic.
  Last Modifications : Mathieu Renouf , Pascal Denoyelle, Franck Perignon
*/

/*! \page FC2DSolvers Friction-Contact 2D problems Solvers

This page gives an overview of the available solvers for friction-contact (2D) problems and their required parameters.

For each solver, the input argument are:
- a FrictionContactProblem
- the unknowns (reaction,velocity)2
- info, the termination value (0: convergence, >0 problem which depends on the solver)
- a SolverOptions structure, which handles iparam and dparam

\section fc2Dcpg CPG Solver

\bf function: pfc_2D_cpg()
\bf parameters:
  - iparam[0] (in), the maximum number of iterations allowed,
  - iparam[1] (out), the number of iterations performed by the algorithm.
  - dparam[0] (in), the tolerance required,
  - dparam[1] (out), the residu.

*/

#include "FrictionContactProblem.h"
#include "NumericsOptions.h"
#include "SolverOptions.h"
#include "Friction_cst.h"
#ifdef __cplusplus
extern "C"
{
#endif

  /** General interface to solvers for friction-contact 2D problem
      \param[in] problem the structure which handles the Friction-Contact problem
      \param[in,out] reaction global vector (n)
      \param[in,out] velocity global vector (n)
      \param[in,out] options structure used to define the solver(s) and their parameters
      \param[in] global_options for Numerics (verbose mode ...)
      \return result (0 if successful otherwise 1).
  */
  int frictionContact2D_driver(FrictionContactProblem* problem, double *reaction , double *velocity, SolverOptions* options, NumericsOptions* global_options);


  /** set the default solver parameters and perform memory allocation for FrictionContact3D
      \param options   the pointer to the options to set
      \param solverId  the identifier of the solver
  */
  int frictionContact2D_setDefaultSolverOptions(SolverOptions* options, int solverId);

  /**  cpg (conjugated projected gradient) solver for primal contact problems with friction (2D)
       \param[in]  problem the friction-contact problem
       \param[out] reaction vector
       \param[out] velocity vector
       \param[in,out] info termination value
       \param[in,out] options structure for options
  */
  void FrictionContact2D_cpg(FrictionContactProblem* problem , double *reaction , double *velocity , int *info, SolverOptions* options);

  /** set the default solver parameters and perform memory allocation for CPG
   \param  options SolverOptions * the pointer to the options to set
   */
  int frictionContact2D_cpg_setDefaultSolverOptions(SolverOptions* options);
  /**  Non Linear Gauss Seidel solver for primal contact problem with friction in 2D case.
       \param[in] problem the friction-contact problem
       \param[out] reaction vector
       \param[out] velocity vector
       \param[in,out] info termination value
       \param[in,out] options structure
  */
  void FrictionContact2D_nsgs(FrictionContactProblem* problem , double *reaction , double *velocity , int *info, SolverOptions* options);

  /** set the default solver parameters and perform memory allocation for LATIN
  \param options  the pointer to the options to set
  */
  int frictionContact2D_nsgs_setDefaultSolverOptions(SolverOptions* options);

  /**  latin solver for primal contact problem with friction in the 2D case.
       \param[in] problem the friction-contact problem
       \param[out] reaction global vector
       \param[out] velocity global vector
       \param[in,out] info termination value
       \param[in,out] options SolverOptions structure
       \author Nineb Sheherazade.
  */
  void FrictionContact2D_latin(FrictionContactProblem* problem , double *reaction , double *velocity , int *info, SolverOptions* options);

  /** set the default solver parameters and perform memory allocation for LATIN
  \param  options the pointer to the options to set
  */
  int frictionContact2D_latin_setDefaultSolverOptions(SolverOptions* options);

  /** FrictionContact2D_projc is a specific projection operator related to CPG (conjugated projected gradient) algorithm for primal contact problem with friction.\n
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
  void FrictionContact2D_projc(double* xi, int* n, int* statusi, double* p, double* fric, double *reaction, int *status);

  /** FrictionContact2D_projf is a specific projection operator related to CPG (conjugated projected gradient) algorithm
   *              for primal contact problem with friction.\n
   *
   *
   * \param[in] statusi  parameter which represents the status vector.
   * \param[in] n      parameter which represents the dimension of the system.
   * \param[in] y    parameter which contains the components of the residue or descent direction vector.
   * \param[in] fric   parameter which contains the friction coefficient.
   * \param[out] projf1 parameter which contains the projected residue or descent direction.
   *
   */
  void FrictionContact2D_projf(int* statusi, int* n , double *y , double *fric, double *projf1);

  /** */
  void frictionContact2D_sparse_nsgs(FrictionContactProblem* problem, double *z, double *w, int *info, SolverOptions* options) ;

  /** set the default solver parameters and perform memory allocation for NSGS
  \param options the pointer to the options to set
  */
  int frictionContact2D_sparse_nsgs_setDefaultSolverOptions(SolverOptions* options);


#ifdef __cplusplus
}
#endif

#endif
