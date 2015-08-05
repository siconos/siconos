/* Siconos-Numerics, Copyright INRIA 2005-2012.
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
#ifndef FRICTIONCONTACT3D_localFischerBurmeister_H
#define FRICTIONCONTACT3D_localFischerBurmeister_H

/*!\file FrictionContact3D_localFischerBurmeister.h

  \brief Typedef and functions declarations related to Fischer-Burmeister
  formulation for 3 dimension frictional contact problems in Local Coordinates.

  Subroutines used when the friction-contact 3D problem is written
  using Fischer-Burmeister formulation:

  [...]

  where M is an n by n matrix, q an n-dimensional vector, reaction an
  n-dimensional vector and velocity an n-dimensional vector.

  We consider a "global" (ie for several contacts) problem, used to
  initialize the static global variables.


 */
#include "NumericsConfig.h"
#include "FrictionContactProblem.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** The Fischer Burmeister function signature for a 3x3 block.
   */
  typedef void (*FischerBurmeisterFun3x3Ptr)(double* reaction,
                                             double* velocity,
                                             double mu,
                                             double* rho,
                                             double* F,
                                             double* A,
                                             double* B);

  /** Nonsmooth Newton solver based on the Fischer--Burmeister function for the
   * local (reduced) frictional contact problem in the dense form
   * \param problem the problem to solve in dense form
   * \param reaction solution and initial guess for reaction
   * \param velocity solution and initial guess for velocity
   * \param info returned info
   * \param options  the solver options
   */
  void frictionContact3D_localFischerBurmeister(
    FrictionContactProblem* problem,
    double *reaction,
    double *velocity,
    int *info,
    SolverOptions *options);

  /** Nonsmooth Newton solver based on the Fischer--Burmeister function for the
   * local (reduced) frictional contact problem in the sparse block form
   * \param problem the problem to solve in sparse block form
   * \param reaction solution and initial guess for reaction
   * \param velocity solution and initial guess for velocity
   * \param info returned info
   * \param options  the solver options
   */
  void frictionContact3D_sparseLocalFischerBurmeister(
    FrictionContactProblem* problem,
    double *reaction,
    double *velocity,
    int *info,
    SolverOptions *options);

  /** Init the nonsmooth Newton solver based on the Fischer--Burmeister function
   * for the local (reduced) frictional contact problem in the sparse block form
   * Mainly, initialize the sparse linear solver
   * \param options  the solver options
   */
  void frictionContact3D_sparseLocalFischerBurmeisterInit(
    SolverOptions *options);

  /** The Fischer & Burmeister function for several contacts.
      On each contact, the specified Fischer Burmeister function in iparam[9] is called.
      \param problemSize the number of contacts.
      \param computeACFun3x3 the block 3x3 Fischer & Burmeister function.
      \param reaction3D the reactions at each contact (size: 3 x problemSize)
      \param velocity3D the velocities at each contact (size: 3 x problemSize)
      \param mu the mu parameter (size : problemSize)
      \param rho3D the rho parameters (size : 3 x problemSize)
      \param output_blocklist3 the computed Fischer & Burmeister function (size : 3 x problemSize)
      \param output_blocklist3x3_1 the computed A part of gradient (size : 9 x problemSize)
      \param output_blocklist3x3_2 the computed B param of gradient (size : 9 x problemSize)
  */
  void frictionContact3D_FischerBurmeisterFunction(
    unsigned int problemSize,
    FischerBurmeisterFun3x3Ptr computeACFun3x3,
    double *reaction3D,
    double *velocity3D,
    double *mu,
    double *rho3D,
    double *output_blocklist3,
    double *output_blocklist3x3_1,
    double *output_blocklist3x3_2);

  /* Set the default solver options for the LOCALAC Solver
   * Some default values:
   * options.iparam[0] = 200 is the maximum number of iterations.
   * options.iparam[3] = 100000 is the nzmax parameter for sparse matrices.
   * options.dparam[0] = 1e-3 precision.
   * \param options  the solver options
   */
  int frictionContact3D_FischerBurmeister_setDefaultSolverOptions(
    SolverOptions* options);


  int frictionContact3D_FischerBurmeister_compute_error(
    FrictionContactProblem* problem,
    double *z , double *w, double tolerance,
    SolverOptions * options, double * error);


  void nonsmoothEqnFischerBurmeisterFun(void* arg,
                                        unsigned int problemSize,
                                        double* reaction,
                                        double* velocity,
                                        double* mu,
                                        double* rho,
                                        double* result,
                                        double* A,
                                        double* B);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
