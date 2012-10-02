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

#ifndef FrictionContact3D_compute_error_H
#define FrictionContact3D_compute_error_H

/*!\file FrictionContact3D_compute_error.h
  \brief functions related to error computation for friction-contact problems

  \author Vincent Acary, 26/05/2008

*/

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** Error computation for friction-contact 3D problem
      \param problem the structure which defines the friction-contact problem
      \param z vector
      \param w vector
      \param tolerance value for error computation
      \param options
      \param[in,out] error value
   */
  int FrictionContact3D_compute_error(FrictionContactProblem* problem, double *z , double *w, double tolerance, SolverOptions * options, double * error);

  /** Error computation for one friction-contact 3D problem
      \param z vector
      \param w vector
      \param mu coeficient of friction
      \param[in,out] error value
   */
  void FrictionContact3D_unitary_compute_and_add_error(double z[3] , double w[3], double mu, double * error);
  /** Error computation for friction-contact 3D problem
      \param problem the structure which defines the friction-contact problem
      \param z vector
      \param w vector
      \param options
      \param tolerance value for error computation
      \param[in,out] error value
   */
  int FrictionContact3D_compute_error_velocity(FrictionContactProblem* problem, double *z , double *w, double tolerance, SolverOptions * options, double * error);



  /** Error computation for friction-contact 3D problem with Tresca Friction
      \param problem the structure which defines the friction-contact problem
      \param z vector
      \param w vector
      \param tolerance value for error computation
      \param options
      \param[in,out] error value
   */
  int FrictionContact3D_Tresca_compute_error(FrictionContactProblem* problem, double *z , double *w, double tolerance, SolverOptions * options, double * error);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
