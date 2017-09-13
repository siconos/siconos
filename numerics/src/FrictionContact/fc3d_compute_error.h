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

#ifndef fc3d_compute_error_H
#define fc3d_compute_error_H

/*!\file fc3d_compute_error.h
  \brief functions related to error computation for friction-contact problems

  \author Vincent Acary, 26/05/2008

*/

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** Error computation (using the normal map residual) for friction-contact 3D problem
      \param problem the structure which defines the friction-contact problem
      \param z vector
      \param w vector
      \param tolerance value for error computation
      \param options
      \param norm norm of a vector (problem->q) for relative error
      \param[in,out] error value
      \return 0 if ok
   */
  int fc3d_compute_error(FrictionContactProblem* problem, double *z , double *w, double tolerance, SolverOptions * options, double norm, double * error);

  /** Error computation (using the normal map residual) for one friction-contact 3D problem
      \param r the reaction force
      \param u the local velocity
      \param mu coeficient of friction
      \param worktmp work vector
      \param[in,out] error value
   */
  void fc3d_unitary_compute_and_add_error(double r[3] , double u[3], double mu, double * error, double * worktmp);

  /** Error computation for a friction-contact 3D problem
      \param problem the structure which defines the friction-contact problem
      \param z vector
      \param w vector
      \param options
      \param tolerance value for error computation
      \param[in,out] error value
      \return 0 if ok
   */
  int fc3d_compute_error_velocity(FrictionContactProblem* problem, double *z , double *w, double tolerance, SolverOptions * options, double * error);

  /** Error computation for one friction-contact 3D problem
      \param z vector
      \param w vector
      \param R radius of the cylinder
      \param worktmp work vector
      \param[in,out] error value
   */
  void fc3d_Tresca_unitary_compute_and_add_error(double z[3] , double w[3], double R, double * error, double *worktmp);


  /** Error computation for friction-contact 3D problem with Tresca Friction
      \param problem the structure which defines the friction-contact problem
      \param z vector
      \param w vector
      \param tolerance value for error computation
      \param options
      \param[in,out] error value
      \return 0 if ok
   */
  int fc3d_Tresca_compute_error(FrictionContactProblem* problem, double *z , double *w, double tolerance, SolverOptions * options,  double norm, double * error);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
