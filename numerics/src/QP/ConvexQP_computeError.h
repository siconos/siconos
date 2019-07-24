/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

#ifndef ConvexQP_compute_error_H
#define ConvexQP_compute_error_H

/*!\file ConvexQP_computeError.h
  \brief functions related to error computation for friction-contact 
*/

#ifdef __cplusplus
#undef restrict
#define restrict __restrict
#endif

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** Error computation for a reduced ConvexQP problem (A=I, b=0).
      This function requires dWork to point to
      at least 2*n+m double of allocated memory or it mallocs this memory.
   \param problem the structure which defines the ConvexQP problem
   \param z vector
   \param w vector
   \param tolerance value for error computation
   \param options solver options
   \param norm coeff to normalize error
   \param[in,out] error value
   \return 0 if ok
   */
  int convexQP_compute_error_reduced(ConvexQP* problem, double *z , double *w, double tolerance, SolverOptions * options, double norm, double * error);

  /** Error computation for a ConvexQP problem;
      this function requires dWork to point to
      at least 2*n+m double of allocated memory or it mallocs this memory.
      \param problem the structure which defines the ConvexQP problem
      \param z vector
      \param xi multiplier vector
      \param w vector (w = s A^T xi)
      \param u (=Az+b) constraints vector
      \param tolerance value for error computation
      \param scaling  parameter s applied on the multiplier xi
      \param options solver options
      \param norm coeff to normalize error
      \param[in,out] error value
      \return 0 if ok
   */
  int convexQP_compute_error(ConvexQP* problem,
                                 double *z , double *xi,
                                 double* w, double * u,
                                 double tolerance,
                                 double scaling,
                                 SolverOptions * options, double norm,
                                 double * error);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
