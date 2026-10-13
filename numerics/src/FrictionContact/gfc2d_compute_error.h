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

#ifndef gfc2d_compute_error_H
#define gfc2d_compute_error_H

/*!\file gfc2d_compute_error.h
 * \brief Error computation for the 2D global friction contact problem.
 */

#include "NumericsFwd.h"  // for SolverOptions, GlobalFrictionContactProblem

#if defined(__cplusplus)
extern "C" {
#endif

/**
 * \brief Compute a relative error for the 2D global friction contact problem.
 *
 * The error measures the violation of the dynamics equation
 * \f$M v - H r - q = 0\f$.
 *
 * \param[in] problem the global friction contact problem
 * \param[in] reaction contact reaction (size: 2 * numberOfContacts)
 * \param[in] velocity local contact velocity (size: 2 * numberOfContacts)
 * \param[in] globalVelocity global velocity (size: M->size0)
 * \param[in] tolerance requested tolerance
 * \param[in,out] options solver options; may allocate/reallocate options->dWork
 * \param[in] norm_q norm of q used for relative scaling
 * \param[in] norm_b norm of b used for relative scaling
 * \param[out] error computed relative error
 * \return 0 on success, nonzero on failure
 */
int gfc2d_compute_error(GlobalFrictionContactProblem* problem, double* reaction,
                        double* velocity, double* globalVelocity, double tolerance,
                        SolverOptions* options, double norm_q, double norm_b, double* error);

#if defined(__cplusplus)
}
#endif

#endif
