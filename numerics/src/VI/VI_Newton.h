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

#ifndef VI_NEWTON_H
#define VI_NEWTON_H

/*!\file VI_Newton.h
 * \brief Functions for solving VI using Newton method
 */

#include "SiconosConfig.h"
#include "NumericsMatrix.h"

#if defined(__cplusplus)
#undef restrict
#define restrict __restrict
#endif

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  void VI_compute_F(void* data_opaque, double* x, double* F);
  void VI_compute_error_box(void* data_opaque, double* x, double* F, double* Jac_F_merit, double tol, double* err);
  void VI_compute_F_box_Qi(void* data_opaque, double* x, double* F, double* Fbox);
  void VI_compute_H_box_Qi(void* data_opaque, double* x, double* F, double* workV1, double* workV2, NumericsMatrix* H);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
