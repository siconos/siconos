/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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

#ifndef NM_MA57_h
#define NM_MA57_h

#include "NumericsFwd.h"  // for NumericsMatrix
#include "SiconosConfig.h" // for BUILD_AS_CPP, WITH_MA57 // IWYU pragma: keep

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

#ifdef WITH_MA57
#include <lbl.h>

  /** Free the config data for MA57
   * \param param p a pointer on the linear solver parameters
   */
  void NM_MA57_free(void* p);

#endif /* WITH_MA57 */

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
