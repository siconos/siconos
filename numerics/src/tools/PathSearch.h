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
/*!\file PathSearch.h
 * \brief Path search related functions and data
 */

#ifndef PATHSEARCH_H
#define PATHSEARCH_H

#include "NMS.h"             // for NMS_data
#include "Newton_methods.h"  // for functions_LSA
#include "NumericsFwd.h"     // for SolverOptions
#include "SiconosConfig.h" // for BUILD_AS_CPP // IWYU pragma: keep

/** struct ncp_pathsearch_data NCP_PathSearch.h
 * solver specific data
 */
typedef struct
{
  NMS_data* data_NMS; /**< struct for the NMS scheme */
  functions_LSA* lsa_functions; /**< functions for the search */
} pathsearch_data;

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** set some default value for the solver option when the path search
   * algorithm is used
   * \param options the structure to be modified
   */
  void pathsearch_set_default(SolverOptions* options);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
