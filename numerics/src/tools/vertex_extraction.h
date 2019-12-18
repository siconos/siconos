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

/*!\file vertex_extraction.h
 * \brief interface to various LP solvers
 */

#ifndef vertex_extraction_h
#define vertex_extraction_h

#include "SiconosConfig.h" // for BUILD_AS_CPP // IWYU pragma: keep
#include "SiconosLapack.h"  // for lapack_int
#include "SiconosSets.h"    // for polyhedron

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  void siconos_find_vertex(const polyhedron* P, unsigned size, lapack_int* basis);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif 
