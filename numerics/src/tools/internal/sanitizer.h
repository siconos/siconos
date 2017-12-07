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


#ifndef SN_SANITIZER_H
#define SN_SANITIZER_H

/*!\file sanitizer.h
 * \brief sanitizer helpers
 */

#include "SiconosBlas.h"

#if defined(__has_feature) 
#if  __has_feature(memory_sanitizer)
#define MSAN_INIT_VAR(v, n) for (size_t ii = 0; ii < (size_t)n; ++ii) v[ii] = 0.;

#else

#define MSAN_INIT_VAR(v, n)

#endif

#else

#define MSAN_INIT_VAR(v, n)

#endif

/** wrapper around cblas_dcopy which initialize the destination in order to make msan happy
 * \param n length of the data
 * \param src source
 * \param inc_src increment on the source
 * \param dest destination
 * \param inc_src increment on the destination
 */
static inline void cblas_dcopy_msan(int n, double* src, int inc_src, double* dest, int inc_dest)
{
  MSAN_INIT_VAR(dest, n);
  cblas_dcopy(n, src, inc_src, dest, inc_dest);
}

#endif
