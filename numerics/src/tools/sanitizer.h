/* Siconos-Numerics, Copyright INRIA 2005-2015
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


#ifndef SN_SANITIZER_H
#define SN_SANITIZER_H

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
