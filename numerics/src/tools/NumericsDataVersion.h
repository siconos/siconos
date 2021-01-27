/* Siconos, Copyright INRIA 2005-2020.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
#ifndef NumericsDataVersion_h
#define NumericsDataVersion_h

#include <limits.h>
#include <stdint.h>
#include <assert.h>

typedef uint64_t version_t;

/** \struct DataVersioning data */
typedef struct
{
  version_t number;
} NumericsDataVersion;

static inline version_t NDV_value(const NumericsDataVersion* v)
{
  return v->number;
}

static inline void NDV_set_value(NumericsDataVersion* v, version_t value)
{
  v->number = value;
}

static inline void NDV_reset(NumericsDataVersion* v)
{
  v->number = 0;
};

static inline void NDV_inc(NumericsDataVersion* v)
{
  assert (v->number < UINT64_MAX);
  v->number += 1;
};

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif


