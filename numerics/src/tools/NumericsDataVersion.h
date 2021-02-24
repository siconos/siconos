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


