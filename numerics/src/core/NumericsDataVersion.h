// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef NumericsDataVersion_h
#define NumericsDataVersion_h

#include <assert.h>
#include <limits.h>
#include <stdint.h>

/** @brief this is a test */
typedef uint64_t version_t;

/** \struct DataVersioning data */
typedef struct {
  version_t number;
} NumericsDataVersion;

static inline version_t NDV_value(const NumericsDataVersion* v) { return v->number; }

static inline void NDV_set_value(NumericsDataVersion* v, version_t value) {
  v->number = value;
}

static inline void NDV_reset(NumericsDataVersion* v) { v->number = 0; }

void NDV_inc(NumericsDataVersion* v);

#if defined(__cplusplus)
extern "C" {
#endif

#if defined(__cplusplus)
}
#endif

#endif
