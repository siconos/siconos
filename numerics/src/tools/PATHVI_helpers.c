// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "PATHVI_helpers.h"

#ifdef HAVE_PATHVI

#include <stdio.h>
#include <string.h>

#include "numerics_errors.h"
#include "numerics_verbose.h"

int pathvi_get_z(struct vi_desc* desc, double* z) {
  SN_generic_pathvi_env* env = (SN_generic_pathvi_env*)vi_desc_get_controller(desc);
  memcpy(z, env->z, env->n * sizeof(double));

  return 0;
}

int pathvi_set_z(struct vi_desc* desc, double* z) {
  SN_generic_pathvi_env* env = (SN_generic_pathvi_env*)vi_desc_get_controller(desc);
  memcpy(env->z, z, env->n * sizeof(double));

  return 0;
}

int pathvi_get_F(struct vi_desc* desc, double* F) {
  SN_generic_pathvi_env* env = (SN_generic_pathvi_env*)vi_desc_get_controller(desc);
  memcpy(F, env->F, env->n * sizeof(double));

  return 0;
}

int pathvi_set_F(struct vi_desc* desc, double* F) {
  SN_generic_pathvi_env* env = (SN_generic_pathvi_env*)vi_desc_get_controller(desc);
  memcpy(env->F, F, env->n * sizeof(double));

  return 0;
}

int pathvi_get_lambda(struct vi_desc* desc, double* lambda) {
  SN_generic_pathvi_env* env = (SN_generic_pathvi_env*)vi_desc_get_controller(desc);
  memcpy(lambda, env->lambda, env->m * sizeof(double));

  return 0;
}

int pathvi_set_lambda(struct vi_desc* desc, double* lambda) {
  SN_generic_pathvi_env* env = (SN_generic_pathvi_env*)vi_desc_get_controller(desc);
  memcpy(env->lambda, lambda, env->m * sizeof(double));

  return 0;
}

int pathvi_get_row_name(struct vi_desc* desc, int i, char* name, int len) {
  snprintf(name, len, "r%d", i);
  return 0;
}

int pathvi_get_col_name(struct vi_desc* desc, int j, char* name, int len) {
  snprintf(name, len, "c%d", j);
  return 0;
}

void pathvi_print(unsigned mode, const char* buf) { numerics_printf(buf); }

#endif /* HAVE_PATHVI */
