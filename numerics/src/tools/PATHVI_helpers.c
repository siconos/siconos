/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2017 INRIA.
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

#include "PATHVI_helpers.h"

#ifdef HAVE_PATHVI

#include "numerics_verbose.h"
#include <string.h>
#include <stdio.h>

int pathvi_get_z(struct vi_desc *desc, double *z)
{
  SN_generic_pathvi_env* env = vi_desc_get_controller(desc);
  memcpy(z, env->z, env->n * sizeof(double));

  return 0;
}

int pathvi_set_z(struct vi_desc *desc, double *z)
{
  SN_generic_pathvi_env* env = vi_desc_get_controller(desc);
  memcpy(env->z, z, env->n * sizeof(double));

  return 0;
}

int pathvi_get_F(struct vi_desc *desc, double *F)
{
  SN_generic_pathvi_env* env = vi_desc_get_controller(desc);
  memcpy(F, env->F, env->n * sizeof(double));

  return 0;
}

int pathvi_set_F(struct vi_desc *desc, double *F)
{
  SN_generic_pathvi_env* env = vi_desc_get_controller(desc);
  memcpy(env->F, F, env->n * sizeof(double));

  return 0;
}

int pathvi_get_lambda(struct vi_desc *desc, double *lambda)
{
  SN_generic_pathvi_env* env = vi_desc_get_controller(desc);
  memcpy(lambda, env->lambda, env->m * sizeof(double));

  return 0;
}

int pathvi_set_lambda(struct vi_desc *desc, double *lambda)
{
  SN_generic_pathvi_env* env = vi_desc_get_controller(desc);
  memcpy(env->lambda, lambda, env->m * sizeof(double));

  return 0;
}

int pathvi_get_row_name(struct vi_desc *desc, int i, char *name, int len)
{
  snprintf(name, len, "r%d", i);
  return 0;
}

int pathvi_get_col_name(struct vi_desc *desc, int j, char *name, int len)
{
  snprintf(name, len, "c%d", j);
  return 0;
}

void pathvi_print(unsigned mode, const char *buf)
{
  numerics_printf(buf);
}

#endif /* HAVE_PATHVI */

