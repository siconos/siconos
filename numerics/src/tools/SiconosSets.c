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


#include "SiconosSets.h"
#include "NumericsMatrix.h"
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

#ifdef __cplusplus
#undef restrict
#define restrict __restrict
#endif

void project_on_set(int n, double* restrict x, void* restrict set)
{
  assert(set);
  assert(x);

  assert(n > 0);

  int set_id = ((generic_set*)set)->id;
  switch(set_id)
  {
    case SICONOS_SET_POSITIVE_ORTHANT:
      for (int i = 0; i < n; ++i)
      {
        if (x[i] < 0.0)
          x[i] = 0.0;
      }
      break;
    case SICONOS_SET_BOX:
      {
      box_constraints* box = (box_constraints*) set;
      assert(box->lb);
      assert(box->ub);
      for (int i = 0; i < n; ++i)
        x[i] = x[i] >= box->lb[i] ? (x[i] <= box->ub[i] ? x[i] : box->ub[i]) : box->lb[i];
      break;
      }
    default:
      printf("project_on_set :: not implemented for set of type %d\n", set_id);
      exit(EXIT_FAILURE);
  }
}

void free_siconos_set(void* set)
{
  assert(set);

  int set_id = ((generic_set*)set)->id;
  switch(set_id)
  {
    case SICONOS_SET_POSITIVE_ORTHANT:
      break;
    case SICONOS_SET_BOX:
      free_box((box_constraints*) set);
      break;
    case SICONOS_SET_POLYHEDRON:
      free_polyhedron((polyhedron*) set);
      break;
    case SICONOS_SET_POLYHEDRON_UNIFIED:
      free_polyhedron_unified((polyhedron_unified*) set);
      break;
    default:
      printf("free_siconos_set :: not implemented for set of type %d\n", set_id);
      exit(EXIT_FAILURE);
  }
}

void free_box(box_constraints* b)
{
  assert(b);
  if (b->lb) { free(b->lb); b->lb = NULL; }
  if (b->ub) { free(b->ub); b->ub = NULL; }
}

void free_polyhedron(polyhedron* poly)
{
  assert(poly);
  if (poly->H)
  {
    NM_free(poly->H);
    free(poly->H);
    poly->H = NULL;
  }
  if (poly->K)
  {
    free(poly->K);
    poly->K = NULL;
  }
  if (poly->Heq)
  {
    NM_free(poly->Heq);
    free(poly->Heq);
    poly->Heq = NULL;
  }
  if (poly->Keq)
  {
    free(poly->Keq);
    poly->Keq = NULL;
  }
}

void free_polyhedron_unified(polyhedron_unified* poly)
{
  assert(poly);
  if (poly->A)
  {
    NM_free(poly->A);
    free(poly->A);
    poly->A = NULL;
  }
  if (poly->b) { free(poly->b); poly->b = NULL; }
  if (poly->type) { free(poly->type); poly->type = NULL; }

}
