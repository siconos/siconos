/* Siconos-Numerics, Copyright INRIA 2005-2014
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


#include "SiconosSets.h"
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
      for (unsigned i = 0; i < n; ++i)
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
      for (unsigned i = 0; i < n; ++i)
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
    default:
      printf("free_siconos_set :: not implemented for set of type %d\n", set_id);
      exit(EXIT_FAILURE);
  }
}

void free_box(box_constraints* b)
{
  assert(b);
  if (b->lb)
  {
    free(b->lb);
    b->lb = NULL;
  }
  if (b->ub)
  {
    free(b->ub);
    b->ub = NULL;
  }
}

void free_polyhedron(polyhedron* poly)
{
  assert(poly);
  if (poly->H)
  {
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
    free(poly->Heq);
    poly->Heq = NULL;
  }
  if (poly->Keq)
  {
    free(poly->Keq);
    poly->Keq = NULL;
  }
}
