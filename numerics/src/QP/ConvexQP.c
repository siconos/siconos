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
#include <stdlib.h>
#include <assert.h>
#include "ConvexQP.h"
#include "numerics_verbose.h"
#include "NumericsMatrix.h"
#include "NumericsVector.h"


void convexQP_display(ConvexQP* cqp)
{
  assert(cqp);
  printf("ConvexQP Display :\n-------------\n");
  printf("size :%d \n", cqp->size);
  printf("m:%d \n", cqp->m);
  printf("M matrix:\n");
  if (cqp->M)
  {

    NM_display(cqp->M);
  }
  else
  {
    printf("cqp->M is NULL\n");
  }
  printf("q vector:\n");
  if (cqp->q)
  {
    NV_display(cqp->q, cqp->size);
  }
  else
  {
    printf("cqp->q is NULL\n");
  }
  printf("A matrix:\n");
  if (cqp->A)
    NM_display(cqp->A);
  else
    printf("cqp->A is NULL\n");
  printf("b vector:\n");
  if (cqp->b)
    NV_display(cqp->b, cqp->m);
  else
    printf("cqp->b is NULL\n");

}

int convexQP_printInFile(ConvexQP*  cqp, FILE* file)
{
  if (! cqp)
  {
    fprintf(stderr, "Numerics, ConvexQP printInFile failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }

  return 0;
}

int convexQP_newFromFile(ConvexQP* cqp, FILE* file)
{
  return 0;
}

void convexQP_free(ConvexQP* cqp)
{
  if (cqp->M)
  {
    NM_free(cqp->M);
    free(cqp->M);
    cqp->M = NULL;
  }
  if (cqp->q)
  {
    free(cqp->q);
    cqp->q = NULL;
  }
  if (cqp->A)
  {
    NM_free(cqp->A);
    free(cqp->A);
    cqp->A = NULL;
  }
  if (cqp->b)
  {
    free(cqp->b);
    cqp->b = NULL;
  }
  free(cqp);
}

void convexQP_clear(ConvexQP* cqp)
{
  assert(cqp);
  cqp->size = 0;
  cqp->m=0;
  cqp->env = NULL;
  cqp->M=NULL;
  cqp->q=NULL;
  cqp->A=NULL;
  cqp->b=NULL;
  cqp->ProjectionOnC = NULL;
  cqp->normConvexQP = 0.;
  cqp->istheNormConvexQPset =0.;
}

ConvexQP* convexQP_new(int size)
{
  ConvexQP* fcqp = (ConvexQP*) malloc(sizeof(ConvexQP)); 
  convexQP_clear(fcqp);
  fcqp->size = size;
  return fcqp;
}

