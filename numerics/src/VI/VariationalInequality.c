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
#include "VariationalInequality.h"
#include "numerics_verbose.h"
#include "NumericsMatrix.h"


void variationalInequality_display(VariationalInequality* vi)
{

  assert(vi);

}

int variationalInequality_printInFile(VariationalInequality*  vi, FILE* file)
{
  if (! vi)
  {
    fprintf(stderr, "Numerics, VariationalInequality printInFile failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }

  return 0;
}

int variationalInequality_newFromFile(VariationalInequality* vi, FILE* file)
{
  return 0;
}

void freeVariationalInequalityProblem(VariationalInequality* vi)
{
  if (vi->nabla_F)
  {
    NM_free(vi->nabla_F);
    free(vi->nabla_F);
    vi->nabla_F = NULL;
  }
  free(vi);
}

void variationalInequality_clear(VariationalInequality* vi)
{
  vi->size = 0;
  vi->env = NULL;
  vi->F = NULL;
  vi->compute_nabla_F = NULL;
  vi->ProjectionOnX = NULL;
  vi->normVI = 0.;
  vi->istheNormVIset =0.;
  vi->set = NULL;
  vi->nabla_F = NULL;
}

VariationalInequality* variationalInequality_new(int size)
{
  VariationalInequality* fvi = (VariationalInequality*) malloc(sizeof(VariationalInequality)); 
  variationalInequality_clear(fvi);
  fvi->size = size;

  return fvi;
}

VariationalInequality* newVI(void)
{
  VariationalInequality* vi = (VariationalInequality*) malloc(sizeof(VariationalInequality));
  variationalInequality_clear(vi);

  return vi;

}
