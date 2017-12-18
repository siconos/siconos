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

#include "NonlinearComplementarityProblem.h"
#include "NumericsMatrix.h"


void freeNCP(NonlinearComplementarityProblem* ncp)
{
  if (ncp->nabla_F)
  {
    NM_free(ncp->nabla_F);
    free(ncp->nabla_F);
    ncp->nabla_F = NULL;
  }

  free(ncp);
}

NonlinearComplementarityProblem* newNCP(void)
{
  NonlinearComplementarityProblem* ncp = (NonlinearComplementarityProblem*) malloc(sizeof(NonlinearComplementarityProblem));

  ncp->n = 0;
  ncp->compute_F = NULL;
  ncp->compute_nabla_F = NULL;
  ncp->nabla_F = NULL;
  ncp->env = NULL;

  return ncp;
}
