/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
#include <stdio.h>
#include <math.h>
#include "Standalone_Path.h"
#include "NCP_Path.h"

static int fill_structure; /* Do we need to fill in the structure of    */
/* the Jacobian?                             */

int NCP_Path(int n, double* z, FuncEvalPtr F, JacEvalPtr jacobianF, int* iparam, double* dparam)
{
  /* Lower bounds on the variables = 0 for NCP */
  double *lb  = (double *)malloc(sizeof(double) * n);
  /* Upper bounds on the variables = +inf for NCP */
  double *ub = (double *)malloc(sizeof(double) * n);
  /* Function evaluation */
  double *f = (double *)malloc(sizeof(double) * n);;
  /* Number of nonzeros */
  int nnz = n * n;
  /* Termination status from PATH */
  int status;
  int j;

  /* Connect F and its jacobian to input functions */
  setFuncEval(F);
  setJacEval(jacobianF);

  /**************************************/
  /* Fill in the lower and upper bounds */
  /**************************************/

  for(j = 0; j < n; j++)
  {
    lb[j] = 0;
    ub[j] = 1e20;
  }

  /************************************************************************/
  /* Set fill_structure to true (we ALWAYS need to fill in the structure  */
  /* of the jacobian at least once per call to PATH).                     */
  /************************************************************************/

  fill_structure = 1;

  /************************************************************************/
  /* Call PATH.                                                           */
  /************************************************************************/

  pathMain(n, nnz, &status, z, f, lb, ub);

  /************************************************************************/
  /* Deallocate memory.                                                   */
  /************************************************************************/

  free(lb);
  free(ub);
  free(f);
  return 0;
}

