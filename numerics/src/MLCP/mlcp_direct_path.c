/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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

#include "mlcp_direct_path.h"
#include "MLCP_Solvers.h"                       // for mixedLinearComplement...
#include "MixedLinearComplementarityProblem.h"  // for MixedLinearComplement...
#include "mlcp_direct.h"                        // for mlcp_direct_addConfig...

static int sN;
static int sM;


void mlcp_direct_path_init(MixedLinearComplementarityProblem* problem, SolverOptions* options)
{
  sN = problem->n;
  sM = problem->m;
  mlcp_direct_init(problem, options);
  //mlcp_path_init(problem, options);

}
void mlcp_direct_path_reset()
{
  mlcp_direct_reset();
  //mlcp_path_reset();
}

void mlcp_direct_path(MixedLinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options)
{
  /*First, try direct solver*/
  mlcp_direct(problem, z, w, info, options);
  if(*info)
  {
    /*solver direct failed, so run the path solver.*/
    mlcp_path(problem, z, w, info, options);
    if(!(*info))
    {
      /*       for (i=0;i<problem->n+problem->m;i++){ */
      /*  printf("w[%d]=%f z[%d]=%f\t",i,w[i],i,z[i]);  */
      /*       } */
      mlcp_direct_addConfigFromWSolution(problem, w + sN);
    }
  }
}

