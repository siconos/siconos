/* Siconos-Numerics, Copyright INRIA 2005-2011.
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

#include <stdio.h>
#include <stdlib.h>
#include "NonSmoothDrivers.h"
#include "mcp_cst.h"

static double M[4] = {2.0, 1.0, 1.0, 2.0};
static double q[4] = { -5.0, -6.0};

void testF(int size, double *z, double * F)
{
  /* printf("call to MCP function F(z) ...\n");   */
  /* for (int i =0 ; i <size; i++) */
  /* { */
  /*   printf("z[%i]= %lf\t",i,z[i]); */
  /* } */
  /* printf("\n"); */
  for (int i = 0 ; i < size; i++)
  {
    F[i] = q[i];
    for (int j = 0 ; j < size ;  j++)
    {
      F[i] += M[i + j * size] * z[j] ;
    }
  }
  /* for (int i =0 ; i <size; i++) */
  /* { */
  /*   printf("F[%i]= %lf\t",i,F[i]); */
  /* } */
  /* printf("\n"); */
  /* printf("End call to MCP function F(z) ...\n");   */
}

void testNablaF(int size, double *z, double *nablaF)
{
  /* printf("call to MCP function nablaF(z) ...\n"); */

  for (int i = 0 ; i < size; i++)
  {
    for (int j = 0 ; j < size ;  j++)
    {
      nablaF[i + j * size] = M[i + j * size];
    }
  }


}

int main(void)
{
  printf(" Start tests for MCP solvers.\n");

  int info = 0 ;

  /* Set solver options */
  SolverOptions options;

  /* FB solver */
  options.solverId = SICONOS_MCP_FB;



  NumericsOptions global_options;
  global_options.verboseMode = 1;


  /* Create a MixedComplementarityProblem */
  MixedComplementarityProblem* problem = (MixedComplementarityProblem *)malloc(sizeof(MixedComplementarityProblem));

  problem->sizeEqualities = 1;
  problem->sizeInequalities = 1;
  problem->computeFmcp = &testF ;
  problem->computeNablaFmcp = &testNablaF ;
  problem->Fmcp = NULL;
  problem->nablaFmcp = NULL;




  mixedComplementarity_setDefaultSolverOptions(problem, &options);

  int size = problem->sizeEqualities + problem->sizeInequalities ;
  double * z = malloc(size * sizeof(double));
  double * w = malloc(size * sizeof(double));
  for (int i = 0 ; i < size; i++)
  {
    z[i] = 0.0;
    w[i] = 0.0;
  }

  options.dparam[0] = 1e-10;
  options.iparam[0] = 20;


  /* Initialize the solver */
  mcp_driver_init(problem, &options) ;
  info = mcp_driver(problem, z , w,  &options, &global_options);
  mcp_driver_reset(problem, &options) ;
  /// TODO : write a real test ... ////


  for (int i = 0 ; i < size; i++)
  {
    printf("z[%i]= %lf\t", i, z[i]);
  }
  printf("\n");
  for (int i = 0 ; i < size; i++)
  {
    printf("w[%i]= %lf\t", i, w[i]);
  }
  printf("\n");
  free(options.iparam);
  free(options.dparam);
  free(z);
  free(w);
  free(problem);
  printf("End of MCP solvers test. \n");

  return info;
}
