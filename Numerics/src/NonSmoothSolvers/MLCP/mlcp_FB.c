/* Siconos-Numerics, Copyright INRIA 2005-2010.
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
#include "LA.h"
#include "NumericsConfig.h"
#include "MLCP_Solvers.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "NonSmoothNewtonNeighbour.h"

#include "MCP/MCP_FischerBurmeister.h"

static int sN = 0;
static int sM = 0;
static MixedLinearComplementarity_Problem* sProblem;
static double* sFz = 0;
int mixedLinearComplementarity_fb_setDefaultSolverOptions(MixedLinearComplementarity_Problem* problem, SolverOptions* pSolver)
{
  mixedLinearComplementarity_default_setDefaultSolverOptions(problem, pSolver);
  return 0;
}

/*
Warning: this function requires MLCP with M and q, not (A,B,C,D).
The input structure MixedLinearComplementarity_Problem is supposed to fit with this form.
*/
int mlcp_FB_getNbIWork(MixedLinearComplementarity_Problem* problem, SolverOptions* options)
{
  return nonSmoothNewtonNeigh_getNbIWork(problem->n, problem->m);

}
int mlcp_FB_getNbDWork(MixedLinearComplementarity_Problem* problem, SolverOptions* options)
{
  return problem->n + problem->m + nonSmoothNewtonNeigh_getNbDWork(problem->n, problem->m);
}


void computeFz(double* z)
{
  int incx = 1, incy = 1;
  int size = sN + sM;
  //F(z)=Mz+q
  DCOPY(size , sProblem->q , incx , sFz , incy);
  prodNumericsMatrix(size, size, 1.0, sProblem->M, z, 1.0, sFz);
}
void F_MCPFischerBurmeister(int size, double* z, double* FBz, int a)
{

  computeFz(z);
  phi_MCP_FB(sN, sM, z, sFz, FBz);
}

/** writes \f$ \nabla_z F(z) \f$  using MLCP formulation and the Fischer-Burmeister function.
 */
void jacobianF_MCPFischerBurmeister(int size, double* z, double* jacobianFMatrix, int a)
{

  computeFz(z);
  jacobianPhi_MCP_FB(sN, sM, z, sFz, sProblem->M->matrix0, jacobianFMatrix);
}



void mlcp_FB_init(MixedLinearComplementarity_Problem* problem, SolverOptions* options)
{

  /*
     Initialize solver (Connect F and its jacobian, set local size ...) according to the chosen formulation.
  */
  sProblem = problem;
  sN = problem->n;
  sM = problem->m;
  sFz = options->dWork;
  double * last = nonSmoothNewtonNeighInitMemory(sN + sM, options->dWork + sN + sM, options->iWork);
  double * tlast = options->dWork + mlcp_FB_getNbDWork(problem, options);
  if (last > tlast)
  {
    printf("internal error.");
    exit(1);
  }

}

void mlcp_FB_reset()
{
  NSNN_reset();
  /*free(sFz) ;*/
}


void mlcp_FB(MixedLinearComplementarity_Problem* problem, double *z, double *w, int *info, SolverOptions* options)
{
  *info = 1;

  NewtonFunctionPtr F = &F_MCPFischerBurmeister;
  NewtonFunctionPtr jacobianF = &jacobianF_MCPFischerBurmeister;
  double err;
  double tol = options->dparam[0];
  int i;
  /*only for debug
  double * zz = (double *)malloc((sN+sM)*sizeof(double));
  memcpy(zz,z,(sN+sM)*sizeof(double));*/


  *info = nonSmoothNewtonNeigh(sN + sM, z, &F, &jacobianF, options->iparam, options->dparam);
  if (*info > 0)
  {
    fprintf(stderr, "Numerics, mlcp_FB failed, reached max. number of iterations without convergence. Error = %f\n", options->dparam[1]);
    /*ONLY FOR DEBUG
      displayMLCP(problem);
    printf("with z init;\n");
    for (i=0;i<sN+sM;i++)
    printf("%.32e \n",zz[i]);
    exit(1);*/
  }
  /*  free(zz);*/
  mlcp_compute_error(problem, z, w, tol, &err);
  for (i = 0; i < sM; i++)
  {
    if (z[sN + i] > w[sN + i])
      w[sN + i] = 0;
  }

  if (verbose)
    printf("FB : MLCP Solved, error %10.7f.\n", err);

  return;
}
