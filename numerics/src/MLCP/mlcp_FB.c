/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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
#include "mlcp_FB.h"  // for mlcp_FB_getNbDWork

#include <stdio.h>   // for printf, fprintf, stderr
#include <stdlib.h>  // for exit

#include "FischerBurmeister.h"                  // for jacobianPhi_Mixed_FB
#include "MLCP_Solvers.h"                       // for mixedLinearComplement...
#include "MixedLinearComplementarityProblem.h"  // for MixedLinearComplement...
#include "NonSmoothNewton.h"                    // for NewtonFunctionPtr
#include "NonSmoothNewtonNeighbour.h"           // for NSNN_reset, nonSmooth...
#include "NumericsFwd.h"                        // for MixedLinearComplement...
#include "NumericsMatrix.h"                     // for NM_gemv, NumericsMatrix
#include "SiconosBlas.h"                        // for cblas_dcopy
#include "SolverOptions.h"                      // for SolverOptions
#include "numerics_verbose.h"                   // for verbose

static int sN = 0;
static int sM = 0;
static MixedLinearComplementarityProblem* sProblem;
static double* sFz = 0;
static double sMaxError = 0;

static void computeFz(double* z);
static void F_MCPFischerBurmeister(int size, double* z, double* FBz, int a);
static void jacobianF_MCPFischerBurmeister(int size, double* z, double* jacobianFMatrix,
                                           int a);

/*
Warning: this function requires MLCP with M and q, not (A,B,C,D).
The input structure MixedLinearComplementarityProblem is supposed to fit with this form.
*/
int mlcp_FB_getNbIWork(MixedLinearComplementarityProblem* problem, SolverOptions* options) {
  return nonSmoothNewtonNeigh_getNbIWork(problem->n, problem->m);
}
int mlcp_FB_getNbDWork(MixedLinearComplementarityProblem* problem, SolverOptions* options) {
  return problem->n + problem->m + nonSmoothNewtonNeigh_getNbDWork(problem->n, problem->m);
}

void computeFz(double* z) {
  int incx = 1, incy = 1;
  int size = sN + sM;
  // F(z)=Mz+q
  cblas_dcopy(size, sProblem->q, incx, sFz, incy);
  NM_gemv(1.0, sProblem->M, z, 1.0, sFz);
}

void F_MCPFischerBurmeister(int size, double* z, double* FBz, int a) {
  computeFz(z);
  phi_Mixed_FB(sN, sM, z, sFz, FBz);
}

/** writes \f$ \nabla_z F(z) \f$  using MLCP formulation and the Fischer-Burmeister function.
 */
void jacobianF_MCPFischerBurmeister(int size, double* z, double* jacobianFMatrix, int a) {
  computeFz(z);
  jacobianPhi_Mixed_FB(sN, sM, z, sFz, sProblem->M->matrix0, jacobianFMatrix);
}

void mlcp_FB_init(MixedLinearComplementarityProblem* problem, SolverOptions* options) {
  /*
     Initialize solver (Connect F and its jacobian, set local size ...) according to the chosen
     formulation.
  */
  sProblem = problem;
  sN = problem->n;
  sM = problem->m;
  sFz = options->dWork;
  double* last =
      nonSmoothNewtonNeighInitMemory(sN + sM, options->dWork + sN + sM, options->iWork);
  double* tlast = options->dWork + mlcp_FB_getNbDWork(problem, options);
  if (last > tlast) {
    printf("internal error.");
    exit(1);
  }
}

void mlcp_FB_reset() {
  NSNN_reset();
  /*free(sFz) ;*/
}

void mlcp_FB(MixedLinearComplementarityProblem* problem, double* z, double* w, int* info,
             SolverOptions* options) {
  *info = 1;

  NewtonFunctionPtr F = &F_MCPFischerBurmeister;
  NewtonFunctionPtr jacobianF = &jacobianF_MCPFischerBurmeister;
  double err;
  double tol = options->dparam[SICONOS_DPARAM_TOL];
  int i;
  /*only for debug
  double * zz = (double *)malloc((sN+sM)*sizeof(double));
  memcpy(zz,z,(sN+sM)*sizeof(double));*/

  *info = nonSmoothNewtonNeigh(sN + sM, z, &F, &jacobianF, options->iparam, options->dparam);
  if (*info > 0) {
    numerics_printf_verbose(0,
                            "mlcp_FB failed, reached max. number of iterations without "
                            "convergence. Residual = %f\n",
                            options->dparam[SICONOS_DPARAM_RESIDU]);
    /*ONLY FOR DEBUG
      mixedLinearComplementarity_display(problem);
    printf("with z init;\n");
    for (i=0;i<sN+sM;i++)
    printf("%.32e \n",zz[i]);
    exit(1);*/
  }
  /*  free(zz);*/
  mlcp_compute_error(problem, z, w, tol, &err);
  for (i = 0; i < sM; i++) {
    if (z[sN + i] > w[sN + i]) w[sN + i] = 0;
  }

  if (err > sMaxError) sMaxError = err;
  numerics_printf_verbose(
      0, "mlcp_FB : succeed in %i iterations, error %14.8e   and max error  %14.8e \n",
      options->iparam[SICONOS_IPARAM_ITER_DONE], err, sMaxError);
  options->dparam[SICONOS_DPARAM_RESIDU] = err;

  return;
}
