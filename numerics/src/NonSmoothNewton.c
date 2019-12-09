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
#include "numerics_verbose.h"
#include "NonSmoothNewton.h"
#include "SolverOptions.h"
#include "SiconosLapack.h"
#include "NumericsMatrix.h"


#include "Newton_methods.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"

/* #define DEBUG_MESSAGES */
#include "debug.h"


void linesearch_Armijo(int n, double *z, double* dir, double psi_k,
                       double descentCondition, NewtonFunctionPtr* phi)
{
  double * phiVector = (double*)malloc(n * sizeof(*phiVector));
  if (phiVector == NULL)
  {
    fprintf(stderr, "NonSmoothNewton: linesearch_Armijo, memory allocation failed for phiVector\n");
    exit(EXIT_FAILURE);
  }

  /* IN :
     psi_k (merit function for current iteration)
     jacobian_psi_k (jacobian of the merit function)
     dk: descent direction

     OUT: tk, z
  */

  double sigma = 1e-4;
  double tk = 1;
  int incx = 1, incy = 1;
  double merit, merit_k;
  double tmin = 1e-12;

  /* z1 = z0 + dir */
  cblas_daxpy(n , 1.0 , dir , incx , z , incy);

  while (tk > tmin)
  {
    numerics_printf_verbose(2,"Non Smooth Newton:\t\tlinesearch_Armijo. try tk = %e", tk);
    /* Computes merit function = 1/2*norm(phi(z_{k+1}))^2 */
    (*phi)(n, z, phiVector, 0);
    merit =  cblas_dnrm2(n, phiVector , incx);
    merit = 0.5 * merit * merit;
    merit_k = psi_k + sigma * tk * descentCondition;
    if (merit < merit_k) break;
    tk = tk * 0.5;
    /* Computes z_k+1 = z0 + tk.dir
    warning: (-tk) because we need to start from z0 at each step while z_k+1 is saved in place of z_k ...*/
    cblas_daxpy(n , -tk , dir , incx , z , incy);
  }
  free(phiVector);
  if (tk <= tmin)
    numerics_printf("Non Smooth Newton:\t\t linesearch_Armijo warning, resulting tk < tmin, linesearch stopped.");
  else
    numerics_printf_verbose(2,"Non Smooth Newton:\t\tlinesearch_Armijo succeeded with tk = %e", tk);

}

int nonSmoothNewton(
  int n,
  double* z,
  NewtonFunctionPtr* phi,
  NewtonFunctionPtr* jacobianPhi,
  SolverOptions * options)
{
  if (phi == NULL || jacobianPhi == NULL)
  {
    fprintf(stderr, "NonSmoothNewton error: phi or its jacobian function = NULL pointer.\n");
    exit(EXIT_FAILURE);
  }

  int * iparam = options->iparam;
  double * dparam = options->dparam;
  
  int itermax = iparam[SICONOS_IPARAM_MAX_ITER]; // maximum number of iterations allowed
  int niter = 0; // current iteration number
  double tolerance = dparam[SICONOS_DPARAM_TOL];

  numerics_printf("Non Smooth Newton: ============= Starting of Newton process =============");
  numerics_printf("   - tolerance: %14.7e ", tolerance);
  numerics_printf("   - maximum number of iterations: %i", itermax);

  int incx = 1;

  lapack_int infoDGESV;

  /* Memory allocation for phi and its jacobian */
  double * phiVector = (double*)malloc(n * sizeof(double));


  NumericsMatrix * H = NM_create(NM_DENSE, n, n);
  double * jacobianPhiMatrix = H->matrix0;

  /** Memory allocation for the merit function and its gradient */
  double psi;
  double * gradient_psi = (double*)malloc(n * sizeof(*gradient_psi));

  lapack_int* ipiv = (lapack_int *)malloc(n * sizeof(lapack_int));
  if (phiVector == NULL || jacobianPhiMatrix == NULL ||  gradient_psi == NULL || ipiv == NULL)
  {
    fprintf(stderr, "NonSmoothNewton, memory allocation failed.\n");
    exit(EXIT_FAILURE);
  }

  /** The algorithm is alg 4.1 of the paper of Kanzow and Kleinmichel,
   * "A new class of semismooth Newton-type methods for nonlinear complementarity problems",
   * in Computational Optimization and Applications, 11, 227-251 (1998).
   *   We try to keep the same notations
   */

  double rho = 1e-8;
  double descentCondition, criterion, norm_gradient_psi, normPhi;
  double p = 2.1;
  double terminationCriterion = 1;

  /** Iterations ... */
  while ((niter < itermax) && (terminationCriterion > tolerance))
  {
    ++niter;
    /** Computes phi and its jacobian */
    (*phi)(n, z, phiVector, 0);
    (*jacobianPhi)(n, z, jacobianPhiMatrix, 1);

    /* Computes the gradient of the merit function,
     * gradient_psi = transpose(jacobianPhiMatrix).phiVector */

    cblas_dgemv(CblasColMajor, CblasTrans, n, n, 1.0, jacobianPhiMatrix, n,
                phiVector, incx, 0.0, gradient_psi, incx);
    DEBUG_PRINTF("norm 1 of jacobianPhiMatrix = %e\n", NM_norm_1(H));

    norm_gradient_psi = cblas_dnrm2(n, gradient_psi, 1);

    /* Computes norm2(phi) */
    normPhi = cblas_dnrm2(n, phiVector, 1);
    DEBUG_PRINTF("norm of phiVector = %e\n",normPhi )

    /* Computes merit function */
    psi = 0.5 * normPhi * normPhi;

    /* Stops if the termination criterion is satisfied */
    if (options->iparam[SICONOS_IPARAM_STOPPING_CRITERION] == SICONOS_STOPPING_CRITERION_RESIDU)
    {
      terminationCriterion = normPhi;
    }
    else if (options->iparam[SICONOS_IPARAM_STOPPING_CRITERION] == SICONOS_STOPPING_CRITERION_STATIONARITY)
    {
      terminationCriterion = norm_gradient_psi;
    }
    else if (options->iparam[SICONOS_IPARAM_STOPPING_CRITERION] ==
             SICONOS_STOPPING_CRITERION_RESIDU_AND_STATIONARITY)
    {
      terminationCriterion = fmax(normPhi, norm_gradient_psi);
    }
      
    if (terminationCriterion < tolerance)
      break;

    /* Search direction calculation
     *
     * Find a solution dk of jacobianPhiMatrix.d = -phiVector.
     * dk is saved in phiVector.
    */
    cblas_dscal(n , -1.0 , phiVector, incx);
    DGESV(n, 1, jacobianPhiMatrix, n, ipiv, phiVector, n, &infoDGESV);

    /* descentCondition = gradient_psi^T dk */
    descentCondition = cblas_ddot(n, gradient_psi,  1,  phiVector, 1);

    /* Criterion to be satisfied: error < -rho*norm(dk)^p */
    criterion = cblas_dnrm2(n, phiVector, 1);
    criterion = -rho * pow(criterion, p);

    DEBUG_PRINTF("descentcondition = %e\n", descentCondition );
    DEBUG_PRINTF("criterion = %e\n", criterion );
    if (infoDGESV != 0 || descentCondition > criterion)
    {
      numerics_printf("Newton descent direction is not good. Use the gradient direction");
      /* If the linear system is not solved correctly or the descent condition
       * is not satisfied, we fall back to the gradient for the descent direction
       * dk = - gradient_psi (remind that dk is saved in phiVector) */
      cblas_dcopy(n, gradient_psi, 1, phiVector, 1);
      cblas_dscal(n , -1.0 , phiVector, incx);
    }

    /* Step-3 Line search: computes z_k+1 */
    linesearch_Armijo(n, z, phiVector, psi, descentCondition, phi);

    numerics_printf("Non Smooth Newton: iteration number %i, norm merit function = %e, norm grad. merit function = %14.7e .", niter, normPhi, norm_gradient_psi);
  }

  /* Total number of iterations */
  iparam[SICONOS_IPARAM_ITER_DONE] = niter;
  /* Final error */
  dparam[SICONOS_DPARAM_RESIDU] = terminationCriterion;

  /** Free memory*/
  free(phiVector);
  NM_clear(H);
  free(gradient_psi);
  free(ipiv);


  if (dparam[SICONOS_DPARAM_RESIDU] > tolerance)
    numerics_printf("Non Smooth Newton:  warning. no convergence after %i iterations" , niter);

  else
    numerics_printf("Non Smooth Newton: convergence after %i iterations" , niter);
  numerics_printf("Non Smooth Newton:  residual = : %e ", dparam[1]);


  if (dparam[SICONOS_DPARAM_RESIDU] > tolerance)
    return 1;
  else return 0;
}

int nonSmoothDirectNewton(
  int n,
  double* z,
  NewtonFunctionPtr* phi,
  NewtonFunctionPtr* jacobianPhi,
  SolverOptions * options)
{
  if (phi == NULL || jacobianPhi == NULL)
  {
    fprintf(stderr, "NonSmoothNewton error: phi or its jacobian function = NULL pointer.\n");
    exit(EXIT_FAILURE);
  }
  int * iparam = options->iparam;
  double * dparam = options->dparam;
  int itermax = iparam[SICONOS_IPARAM_MAX_ITER]; // maximum number of iterations allowed
  int niter = 0; // current iteration number
  double tolerance = dparam[SICONOS_DPARAM_TOL];
  numerics_printf("Non Smooth Newton: ============= Starting of Newton process =============");
  numerics_printf("   - tolerance: %14.7e ", tolerance);
  numerics_printf("   - maximum number of iterations: %i", itermax);

  int incx = 1;
  lapack_int infoDGESV = 0;
  int n2 = n*n;
  /* Memory allocation for phi and its jacobian */
  double * phiVector = (double*)malloc(n * sizeof(*phiVector));
  double *jacobianPhiMatrix = (double*)malloc(n2 * sizeof(*jacobianPhiMatrix));
  /** merit function and its jacobian */
  double *gradient_psi = (double*)malloc(n * sizeof(*gradient_psi));
  lapack_int* ipiv = (lapack_int *)malloc(n * sizeof(lapack_int));
  if (phiVector == NULL || jacobianPhiMatrix == NULL ||  gradient_psi == NULL || ipiv == NULL)
  {
    fprintf(stderr, "NonSmoothNewton, memory allocation failed.\n");
    exit(EXIT_FAILURE);
  }

  /** The algorithm is alg 4.1 of the paper of Kanzow and Kleinmichel, "A new class of semismooth Newton-type methods
   *  for nonlinear complementarity problems", in Computational Optimization and Applications, 11, 227-251 (1998).
   *
   *   We try to keep the same notations
   */

  double norm_gradient_psi;
  double terminationCriterion = 1;
  if (gradient_psi == NULL)
  {
    fprintf(stderr, "NonSmoothNewton, memory allocation failed for gradient_psi.\n");
    exit(EXIT_FAILURE);
  }

  /** Iterations ... */
  while ((niter < itermax) && (terminationCriterion > tolerance))
  {
    ++niter;
    /** Computes phi and its jacobian */
    (*phi)(n, z, phiVector, 0);
    (*jacobianPhi)(n, z, jacobianPhiMatrix, 1);
    /* Computes the jacobian of the merit function, gradient_psi = transpose(jacobianPhiMatrix).phiVector */
    cblas_dgemv(CblasColMajor,CblasTrans, n, n, 1.0, jacobianPhiMatrix, n, phiVector, incx, 0.0, gradient_psi, incx);
    norm_gradient_psi = cblas_dnrm2(n, gradient_psi, 1);

    /* Computes norm2(phi) */
    // normPhi = cblas_dnrm2(n, phiVector, 1);
    /* Computes merit function */
    //psi = 0.5 * normPhi * normPhi;

    /* Stops if the termination criterion is satisfied */
    terminationCriterion = norm_gradient_psi;
    if (terminationCriterion < tolerance)
      break;

    /* Search direction calculation
    Find a solution dk of jacobianPhiMatrix.d = -phiVector.
    dk is saved in phiVector.
    */
    cblas_dscal(n , -1.0 , phiVector, incx);
    DGESV(n, 1, jacobianPhiMatrix, n, ipiv, phiVector, n, &infoDGESV);

    double tk = -1;

    cblas_daxpy(n , tk , phiVector , 1 , z , 1);

    numerics_printf("Non Smooth Newton: iteration number %i, norm of the merit function = %14.7e .", niter, terminationCriterion);

  }

  /* Total number of iterations */
  iparam[SICONOS_IPARAM_ITER_DONE] = niter;
  /* Final error */
  dparam[SICONOS_DPARAM_RESIDU] = terminationCriterion;

  /** Free memory*/
  free(phiVector);
  free(jacobianPhiMatrix);
  free(gradient_psi);
  free(ipiv);

  if (dparam[SICONOS_DPARAM_RESIDU] > tolerance)
    numerics_printf("Non Smooth Newton:  warning. no convergence after %i iterations" , niter);

  else
    numerics_printf("Non Smooth Newton: convergence after %i iterations" , niter);
  numerics_printf("Non Smooth Newton:  residual = : %e ", dparam[SICONOS_DPARAM_RESIDU]);

  if (dparam[SICONOS_DPARAM_RESIDU] > tolerance)
    return 1;
  else return 0;
}


void nonSmoothNewton_setDefaultSolverOptions(SolverOptions* options)
{

  numerics_printf_verbose(1,"nonSmoothNewton_setDefaultSolverOptions");

  options->solverId = SICONOS_NONSMOOTH_NEWTON_LSA;
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 20;
  options->dSize = 20;
  options->iparam = (int *)calloc(options->iSize, sizeof(int));
  options->dparam = (double *)calloc(options->dSize, sizeof(double));
  options->dWork = NULL;
  solver_options_nullify(options);

  options->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
  options->dparam[SICONOS_DPARAM_TOL] = 1e-10;

  options->dparam[SICONOS_IPARAM_STOPPING_CRITERION] = SICONOS_STOPPING_CRITERION_RESIDU;
  
}
