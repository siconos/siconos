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
#include "numerics_verbose.h"
#include "NonSmoothNewton.h"
#include "SiconosLapack.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"

void linesearch_Armijo(int n, double *z, double* dir, double psi_k,
                       double descentCondition, NewtonFunctionPtr* phi)
{
  double * phiVector = (double*)malloc(n * sizeof(*phiVector));
  if (phiVector == NULL)
  {
    fprintf(stderr, "NonSmoothNewton::linesearch_Armijo, memory allocation failed for phiVector\n");
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
    if (verbose > 0)
      printf("NonSmoothNewton::linesearch_Armijo warning, resulting tk < tmin, linesearch stopped.\n");

}

int nonSmoothNewton(int n, double* z, NewtonFunctionPtr* phi, NewtonFunctionPtr* jacobianPhi, int* iparam, double* dparam)
{
  if (phi == NULL || jacobianPhi == NULL)
  {
    fprintf(stderr, "NonSmoothNewton error: phi or its jacobian function = NULL pointer.\n");
    exit(EXIT_FAILURE);
  }

  int itermax = iparam[0]; // maximum number of iterations allowed
  int niter = 0; // current iteration number
  double tolerance = dparam[0];
  if (verbose > 0)
  {
    printf(" ============= Starting of Newton process =============\n");
    printf(" - tolerance: %14.7e\n - maximum number of iterations: %i\n", tolerance, itermax);
  }

  int incx = 1;
  int n2 = n * n;
  lapack_int infoDGESV;

  /* Memory allocation for phi and its jacobian */
  double * phiVector = (double*)malloc(n * sizeof(*phiVector));
  double *jacobianPhiMatrix = (double*)malloc(n2 * sizeof(*jacobianPhiMatrix));
  /** merit function and its jacobian */
  double psi;
  double *jacobian_psi = (double*)malloc(n * sizeof(*jacobian_psi));
  lapack_int* ipiv = (lapack_int *)malloc(n * sizeof(lapack_int));
  if (phiVector == NULL || jacobianPhiMatrix == NULL ||  jacobian_psi == NULL || ipiv == NULL)
  {
    fprintf(stderr, "NonSmoothNewton, memory allocation failed.\n");
    exit(EXIT_FAILURE);
  }

  /** The algorithm is alg 4.1 of the paper of Kanzow and Kleinmichel, "A new class of semismooth Newton-type methods
      for nonlinear complementarity problems", in Computational Optimization and Applications, 11, 227-251 (1998).

      We try to keep the same notations
  */

  double rho = 1e-8;
  double descentCondition, criterion, norm_jacobian_psi, normPhi;
  double p = 2.1;
  double terminationCriterion = 1;
  if (jacobian_psi == NULL)
  {
    fprintf(stderr, "NonSmoothNewton, memory allocation failed for jacobian_psi.\n");
    exit(EXIT_FAILURE);
  }

  /** Iterations ... */
  while ((niter < itermax) && (terminationCriterion > tolerance))
  {
    ++niter;
    /** Computes phi and its jacobian */
    (*phi)(n, z, phiVector, 0);
    (*jacobianPhi)(n, z, jacobianPhiMatrix, 1);
    /* Computes the jacobian of the merit function, jacobian_psi = transpose(jacobianPhiMatrix).phiVector */
    cblas_dgemv(CblasColMajor,CblasTrans, n, n, 1.0, jacobianPhiMatrix, n, phiVector, incx, 0.0, jacobian_psi, incx);
    norm_jacobian_psi = cblas_dnrm2(n, jacobian_psi, 1);

    /* Computes norm2(phi) */
    normPhi = cblas_dnrm2(n, phiVector, 1);
    /* Computes merit function */
    psi = 0.5 * normPhi * normPhi;

    /* Stops if the termination criterion is satisfied */
    terminationCriterion = norm_jacobian_psi;
    if (terminationCriterion < tolerance)
      break;

    /* Search direction calculation
    Find a solution dk of jacobianPhiMatrix.d = -phiVector.
    dk is saved in phiVector.
    */
    cblas_dscal(n , -1.0 , phiVector, incx);
    DGESV(n, 1, jacobianPhiMatrix, n, ipiv, phiVector, n, &infoDGESV);

    /* descentCondition = jacobian_psi.dk */
    descentCondition = cblas_ddot(n, jacobian_psi,  1,  phiVector, 1);

    /* Criterion to be satisfied: error < -rho*norm(dk)^p */
    criterion = cblas_dnrm2(n, phiVector, 1);
    criterion = -rho * pow(criterion, p);

    if (infoDGESV != 0 || descentCondition > criterion)
    {
      /* dk = - jacobian_psi (remind that dk is saved in phiVector) */
      cblas_dcopy(n, jacobian_psi, 1, phiVector, 1);
      cblas_dscal(n , -1.0 , phiVector, incx);
    }

    /* Step-3 Line search: computes z_k+1 */
    linesearch_Armijo(n, z, phiVector, psi, descentCondition, phi);

    if (verbose > 0)
    {
      printf("Non Smooth Newton, iteration number %i, error equal to %14.7e .\n", niter, terminationCriterion);
      printf(" -----------\n");
    }
  }

  /* Total number of iterations */
  iparam[1] = niter;
  /* Final error */
  dparam[1] = terminationCriterion;

  /** Free memory*/
  free(phiVector);
  free(jacobianPhiMatrix);
  free(jacobian_psi);
  free(ipiv);

  if (verbose > 0)
  {
    if (dparam[1] > tolerance)
      printf("Non Smooth Newton warning: no convergence after %i iterations\n" , niter);

    else
      printf("Non Smooth Newton: convergence after %i iterations\n" , niter);
    printf(" The residue is : %e \n", dparam[1]);
  }

  if (dparam[1] > tolerance)
    return 1;
  else return 0;
}

int nonSmoothDirectNewton(int n, double* z, NewtonFunctionPtr* phi, NewtonFunctionPtr* jacobianPhi, int* iparam, double* dparam)
{
  if (phi == NULL || jacobianPhi == NULL)
  {
    fprintf(stderr, "NonSmoothNewton error: phi or its jacobian function = NULL pointer.\n");
    exit(EXIT_FAILURE);
  }

  int itermax = iparam[0]; // maximum number of iterations allowed
  int niter = 0; // current iteration number
  double tolerance = dparam[0];
  if (verbose > 0)
  {
    printf(" ============= Starting of Newton process =============\n");
    printf(" - tolerance: %14.7e\n - maximum number of iterations: %i\n", tolerance, itermax);
  }

  int incx = 1;
  int n2 = n * n;
  lapack_int infoDGESV = 0;

  /* Memory allocation for phi and its jacobian */
  double * phiVector = (double*)malloc(n * sizeof(*phiVector));
  double *jacobianPhiMatrix = (double*)malloc(n2 * sizeof(*jacobianPhiMatrix));
  /** merit function and its jacobian */
  double *jacobian_psi = (double*)malloc(n * sizeof(*jacobian_psi));
  lapack_int* ipiv = (lapack_int *)malloc(n * sizeof(lapack_int));
  if (phiVector == NULL || jacobianPhiMatrix == NULL ||  jacobian_psi == NULL || ipiv == NULL)
  {
    fprintf(stderr, "NonSmoothNewton, memory allocation failed.\n");
    exit(EXIT_FAILURE);
  }

  /** The algorithm is alg 4.1 of the paper of Kanzow and Kleinmichel, "A new class of semismooth Newton-type methods
      for nonlinear complementarity problems", in Computational Optimization and Applications, 11, 227-251 (1998).

      We try to keep the same notations
  */

  double norm_jacobian_psi;
  double terminationCriterion = 1;
  if (jacobian_psi == NULL)
  {
    fprintf(stderr, "NonSmoothNewton, memory allocation failed for jacobian_psi.\n");
    exit(EXIT_FAILURE);
  }

  /** Iterations ... */
  while ((niter < itermax) && (terminationCriterion > tolerance))
  {
    ++niter;
    /** Computes phi and its jacobian */
    (*phi)(n, z, phiVector, 0);
    (*jacobianPhi)(n, z, jacobianPhiMatrix, 1);
    /* Computes the jacobian of the merit function, jacobian_psi = transpose(jacobianPhiMatrix).phiVector */
    cblas_dgemv(CblasColMajor,CblasTrans, n, n, 1.0, jacobianPhiMatrix, n, phiVector, incx, 0.0, jacobian_psi, incx);
    norm_jacobian_psi = cblas_dnrm2(n, jacobian_psi, 1);

    /* Computes norm2(phi) */
    // normPhi = cblas_dnrm2(n, phiVector, 1);
    /* Computes merit function */
    //psi = 0.5 * normPhi * normPhi;

    /* Stops if the termination criterion is satisfied */
    terminationCriterion = norm_jacobian_psi;
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



    if (verbose > 0)
    {
      printf("Non Smooth Newton, iteration number %i, error equal to %14.7e .\n", niter, terminationCriterion);
      printf(" -----------\n");
    }
  }

  /* Total number of iterations */
  iparam[1] = niter;
  /* Final error */
  dparam[1] = terminationCriterion;

  /** Free memory*/
  free(phiVector);
  free(jacobianPhiMatrix);
  free(jacobian_psi);
  free(ipiv);

  if (verbose > 0)
  {
    if (dparam[1] > tolerance)
      printf("Non Smooth Newton warning: no convergence after %i iterations\n" , niter);

    else
      printf("Non Smooth Newton: convergence after %i iterations\n" , niter);
    printf(" The residue is : %e \n", dparam[1]);
  }

  if (dparam[1] > tolerance)
    return 1;
  else return 0;
}
