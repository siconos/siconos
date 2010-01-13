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
#include "NonSmoothNewton.h"
#include "Numerics_Options.h"
#include "LA.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"

/* Linesearch */
void linesearch_Armijo(int n, double *z, double* dir, double psi_k, double descentCondition, NewtonFunctionPtr* phi)
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
  DAXPY(n , 1.0 , dir , incx , z , incy);

  while (tk > tmin)
  {
    /* Computes merit function = 1/2*norm(phi(z_{k+1}))^2 */
    (*phi)(n, z, phiVector, 0);
    merit =  DNRM2(n, phiVector , incx);
    merit = 0.5 * merit * merit;
    merit_k = psi_k + sigma * tk * descentCondition;
    if (merit < merit_k) break;
    tk = tk * 0.5;
    /* Computes z_k+1 = z0 + tk.dir
    warning: (-tk) because we need to start from z0 at each step while z_k+1 is saved in place of z_k ...*/
    DAXPY(n , -tk , dir , incx , z , incy);
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
  int infoDGESV;

  /* Memory allocation for phi and its jacobian */
  double * phiVector = (double*)malloc(n * sizeof(*phiVector));
  double *jacobianPhiMatrix = (double*)malloc(n2 * sizeof(*jacobianPhiMatrix));
  /** merit function and its jacobian */
  double psi;
  double *jacobian_psi = (double*)malloc(n * sizeof(*jacobian_psi));
  int* ipiv = (int *)malloc(n * sizeof(*ipiv));
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
    DGEMV(LA_TRANS, n, n, 1.0, jacobianPhiMatrix, n, phiVector, incx, 0.0, jacobian_psi, incx);
    norm_jacobian_psi = DNRM2(n, jacobian_psi, 1);

    /* Computes norm2(phi) */
    normPhi = DNRM2(n, phiVector, 1);
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
    DSCAL(n , -1.0 , phiVector, incx);
    DGESV(n, 1, jacobianPhiMatrix, n, ipiv, phiVector, n, infoDGESV);

    /* descentCondition = jacobian_psi.dk */
    descentCondition = DDOT(n, jacobian_psi,  1,  phiVector, 1);

    /* Criterion to be satisfied: error < -rho*norm(dk)^p */
    criterion = DNRM2(n, phiVector, 1);
    criterion = -rho * pow(criterion, p);

    if (infoDGESV != 0 || descentCondition > criterion)
    {
      /* dk = - jacobian_psi (remind that dk is saved in phiVector) */
      DCOPY(n, jacobian_psi, 1, phiVector, 1);
      DSCAL(n , -1.0 , phiVector, incx);
    }

    /* Step-3 Line search: computes z_k+1 */
    linesearch_Armijo(n, z, phiVector, psi, descentCondition, phi);

    if (verbose > 0)
    {
      printf("Non Smooth Newton, iteration number %i, error equal to %14.7e .\n", niter, terminationCriterion);
      printf(" -----------------------------------------------------------------------\n");
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
  int infoDGESV;

  /* Memory allocation for phi and its jacobian */
  double * phiVector = (double*)malloc(n * sizeof(*phiVector));
  double *jacobianPhiMatrix = (double*)malloc(n2 * sizeof(*jacobianPhiMatrix));
  /** merit function and its jacobian */
  double psi;
  double *jacobian_psi = (double*)malloc(n * sizeof(*jacobian_psi));
  int* ipiv = (int *)malloc(n * sizeof(*ipiv));
  if (phiVector == NULL || jacobianPhiMatrix == NULL ||  jacobian_psi == NULL || ipiv == NULL)
  {
    fprintf(stderr, "NonSmoothNewton, memory allocation failed.\n");
    exit(EXIT_FAILURE);
  }

  /** The algorithm is alg 4.1 of the paper of Kanzow and Kleinmichel, "A new class of semismooth Newton-type methods
      for nonlinear complementarity problems", in Computational Optimization and Applications, 11, 227-251 (1998).

      We try to keep the same notations
  */

  double norm_jacobian_psi, normPhi;
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
    DGEMV(LA_TRANS, n, n, 1.0, jacobianPhiMatrix, n, phiVector, incx, 0.0, jacobian_psi, incx);
    norm_jacobian_psi = DNRM2(n, jacobian_psi, 1);

    /* Computes norm2(phi) */
    normPhi = DNRM2(n, phiVector, 1);
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
    DSCAL(n , -1.0 , phiVector, incx);
    DGESV(n, 1, jacobianPhiMatrix, n, ipiv, phiVector, n, infoDGESV);

    double tk = -1;

    DAXPY(n , tk , phiVector , 1 , z , 1);



    if (verbose > 0)
    {
      printf("Non Smooth Newton, iteration number %i, error equal to %14.7e .\n", niter, terminationCriterion);
      printf(" -----------------------------------------------------------------------\n");
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
