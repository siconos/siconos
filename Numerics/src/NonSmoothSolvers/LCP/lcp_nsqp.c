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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "LCP_Solvers.h"
#include "QP_Solvers.h"

void lcp_nsqp(LinearComplementarity_Problem* problem, double *z, double *w, int *info , Solver_Options* options)
{
  /* matrix M/vector q of the lcp */
  double * M = problem->M->matrix0;
  double * q = problem->q;
  /* size of the LCP */
  int n = problem->size;

  int i, j;

  int nmax;
  int m, me, mmax, mnn;

  double *Q, *A;
  double *p, *b, *xl, *xu;

  double *lambda;

  int lwar, liwar, iout, un;
  int *iwar;
  double *war;

  double tol = options->dparam[0];

  /* / m :        total number of constraints.*/
  m = n;
  /* / me :       number of equality constraints.*/
  me = 0;
  /* /  mmax :     row dimension of a. mmax must be at least one and greater than m.*/
  mmax = m + 1;
  /* /n :        number of variables.
  //nmax :     row dimension of C. nmax must be greater or equal to n.*/
  nmax = n;
  /* /mnn :      must be equal to m + n + n.  */
  mnn = m + n + n;

  for (i = 0; i < n; i++)
  {
    z[i] = 0.0;
    w[i] = 0.0;
  }



  /* / Creation of objective function matrix Q and the the constant vector of the objective function p*/

  Q = (double *)malloc(nmax * nmax * sizeof(double));
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++) Q[j * n + i] = (M[j * n + i] + M[i * n + j]);
  }
  /* /for (i=0;i<n*n;i++) printf("Q[%i] = %g\n",i,Q[i]);*/

  p = (double *)malloc(nmax * sizeof(double));
  for (i = 0; i < n; i++)
    p[i] = q[i] ;
  /* /for (i=0;i<n;i++) printf("p[%i] = %g\n",i,p[i]);*/

  /* / Creation of the data matrix of the linear constraints, A and  the constant data of the linear constraints b*/
  A = (double *)malloc(mmax * nmax * sizeof(double));
  for (i = 0; i < m; i++)
  {
    for (j = 0; j < n; j++) A[j * mmax + i] = M[j * n + i];
  }

  /* /for (i=0;i<mmax*mmax;i++) printf("A[%i] = %g\n",i,A[i]);*/

  b = (double *)malloc(mmax * sizeof(double));
  for (i = 0; i < m; i++) b[i] = q[i] ;

  /* /for (i=0;i<m;i++) printf("b[%i] = %g\n",i,b[i]);*/

  /* / Creation of the the lower and upper bounds for the variables.*/
  xu = (double *)malloc(n * sizeof(double));
  for (i = 0; i < n; i++) xu[i] = 1e300 ;
  xl = (double *)malloc(n * sizeof(double));
  for (i = 0; i < n; i++) xl[i] = 0.0 ;

  /* / on return, lambda contains the lagrange multipliers.*/
  lambda = (double *)malloc(mnn * sizeof(double));
  for (i = 0; i < mnn; i++) lambda[i] = 0.0 ;

  /* /   integer indicating the desired output unit number,*/
  iout = 6;

  /* /   output control.*/
  un = 1;

  /* / real working array. */
  lwar = 3 * nmax * nmax / 2 + 10 * nmax + 2 * mmax;
  war = (double *)malloc(lwar * sizeof(double));
  /* / integer working array. */
  liwar = n ;
  iwar = (int *)malloc(liwar * sizeof(int));
  iwar[0] = 1;


  /* / call ql0001_*/
  /*   F77NAME(ql0001)(m, me, mmax, n, nmax, mnn, Q, p, A, b, xl, xu, */
  /*    z, lambda, iout, *info , un, war, lwar, iwar, liwar, tol); */
  F77NAME(ql0001)(&m, &me, &mmax, &n, &nmax, &mnn, Q, p, A, b, xl, xu,
                  z, lambda, &iout, info, &un, war, &lwar, iwar, &liwar, &tol);

  /* /    printf("tol = %10.4e\n",*tol);
  // for (i=0;i<mnn;i++)printf("lambda[%i] = %g\n",i,lambda[i]);
  // for (i=0;i<n;i++)printf("z[%i] = %g\n",i,z[i]);

  // getting the multiplier due to the lower bounds*/
  for (i = 0; i < n; i++) w[i] = lambda[m + i] ;

  /* / memory freeing*/

  free(Q);
  free(p);
  free(A);
  free(b);
  free(xu);
  free(xl);
  free(lambda);
  free(iwar);
  free(war);
}
int linearComplementarity_nsqp_setDefaultSolverOptions(Solver_Options** arrayOfSolver_Options)
{
  int i;
  if (verbose > 0)
  {
    printf("Set the Default Solver_Options for the NSQP Solver\n");
  }
  int nbSolvers = 1 ;
  Solver_Options * options = (Solver_Options *)malloc(nbSolvers * sizeof(Solver_Options));
  arrayOfSolver_Options[0] = options;


  strcpy(options->solverName, "NSQP");

  options->numberOfInternalSolvers = 1;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 5;
  options->dSize = 5;
  options->iparam = (int *)malloc(options->iSize * sizeof(int));
  options->dparam = (double *)malloc(options->dSize * sizeof(double));
  options->dWork = NULL;
  options->iWork = NULL;
  for (i = 0; i < 5; i++)
  {
    options->iparam[i] = 0;
    options->dparam[i] = 0.0;
  }
  options->dparam[0] = 1e-6;


  return 0;
}

int linearComplementarity_nsqp_deleteDefaultSolverOptions(Solver_Options** arrayOfSolver_Options)
{

  int i;
  if (verbose > 0)
  {
    printf("Set the Default Solver_Options for the NSQP Solver\n");
  }

  Solver_Options * options = arrayOfSolver_Options[0];

  int nbSolvers = 1 ;
  for (i = 0; i < nbSolvers; i++)
  {
    if (options[i].iparam) free(options[i].iparam);
    options[i].iparam = NULL;
    if (options[i].dparam) free(options[i].dparam);
    options[i].dparam = NULL;
    if (options[i].dWork)  free(options[i].dWork);
    options[i].dWork = NULL;
    if (options[i].iWork)  free(options[i].iWork);
    options[i].iWork = NULL;
  }
  free(options);


  return 0;
}
