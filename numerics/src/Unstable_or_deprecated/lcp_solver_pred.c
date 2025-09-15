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
/*!\file lcp_solver_pred.c

This subroutine allows the resolution of LCP (Linear Complementary Problem).\n
Try \f$(z,w)\f$ such that:\n

\f$
\left\lbrace
\begin{array}{l}
w - M z = q\\
0 \le z \perp  w \ge 0\\
\end{array}
\right.
\f$

M is an (\f$ n \times n\f$ ) matrix, q , w and z n-vector. This system of equalities and
inequalities is solved by a prediction on the non-zero z values and linear system solving or
thanks to  lcp solvers.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifndef MEXFLAG
#include "NonSmoothDrivers.h"
#endif
/* int lcp_solver_pred(double *vec, double *q , int *n , method_lcp *pt , double *z , double *w
 * , */
/*                     int firsttime, int *soltype , int *indic , int *indicop , double
 * *submatlcp , double *submatlcpop , */
/*                     int *ipiv , int *sizesublcp , int *sizesublcpop , */
/*                     double *subq , double *bufz , double *newz , double *workspace) */

int lcp_solver_pred(n, method_lcp *pt, double *z, double *w, int firsttime, int *soltype,
                    int *indic, int *indicop, double *submatlcp, double *submatlcpop,
                    int *ipiv, int *sizesublcp, int *sizesublcpop, double *subq, double *bufz,
                    double *newz, double *workspace) {
  /* subq, bufz, newz, workspace: work vectors*/

  const char lcpkey1[10] = "Lemke", lcpkey2[10] = "PGS", lcpkey3[10] = "CPG";
  const char lcpkey4[10] = "Latin", lcpkey5[10] = "QP", lcpkey6[10] = "NSQP";
  const char lcpkey7[15] = "LexicoLemke", lcpkey8[15] = "NewtonMin";
  const char lcpkey9[15] = "Latin_w", lcpkey10[15] = "NewtonFB", lcpkey11[15] = "PSOR";
  const char lcpkey12[10] = "NLGS";
  const char lcpkey13[10] = "RPGS";

  // Remark: Lemke = LexicoLemke. Only one solver is called: lexicoLemke.

  int i, j, info = 1;

  int iparamLCP[5];
  double dparamLCP[5];

  for (i = 0; i < 5; ++i) iparamLCP[i] = 0;
  for (i = 0; i < 5; ++i) dparamLCP[i] = 0.0;

  *soltype = 0;
  /*  limqpos = -DBL_EPSILON / sqrt((double) *n); */
  if (firsttime == 0) {
    i = 0;
    while ((i < (*n - 1)) && (q[i] >= 0.)) i++;
    if ((i == (*n - 1)) && (q[*n - 1] >= 0.)) {
      /* TRIVIAL CASE : q >= 0
       * z = 0 and w = q is solution of LCP(q,M)
       */
      for (j = 0; j < *n; j++) {
        z[j] = 0.0;
        w[j] = q[j];
      }
      *soltype = 1;
      pt->lcp.iter = 0;
      pt->lcp.err = 0.;
      if (pt->lcp.chat > 0) printf("Trivial case of LCP : positive vector q \n");
      return 0;
    }

    info = predictLCP(q, n, z, w, pt->lcp.tol, indic, indicop, submatlcp, submatlcpop, ipiv,
                      sizesublcp, sizesublcpop, subq, bufz, newz);

    if (info >= 0) {
      *soltype = 2;
      pt->lcp.iter = 1;
      if (pt->lcp.chat > 0) printf("LCP solved by prediction on z,w signs\n");
      return 0;
    } else
      info = 1;
  }

  /* Solver name */
  char *name = options->solverName;

  if (verbose == 1)
    printf(
        " ========================== Call %s solver for Linear Complementarity problem "
        "==========================\n",
        name);

  /****** Lemke algorithm ******/
  /* IN: itermax
     OUT: iter */
  if (strcmp(name, "Lemke") == 0 || strcmp(name, "LexicoLemke") == 0)
    lcp_lexicolemke(problem, z, w, &info, options);

  /****** PGS Solver ******/
  /* IN: itermax, tolerance
     OUT: iter, error */
  else if (strcmp(name, "PGS") == 0)
    lcp_pgs(problem, z, w, &info, options);

  /****** CPG Solver ******/
  /* IN: itermax, tolerance
     OUT: iter, error */
  else if (strcmp(name, "CPG") == 0)
    lcp_cpg(problem, z, w, &info, options);

  /****** Latin Solver ******/
  /* IN: itermax, tolerance, k_latin
     OUT: iter, error */
  else if (strcmp(name, "Latin") == 0)
    lcp_latin(problem, z, w, &info, options);

  /****** Latin_w Solver ******/
  /* IN: itermax, tolerance, k_latin, relax
     OUT: iter, error */
  else if (strcmp(name, "Latin_w") == 0)
    lcp_latin_w(problem, z, w, &info, options);

  /****** QP Solver ******/
  /* IN: tolerance
     OUT:
     We assume that the LCP matrix M is symmetric
  */
  else if (strcmp(name, "QP") == 0)
    lcp_qp(problem, z, w, &info, options);

  /****** NSQP Solver ******/
  /* IN: tolerance
     OUT:
  */
  else if (strcmp(name, "NSQP") == 0)
    lcp_nsqp(problem, z, w, &info, options);

  /****** Newton min ******/
  /* IN: itermax, tolerance
     OUT: iter, error
  */
  else if (strcmp(name, "NewtonMin") == 0)
    lcp_newton_min(problem, z, w, &info, options);

  /****** Newton Fischer-Burmeister ******/
  /* IN: itermax, tolerance
     OUT: iter, error
  */
  else if (strcmp(name, "Newton_FB") == 0)
    lcp_newton_FB(problem, z, w, &info, options);

  /****** PSOR Solver ******/
  /* IN: itermax, tolerance, relax
     OUT: iter, error
  */
  else if (strcmp(name, "PSOR") == 0)
    lcp_psor(problem, z, w, &info, options);

  /****** RPGS (Regularized Projected Gauss-Seidel) Solver ******/
  /* IN: itermax, tolerance, rho
     OUT: iter, error
  */
  else if (strcmp(name, "RPGS") == 0)
    lcp_rpgs(problem, z, w, &info, options);

  /****** PATH (Ferris) Solver ******/
  /* IN: itermax, tolerance, rho
     OUT: iter, error
  */
  else if (strcmp(name, "Path") == 0)
    lcp_path(problem, z, w, &info, options);

  else
    printf("LCP_driver error: unknown solver named: %s\n", pt->lcp.name);

  /*************************************************
   *  3 - Check solution validity
   *************************************************/

  /* Warning: it depends on the chosen solver */

  /* Not done for:  PGS, RPGS */
  if ((strcmp(name, "PGS") != 0) && (strcmp(name, "RPGS") != 0))
    info = filter_result_LCP(problem, z, w, options->dparam[0]);

  if (info == 0) {
    info = extractLCP(problem->M, z, indic, indicop, submatlcp, submatlcpop, ipiv, sizesublcp,
                      sizesublcpop);
    *soltype = 3;
  }

  return info;
}

int extractLCP(NumericsMatrix *MGlobal, double *z, int *indic, int *indicop, double *submatlcp,
               double *submatlcpop, int *ipiv, int *sizesublcp, int *sizesublcpop) {
  if (MGlobal == NULL || z == NULL)
    numerics_error("extractLCP", "Null input for one arg (problem, z, ...)");

  int info;
  /*  double epsdiag = DBL_EPSILON;*/

  /* Extract data from problem */
  if (MGlobal->storageType == 1)
    numerics_error("extractLCP", "Not yet implemented for sparse storage");
  double *M = MGlobal->matrix0;
  int sizelcp = MGlobal->size0;
  if (M == NULL) numerics_error("extractLCP", "Null input matrix M");

  /*  workspace = (double*)malloc(sizelcp * sizeof(double)); */
  /*    printf("recalcul_submat\n");*/

  /* indic = set of indices for which z[i] is positive */
  /* indicop = set of indices for which z[i] is null */

  /* test z[i] sign */
  int i, j = 0, k = 0;
  for (i = 0; i < sizelcp; i++) {
    if (z[i] > w[i]) /* if (z[i] >= epsdiag)*/
    {
      indic[j] = i;
      j++;
    } else {
      indicop[k] = i;
      k++;
    }
  }

  /* size of the sub-matrix that corresponds to indic */
  *sizesublcp = j;
  /* size of the sub-matrix that corresponds to indicop */
  *sizesublcpop = k;

  /* If indic is non-empty, copy corresponding M sub-matrix into submatlcp */
  if (*sizesublcp != 0) {
    for (j = 0; j < *sizesublcp; j++) {
      for (i = 0; i < *sizesublcp; i++)
        submatlcp[(j * (*sizesublcp)) + i] = M[(indic[j] * sizelcp) + indic[i]];
    }

    /* LU factorization and inverse in place for submatlcp */
    DGETRF(*sizesublcp, *sizesublcp, submatlcp, *sizesublcp, ipiv, info);
    if (info != 0) {
      numerics_warning("extractLCP", "LU factorization failed");
      return 1;
    }

    DGETRI(*sizesublcp, submatlcp, *sizesublcp, ipiv, info);
    if (info != 0) {
      numerics_warning("extractLCP", "LU inversion failed");
      return 1;
    }

    /* if indicop is not empty, copy corresponding M sub-matrix into submatlcpop */
    if (*sizesublcpop != 0) {
      for (j = 0; j < *sizesublcp; j++) {
        for (i = 0; i < *sizesublcpop; i++)
          submatlcpop[(j * (*sizesublcpop)) + i] = vec[(indic[j] * sizelcp) + indicop[i]];
      }
    }
  }

  return 0;
}

int predictLCP(int sizeLCP, double *q, double *z, double *w, double tol, int *indic,
               int *indicop, double *submatlcp, double *submatlcpop, int *ipiv,
               int *sizesublcp, int *sizesublcpop, double *subq, double *bufz, double *newz) {
  if (q == NULL || z == NULL || w == NULL)
    numerics_error("predictLCP", "Null input for one arg (problem, q,w ...)");

  int i, sizelcp, info, incx;
  double error, norm_q;
  double zi, wi;
  int incx = 1;

  /* Copy of z into a buffer for restart if predict failed */
  cblas_dcopy(sizeLCP, z, incx, bufz, incx); /* Saving z on enter. */

  /* Sets z and w to 0*/
  for (i = 0; i < sizeLCP; i++) {
    z[i] = 0.;
    w[i] = 0.;
  }

  /* if indic is not empty, computes solution of newz of submatlcp.newz = subq */
  if (*sizesublcp != 0) {
    /* Gets subq */
    for (i = 0; i < *sizesublcp; i++) subq[i] = -q[indic[i]];

    cblas_dgemv(CblasColMajor, CblasNoTrans, *sizesublcp, *sizesublcp, 1.0, submatlcp,
                *sizesublcp, subq, incx, 0.0, newz, incx);

    /* Copy of newz into z for i in indic */
    for (i = 0; i < *sizesublcp; i++) {
      /*        z[indic[i]] = subq[i];*/
      /*        if (newz[i] > 0) z[indic[i]] = newz[i];*/
      z[indic[i]] = newz[i];
    }
  }

  /* if indicop is not empty, computes subw = submatlcpop.newz + subq - subw saved in subq */
  if (*sizesublcpop != 0) {
    /* Gets subq */
    for (i = 0; i < *sizesublcpop; i++) subq[i] = q[indicop[i]];

    if (*sizesublcp != 0)
      cblas_dgemv(CblasColMajor, CblasNoTrans, *sizesublcpop, *sizesublcp, 1.0, submatlcpop,
                  *sizesublcpop, newz, incx, 1.0, subq, incx);

    /* Copy of subq=subw into w for indices in indicop */
    for (i = 0; i < *sizesublcpop; i++) w[indicop[i]] = subq[i];
  }

  /* Error evaluation */
  error = 0.;
  for (i = 0; i < sizeLCP; i++) {
    zi = z[i];
    wi = w[i];
    if (zi < 0.0) {
      error += -zi;
      if (wi < 0.0) error += zi * wi;
    }
    if (wi < 0.0) error += -wi;
    if ((zi > 0.0) && (wi > 0.0)) error += zi * wi;
  }

  norm_q = cblas_dnrm2(sizeLCP, q, incx);
  error = error / norm_q;

  if (error > tol) {
    printf(
        "Numerics warning - predictLCP failed, error = %g > tolerance = %g - Reset z to "
        "starting value.\n",
        error, tol);
    info = -1;
  } else
    info = *sizesublcp;

  /* If failed, reset z to starting value (saved in bufz) */
  if (info < 0) cblas_dcopy(sizeLCP, bufz, incx, z, incx);

  return info;
}
