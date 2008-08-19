/* Siconos-Numerics version 2.0.1, Copyright INRIA 2005-2006.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
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

M is an (\f$ n \times n\f$ ) matrix, q , w and z n-vector. This system of equalities and inequalities
is solved by a prediction on the non-zero z values and linear system solving or thanks to @ref lcp solvers.
*/

/**
   lcp_solver_pred is a generic interface allowing the call of one of the LCP solvers.
   - At first, the signs of q elements are checked to detect the trivial case of positive q.\n
   - Then it tries to find z by assuming that the indices of the non-zero elements
   are the same as the previous solution (previous solutions with trivial case of q positive excepted of course).\n
   - If q is not positive and prediction failed, the regular LCP solver is called.\n

   \param[in] vec          On enter, a (\f$n \times n\f$)-vector of doubles which contains the components of the LCP matrix with a Fortran storage.
   \param[in] q            On enter, a n-vector of doubles which contains the components of the constant right hand side vector.
   \param[in] n            On enter, an integer which represents the dimension of the LCP problem.
   \param[in] pt           On enter, a union containing the LCP structure.
   \n \n
   \param[in,out] z        On enter, an initial guess for iterative LCP solvers.\n
   On return, a n-vector of doubles which contains the solution of the problem.
   \param[out] w           On return, a n-vector of doubles which contains the complementary solution of the problem.
   \n
   \param[in]     firsttime      At 1, forces the regular LCP solver to be used (for initialization purpose).
   \param[out]    soltype        On return, indicates how the solution was found (0 : no sol,1 : q positive,2 : prediction,3 : LCP solver)
   \param[in,out] indic          The set of indices of non-zero z values in ascending order
   \param[in,out] indicop        The complementary set of indices of "indic".
   \param[in,out] submatlcp      The submatrix of M defined by "indic".
   \param[in,out] submatlcpop    The submatrix of M defined by "indicop".
   \param[in,out] ipiv           Pivot indices in LU factorization of "submatlcp".
   \param[in,out] sizesublcp     "submatlcp" size.
   \param[in,out] sizesublcpop   "submatlcpop" size.

   \return integer
   - 0 : successful\n
   - >0 : otherwise (see specific solvers for more information about the log info)

   \author Nineb Sheherazade & Mathieu Renouf & Pascal Denoyelle
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "LA.h"

#ifndef MEXFLAG
#include "NonSmoothDrivers.h"
#endif
int lcp_solver_pred(double *vec, double *q , int *n , method_lcp *pt , double *z , double *w ,
                    int firsttime, int *soltype , int *indic , int *indicop , double *submatlcp , double *submatlcpop ,
                    int *ipiv , int *sizesublcp , int *sizesublcpop ,
                    double *subq , double *bufz , double *newz , double *workspace)
int lcp_solver_pred(n , method_lcp *pt , double *z , double *w ,
                    int firsttime, int *soltype , int *indic , int *indicop , double *submatlcp , double *submatlcpop ,
                    int *ipiv , int *sizesublcp , int *sizesublcpop ,
                    double *subq , double *bufz , double *newz , double *workspace)
{

  /* subq, bufz, newz, workspace: work vectors*/

  const char lcpkey1[10] = "Lemke", lcpkey2[10] = "PGS", lcpkey3[10] = "CPG";
  const char lcpkey4[10] = "Latin", lcpkey5[10] = "QP", lcpkey6[10] = "NSQP";
  const char lcpkey7[15] = "LexicoLemke", lcpkey8[15] = "NewtonMin";
  const char lcpkey9[15] = "Latin_w", lcpkey10[15] = "NewtonFB", lcpkey11[15] = "PSOR";
  const char lcpkey12[10] = "NLGS";
  const char lcpkey13[10] = "RPGS";

  // Remark: Lemke = LexicoLemke. Only one solver is called: lexicoLemke.

  int i, j, info = 1;

  int     iparamLCP[5];
  double  dparamLCP[5];

  for (i = 0 ; i < 5 ; ++i) iparamLCP[i] = 0;
  for (i = 0 ; i < 5 ; ++i) dparamLCP[i] = 0.0;

  *soltype = 0;
  /*  limqpos = -1e-16 / sqrt((double) *n); */
  if (firsttime == 0)
  {
    i = 0;
    while ((i < (*n - 1)) && (q[i] >= 0.)) i++;
    if ((i == (*n - 1)) && (q[*n - 1] >= 0.))
    {
      /* TRIVIAL CASE : q >= 0
       * z = 0 and w = q is solution of LCP(q,M)
       */
      for (j = 0 ; j < *n; j++)
      {
        z[j] = 0.0;
        w[j] = q[j];
      }
      *soltype = 1;
      pt->lcp.iter = 0;
      pt->lcp.err  = 0.;
      if (pt->lcp.chat > 0) printf("Trivial case of LCP : positive vector q \n");
      return 0;
    }

    info = predictLCP(q , n , z , w , pt->lcp.tol,
                      indic , indicop , submatlcp , submatlcpop ,
                      ipiv , sizesublcp , sizesublcpop , subq , bufz , newz);

    if (info >= 0)
    {
      *soltype = 2;
      pt->lcp.iter = 1;
      if (pt->lcp.chat > 0) printf("LCP solved by prediction on z,w signs\n");
      return 0;
    }
    else info = 1;
  }

  /* Solver name */
  char * name = options->solverName;

  if (verbose == 1)
    printf(" ========================== Call %s solver for Linear Complementarity problem ==========================\n", name);

  /****** Lemke algorithm ******/
  /* IN: itermax
     OUT: iter */
  if (strcmp(name, "Lemke") == 0 || strcmp(name, "LexicoLemke") == 0)
    lcp_lexicolemke(problem, z , w , &info , options);

  /****** PGS Solver ******/
  /* IN: itermax, tolerance
     OUT: iter, error */
  else if (strcmp(name, "PGS") == 0)
    lcp_pgs(problem, z , w , &info , options);

  /****** CPG Solver ******/
  /* IN: itermax, tolerance
     OUT: iter, error */
  else if (strcmp(name, "CPG") == 0)
    lcp_cpg(problem, z , w , &info , options);

  /****** Latin Solver ******/
  /* IN: itermax, tolerance, k_latin
     OUT: iter, error */
  else if (strcmp(name, "Latin") == 0)
    lcp_latin(problem, z , w , &info , options);

  /****** Latin_w Solver ******/
  /* IN: itermax, tolerance, k_latin, relax
     OUT: iter, error */
  else if (strcmp(name, "Latin_w") == 0)
    lcp_latin_w(problem, z , w , &info , options);

  /****** QP Solver ******/
  /* IN: tolerance
     OUT:
     We assume that the LCP matrix M is symmetric
  */
  else if (strcmp(name, "QP") == 0)
    lcp_qp(problem, z , w , &info , options);

  /****** NSQP Solver ******/
  /* IN: tolerance
     OUT:
  */
  else if (strcmp(name, "NSQP") == 0)
    lcp_nsqp(problem, z , w , &info , options);

  /****** Newton min ******/
  /* IN: itermax, tolerance
     OUT: iter, error
  */
  else if (strcmp(name, "NewtonMin") == 0)
    lcp_newton_min(problem, z , w , &info , options);

  /****** Newton Fischer-Burmeister ******/
  /* IN: itermax, tolerance
     OUT: iter, error
  */
  else if (strcmp(name, "Newton_FB") == 0)
    lcp_newton_FB(problem, z , w , &info , options);

  /****** PSOR Solver ******/
  /* IN: itermax, tolerance, relax
     OUT: iter, error
  */
  else if (strcmp(name, "PSOR") == 0)
    lcp_psor(problem, z , w , &info , options);

  /****** RPGS (Regularized Projected Gauss-Seidel) Solver ******/
  /* IN: itermax, tolerance, rho
     OUT: iter, error
  */
  else if (strcmp(name, "RPGS") == 0)
    lcp_rpgs(problem, z , w , &info , options);

  /****** PATH (Ferris) Solver ******/
  /* IN: itermax, tolerance, rho
     OUT: iter, error
  */
  else if (strcmp(name, "Path") == 0)
    lcp_path(problem, z , w , &info , options);

  else
    printf("LCP_driver error: unknown solver named: %s\n", pt->lcp.name);

  /*************************************************
   *  3 - Check solution validity
   *************************************************/

  /* Warning: it depends on the chosen solver */

  /* Not done for:  PGS, RPGS */
  if ((strcmp(name, "PGS") != 0) && (strcmp(name, "RPGS") != 0))
    info = filter_result_LCP(problem, z, w, options->dparam[0]);



  if (info == 0)
  {
    info = extractLCP(problem->M , z, indic , indicop , submatlcp , submatlcpop , ipiv , sizesublcp , sizesublcpop);
    *soltype = 3;
  }

  return info;

}

int extractLCP(NumericsMatrix* MGlobal, double *z , int *indic, int *indicop, double *submatlcp , double *submatlcpop,
               int *ipiv , int *sizesublcp , int *sizesublcpop)
{
  if (MGlobal == NULL || z == NULL)
    numericsError("extractLCP", "Null input for one arg (problem, z, ...)");

  int info;
  /*  double epsdiag = 1e-16;*/

  /* Extract data from problem */
  if (MGlobal->storageType == 1)
    numericsError("extractLCP", "Not yet implemented for sparse storage");
  double * M = MGlobal->matrix0;
  int sizelcp = MGlobal->size0;
  if (M == NULL)
    numericsError("extractLCP", "Null input matrix M");

  /*  workspace = (double*)malloc(sizelcp * sizeof(double)); */
  /*    printf("recalcul_submat\n");*/


  /* indic = set of indices for which z[i] is positive */
  /* indicop = set of indices for which z[i] is null */

  /* test z[i] sign */
  int i, j = 0, k = 0;
  for (i = 0; i < sizelcp; i++)
  {
    if (z[i] > w[i]) /* if (z[i] >= epsdiag)*/
    {
      indic[j] = i;
      j++;
    }
    else
    {
      indicop[k] = i;
      k++;
    }
  }

  /* size of the sub-matrix that corresponds to indic */
  *sizesublcp = j;
  /* size of the sub-matrix that corresponds to indicop */
  *sizesublcpop = k;

  /* If indic is non-empty, copy corresponding M sub-matrix into submatlcp */
  if (*sizesublcp != 0)
  {
    for (j = 0; j < *sizesublcp; j++)
    {
      for (i = 0; i < *sizesublcp; i++)
        submatlcp[(j * (*sizesublcp)) + i] = M[(indic[j] * sizelcp) + indic[i]];
    }

    /* LU factorization and inverse in place for submatlcp */
    DGETRF(*sizesublcp, *sizesublcp, submatlcp, *sizesublcp, ipiv, info);
    if (info != 0)
    {
      numericsWarning("extractLCP", "LU factorization failed") ;
      return 1;
    }

    DGETRI(*sizesublcp, submatlcp, *sizesublcp, ipiv , info);
    if (info != 0)
    {
      numericsWarning("extractLCP", "LU inversion failed");
      return 1;
    }

    /* if indicop is not empty, copy corresponding M sub-matrix into submatlcpop */
    if (*sizesublcpop != 0)
    {
      for (j = 0; j < *sizesublcp; j++)
      {
        for (i = 0; i < *sizesublcpop; i++)
          submatlcpop[(j * (*sizesublcpop)) + i] = vec[(indic[j] * sizelcp) + indicop[i]];
      }
    }
  }

  return 0;
}

int predictLCP(int sizeLCP, double* q, double *z , double *w , double tol,
               int *indic , int *indicop , double *submatlcp , double *submatlcpop ,
               int *ipiv , int *sizesublcp , int *sizesublcpop , double *subq , double *bufz , double *newz)
{
  if (q == NULL ||  z == NULL || w == NULL)
    numericsError("predictLCP", "Null input for one arg (problem, q,w ...)");

  int i, sizelcp, info, incx;
  double error, normq;
  double zi, wi;
  int incx = 1;

  /* Copy of z into a buffer for restart if predict failed */
  DCOPY(sizeLCP, z , incx , bufz , incx);   /* Saving z on enter. */

  /* Sets z and w to 0*/
  for (i = 0; i < sizeLCP; i++)
  {
    z[i] = 0.;
    w[i] = 0.;
  }

  /* if indic is not empty, computes solution of newz of submatlcp.newz = subq */
  if (*sizesublcp != 0)
  {
    /* Gets subq */
    for (i = 0; i < *sizesublcp; i++)
      subq[i] = -q[indic[i]];


    DGEMV(LA_NOTRANS, *sizesublcp , *sizesublcp , 1.0 , submatlcp , *sizesublcp , subq , incx , 0.0 , newz , incx);

    /* Copy of newz into z for i in indic */
    for (i = 0; i < *sizesublcp; i++)
    {
      /*        z[indic[i]] = subq[i];*/
      /*        if (newz[i] > 0) z[indic[i]] = newz[i];*/
      z[indic[i]] = newz[i];
    }
  }

  /* if indicop is not empty, computes subw = submatlcpop.newz + subq - subw saved in subq */
  if (*sizesublcpop != 0)
  {
    /* Gets subq */
    for (i = 0; i < *sizesublcpop; i++)
      subq[i] = q[indicop[i]];

    if (*sizesublcp != 0)
      DGEMV(LA_NOTRANS, *sizesublcpop , *sizesublcp , 1.0, submatlcpop, *sizesublcpop, newz, incx, 1.0, subq, incx);

    /* Copy of subq=subw into w for indices in indicop */
    for (i = 0; i < *sizesublcpop; i++)
      w[indicop[i]] = subq[i];
  }

  /* Error evaluation */
  error = 0.;
  for (i = 0 ; i < sizeLCP ; i++)
  {
    zi = z[i];
    wi = w[i];
    if (zi < 0.0)
    {
      error += -zi;
      if (wi < 0.0) error += zi * wi;
    }
    if (wi < 0.0) error += -wi;
    if ((zi > 0.0) && (wi > 0.0)) error += zi * wi;
  }

  normq = DNRM2(sizeLCP, q, incx);
  error = error / normq;

  if (error > tol)
  {
    printf("Numerics warning - predictLCP failed, error = %g > tolerance = %g - Reset z to starting value.\n", error, tol);
    info = -1;
  }
  else info = *sizesublcp;

  /* If failed, reset z to starting value (saved in bufz) */
  if (info < 0)
    DCOPY(sizeLCP , bufz , incx, z, incx);

  return info;

}
