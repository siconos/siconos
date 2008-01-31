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
int lcp_solver_pred(double *vec, double *q , int *n , method *pt , double *z , double *w ,
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

  if (strcmp(pt->lcp.name , lcpkey1) == 0)
  {

    iparamLCP[0] = pt->lcp.itermax;
    iparamLCP[1] = pt->lcp.chat;

    lcp_lexicolemke(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->lcp.iter = iparamLCP[2];
  }
  /* **** Latin Solver **** */

  else if (strcmp(pt->lcp.name , lcpkey4) == 0)
  {

    iparamLCP[0] = pt->lcp.itermax;
    iparamLCP[1] = pt->lcp.chat;
    dparamLCP[0] = pt->lcp.tol;
    dparamLCP[1] = pt->lcp.k_latin;

    lcp_latin(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->lcp.iter = iparamLCP[2];
    pt->lcp.err  = dparamLCP[2];

  }

  /* **** Latin_w Solver **** */

  else if (strcmp(pt->lcp.name , lcpkey9) == 0)
  {

    iparamLCP[0] = pt->lcp.itermax;
    iparamLCP[1] = pt->lcp.chat;
    dparamLCP[0] = pt->lcp.tol;
    dparamLCP[1] = pt->lcp.k_latin;
    dparamLCP[3] = pt->lcp.relax;


    lcp_latin_w(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->lcp.iter = iparamLCP[2];
    pt->lcp.err  = dparamLCP[2];

  }

  /* **** PGS Solver **** */

  else if (strcmp(pt->lcp.name , lcpkey2) == 0)
  {

    iparamLCP[0] = pt->lcp.itermax;
    iparamLCP[1] = pt->lcp.chat;
    dparamLCP[0] = pt->lcp.tol;
    /* dparamLCP[1] = pt->lcp.relax;*/

    lcp_pgs(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->lcp.iter = iparamLCP[2];
    pt->lcp.err  = dparamLCP[2];

  }
  /* **** NLGS Solver **** */

  else if (strcmp(pt->lcp.name , lcpkey12) == 0)
  {

    printf("Warning: NLGS method is obsolete. Use PGS instead.\n");

  }
  /* **** SOR Solver **** */

  else if (strcmp(pt->lcp.name , lcpkey11) == 0)
  {

    iparamLCP[0] = pt->lcp.itermax;
    iparamLCP[1] = pt->lcp.chat;
    dparamLCP[0] = pt->lcp.tol;
    dparamLCP[1] = pt->lcp.relax;

    lcp_psor(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->lcp.iter = iparamLCP[2];
    pt->lcp.err  = dparamLCP[2];

  }
  /* **** CPG Solver **** */

  else if (strcmp(pt->lcp.name , lcpkey3) == 0)
  {

    iparamLCP[0] = pt->lcp.itermax;
    iparamLCP[1] = pt->lcp.chat;
    dparamLCP[0] = pt->lcp.tol;

    lcp_cpg(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->lcp.iter = iparamLCP[2];
    pt->lcp.err  = dparamLCP[1];

  }

  /* ***** QP Solver ***** */

  else if (strcmp(pt->lcp.name , lcpkey5) == 0)
  {

    /* We assume that the LCP matrix M is symmetric*/

    dparamLCP[0] = pt->lcp.tol;

    lcp_qp(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

  }

  /* **** NSQP Solver **** */

  else if (strcmp(pt->lcp.name , lcpkey6) == 0)
  {

    /* We assume that the LCP matrix M is not symmetric*/

    dparamLCP[0] = pt->lcp.tol;

    lcp_nsqp(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

  }
  else if (strcmp(pt->lcp.name , lcpkey7) == 0)
  {

    iparamLCP[0] = pt->lcp.itermax;
    iparamLCP[1] = pt->lcp.chat;

    lcp_lexicolemke(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->lcp.iter = iparamLCP[2];

  }
  else if (strcmp(pt->lcp.name , lcpkey8) == 0)
  {

    iparamLCP[0] = pt->lcp.itermax;
    iparamLCP[1] = pt->lcp.chat;
    dparamLCP[0] = pt->lcp.tol;

    lcp_newton_min(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->lcp.iter = iparamLCP[2];
    pt->lcp.err  = dparamLCP[1];

  }
  else if (strcmp(pt->lcp.name , lcpkey10) == 0)
  {

    iparamLCP[0] = pt->lcp.itermax;
    iparamLCP[1] = pt->lcp.chat;
    dparamLCP[0] = pt->lcp.tol;

    lcp_newton_FB(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->lcp.iter = iparamLCP[2];
    pt->lcp.err  = dparamLCP[1];

  }
  /* **** RPGS (Regularized Projected Gauss-Seidel) Solver **** */

  else if (strcmp(pt->lcp.name , lcpkey13) == 0)
  {

    iparamLCP[0] = pt->lcp.itermax;
    iparamLCP[1] = pt->lcp.chat;
    dparamLCP[0] = pt->lcp.tol;
    dparamLCP[1] = pt->lcp.rho;
    /* dparamLCP[2] = pt->lcp.relax;*/

    lcp_rpgs(n , vec , q , z , w , &info , iparamLCP , dparamLCP);

    pt->lcp.iter = iparamLCP[2];
    pt->lcp.err  = dparamLCP[3];

  }

  else printf("Warning : Unknown solver : %s\n", pt->lcp.name);

  /* Checking validity of z found  */
  /*  if (info == 0) info = filter_result_LCP(*n,vec,q,z,pt->lcp.tol,pt->lcp.chat,w);*/
  if ((strcmp(pt->lcp.name , lcpkey2) != 0) && (strcmp(pt->lcp.name , lcpkey13) != 0))
    info = filter_result_LCP(*n, vec, q, z, pt->lcp.tol, pt->lcp.chat, w);

  if (info == 0)
  {
    info = extractLCP(vec , n , z , w , indic , indicop , submatlcp , submatlcpop , ipiv , sizesublcp , sizesublcpop , workspace);
    *soltype = 3;
  }

  return info;

}

int extractLCP(double *vec, int *n , double *z , double *w ,
               int *indic , int *indicop , double *submatlcp , double *submatlcpop ,
               int *ipiv , int *sizesublcp , int *sizesublcpop , double *workspace)
{

  int i, j, k, sizelcp, info;
  /*  double *workspace; */
  /*  double epsdiag = 1e-16;*/

  sizelcp = *n;
  /*  workspace = (double*)malloc(sizelcp * sizeof(double)); */

  /*    printf("recalcul_submat\n");*/
  j = 0;
  k = 0;
  for (i = 0; i < sizelcp; i++)
  {
    if (z[i] > w[i])
      /*        if (z[i] >= epsdiag)*/
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
  *sizesublcp = j;
  *sizesublcpop = k;

  if (*sizesublcp != 0)
  {
    for (j = 0; j < *sizesublcp; j++)
    {
      for (i = 0; i < *sizesublcp; i++)
      {
        submatlcp[(j * (*sizesublcp)) + i] = vec[(indic[j] * sizelcp) + indic[i]];
      }
    }

    DGETRF(*sizesublcp, *sizesublcp, submatlcp, *sizesublcp, ipiv, info);
    if (info != 0)
    {
      printf("LU pb in extractLCP !\n"); /*free(workspace)*/;
      return 1;
    }
    /*        DGETRI(*sizesublcp, submatlcp, *sizesublcp, ipiv, workspace , *sizesublcp , info); */
    DGETRI(*sizesublcp, submatlcp, *sizesublcp, ipiv , info);
    if (info != 0)
    {
      printf("LU pb in extractLCP !\n"); /*free(workspace)*/;
      return 1;
    }

    if (*sizesublcpop != 0)
    {
      for (j = 0; j < *sizesublcp; j++)
      {
        for (i = 0; i < *sizesublcpop; i++)
        {
          submatlcpop[(j * (*sizesublcpop)) + i] = vec[(indic[j] * sizelcp) + indicop[i]];
        }
      }
    }

  }
  /*    free(workspace);*/
  return 0;
}

int predictLCP(double *q , int *n , double *z , double *w , double tol,
               int *indic , int *indicop , double *submatlcp , double *submatlcpop ,
               int *ipiv , int *sizesublcp , int *sizesublcpop , double *subq , double *bufz , double *newz)
{

  int i, sizelcp, info, incx;
  /*  double *subq;
  double *bufz;
  double *newz;*/
  double a1 = 1;
  double b1 = 0;
  double error, normq;
  double zi, wi;

  sizelcp = *n;
  /*  subq = (double*)malloc(sizelcp * sizeof(double));
    bufz = (double*)malloc(sizelcp * sizeof(double));
    newz = (double*)malloc(sizelcp * sizeof(double));*/

  incx = 1;
  DCOPY(*n , z , incx , bufz , incx);   /* Saving z on enter. */
  for (i = 0; i < sizelcp; i++)
  {
    z[i] = 0.;
    w[i] = 0.;
  }

  if (*sizesublcp != 0)
  {
    for (i = 0; i < *sizesublcp; i++)
    {
      subq[i] = -q[indic[i]];
    }

    DGEMV(LA_NOTRANS, *sizesublcp , *sizesublcp , a1 , submatlcp , *sizesublcp , subq , incx , b1 , newz , incx);

    for (i = 0; i < *sizesublcp; i++)
    {
      /*        z[indic[i]] = subq[i];*/
      /*        if (newz[i] > 0) z[indic[i]] = newz[i];*/
      z[indic[i]] = newz[i];
    }
  }

  if (*sizesublcpop != 0)
  {
    for (i = 0; i < *sizesublcpop; i++)
    {
      subq[i] = q[indicop[i]];
    }

    if (*sizesublcp != 0)
    {
      a1 = 1.;
      b1 = 1.;
      DGEMV(LA_NOTRANS, *sizesublcpop , *sizesublcp , a1, submatlcpop, *sizesublcpop, newz, incx, b1, subq, incx);
    }

    for (i = 0; i < *sizesublcpop; i++)
    {
      w[indicop[i]] = subq[i];
    }
  }

  error = 0.;
  for (i = 0 ; i < sizelcp ; i++)
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

  normq = DNRM2(sizelcp, q, incx);
  error = error / normq;

  if (error > tol) info = -1;
  else info = *sizesublcp;

  if (info < 0) DCOPY(*n , bufz , incx, z, incx);

  /*  free(subq);
    free(bufz);
    free(newz);*/

  return info;

}
