
/*!\file solve_lcp.c


   This subroutine allows the resolution of LCP (Linear Complementary Problem).
   Try \f$(z,w)\f$ such that:

\f$
\left\lbrace
\begin{array}{l}
M z- w=q\\
0 \le z \perp w \ge 0\\
\end{array}
\right.
\f$

  here M is an n by n  matrix, q an n-dimensional vector, w an n-dimensional  vector and z an n-dimensional vector.
  This system of equalities and inequalities is solved thanks to @ref lcp solvers.
  The routine's call is due to the function solve_lcp.c.
*/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "SiconosNumerics.h"


/*!\fn int solve_lcp (double *vec, double *q, int * nn, methode *pt,double z[], double w[])

   solve_lcp is a generic interface allowing the call of one of the LCP solvers.

   \param double* : vec On enter a pointer over doubles containing the components of the double matrix with a f90 allocation.
   \param double* : q On enter a pointer over doubles containing the components of the double vector.
   \param int* : nn On enter a pointer over integers, the dimension of the second member.
   \param methde* : pt On enter a pointer other a structure (::methode).
   \param double[] : z On return double vector, the solution of the problem.
   \param double[] : w On return double vector, the solution of the problem.

  \return On return integer, the termination reason (0 is successful otherwise 1).

   \author Nineb Sheherazade.
 */
int solve_lcp(double *vec, double *q, int *nn, methode *pt, double z[], double w[])
{
  int info = -1, it_end, fail;
  const char mot1[10] = "Lemke", mot2[10] = "Gsnl", mot3[10] = "Gcp", mot4[10] = "Latin";
  double res;
  int n = *nn;

  if (strcmp(pt->lcp.nom_method, mot1) == 0)
  {
    //    lemke_lcp_(vec,q,&n,& pt->lcp.itermax,z,w,&it_end,&res,&info);
  }
  else if (strcmp(pt->lcp.nom_method, mot2) == 0)
    gsnl_lcp(vec, q, &n, & pt->lcp.itermax, & pt->lcp.tol, z, w, &it_end, &res, &info);
  else if (strcmp(pt->lcp.nom_method, mot3) == 0)
    gcp_lcp(vec, q, &n, & pt->lcp.itermax, & pt->lcp.tol, z, w, &it_end, &res, &info);
  else if (strcmp(pt->lcp.nom_method, mot4) == 0)
    latin_lcp(vec, q, &n, & pt->lcp.k_latin, & pt->lcp.itermax, & pt->lcp.tol, z, w, &it_end, &res, &info);
  else printf("Warning : Unknown solving method : %s\n", pt->lcp.nom_method);

  return info;
}
