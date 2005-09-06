
/*!\file solve_cfp.c

  This subroutine allows the primal resolution of contact problems with frictio.

  Try \f$(z,w)\f$ such that:
\f$
\left\lbrace
\begin{array}{l}
M z- w=q\\
0 \le z_n \perp w_n \ge 0\\
-w_t \in \partial\psi_{[-\mu z_n, \mu z_n]}(z_t)\\
\end{array}
\right.
\f$

 here M is an n by n  matrix, q an n-dimensional vector, z an n-dimensional  vector and w an n-dimensional vector.

 This system of equations and inequalities is solved thanks to @ref pfc solvers.
 The routine's call is due to the function solve_cfp.c.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef MEXFLAG
#include "solverpack.h"
#endif
#include <time.h>

/*!\fn   int solve_cfp (double *vec,double *q,int *nn, methode *pt,double z[],double w[])

   solve_cfp is a generic interface allowing the call of one of the PFC solvers.

   \param double* : vec On enter double vector containing the components of the double matrix with a fortran90 allocation.
   \param double* : q On enter a pointer over doubles containing the components of the second member of the system.
   \param int* : nn On enter a pointer over integers, the dimension of the second member.
   \param methode* : pt On enter a pointer other a structure (::methode).
   \param double[] : z On return real vector, the solution of the problem.
   \param double[] : w On return real vector, the solution of the problem.

   \return On return int, the termination reason (0 is successful otherwise 1).

   \author Nineb Sheherazade.
 */

int solve_cfp(double *vec, double *q, int *nn, methode *pt, double z[], double w[])
{
  int choix, it_end, fail;
  char mot1[10] = "Gsnl", mot2[10] = "Gcp", mot3[10] = "Latin";
  double res;
  int info = -1;
  int n = *nn;

  clock_t t1 = clock();

  if (strcmp(pt->cfp.nom_method, mot1) == 0)
    cfp_gsnl(vec, q, & n, & pt->cfp.mu, & pt->cfp.itermax, & pt->cfp.tol, z, w, &it_end, &res, &info);
  else if (strcmp(pt->cfp.nom_method, mot2) == 0)
  {
    cfp_gcp(vec, q, &n, &pt->cfp.mu, &pt->cfp.itermax, &pt->cfp.tol, z, w, &it_end, &res, &info);
  }
  else if (strcmp(pt->cfp.nom_method, mot3) == 0)
    cfp_latin(vec, q, &n, &pt->cfp.k_latin, &pt->cfp.mu, &pt->cfp.itermax, &pt->cfp.tol, z, w, &it_end, &res, &info);
  else printf("Warning : Unknown solving method : %s\n", pt->cfp.nom_method);


  clock_t t2 = clock();
  printf("%.4lf seconds of processing\n", (t2 - t1) / (double)CLOCKS_PER_SEC);


  return info;
}
