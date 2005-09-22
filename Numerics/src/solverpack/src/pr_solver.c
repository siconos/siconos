/*!\file pr_solver.c

  This subroutine allows the primal resolution of relay problems.

  Try \f$(z,w)\f$ such that:
\f$
\left\lbrace
\begin{array}{l}
M z- w=q\\
-w \in \partial\psi_{[-a, a]}(z)\\
\end{array}
\right.
\f$

 here M is an n by n  matrix, q an n-dimensional vector, z an n-dimensional  vector and w an n-dimensional vector.
 This system of equations and inequalities is solved thanks to @ref pr solvers.
 The routine's call is due to the function pr_solver.c.

\fn int pr_solver(double vec[],double *q,int *nn, method *pt,double z[],double w[])

   pr_solver is a generic interface allowing the call of one of the PR solvers.

   \param double*  : vec On enter double vector containing the components of the double matrix with a fortran90 allocation.
   \param double*  : q On enter a pointer over doubles containing the components of the second member of the system.
   \param int*     : nn On enter a pointer over integers, the dimension of the second member.
   \param method*  : pt On enter a pointer other a structure (::method).
   \param double[] : z On return real vector, the solution of the problem.
   \param double[] : w On return real vector, the solution of the problem.

  \return  On return int, the termination reason (0 is successful otherwise 1).

   \author Nineb Sheherazade.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef MEXFLAG
#include "solverpack.h"
#endif
#include <time.h>

int pr_solver(double *vec, double *q, int *nn, method *pt, double z[], double w[])
{


  int info = -1, it_end;

  char prkey1[10] = "NLGS", prkey2[10] = "Latin";

  double res;

  clock_t t1, t2;


  t1 = clock();


  if (strcmp(pt->pr.name , prkey1) == 0)
  {
    pr_nlgs(vec, q, nn, pt->pr.a, pt->pr.b, & pt->pr.itermax, & pt->pr.tol, &pt->pr.chat, z, w, &it_end, &res, &info);

    pt->pr.err = res;
    pt->pr.iter = it_end;

  }
  else if (strcmp(pt->pr.name, prkey2) == 0)
  {
    pr_latin(vec, q, nn, &pt->pr.k_latin, pt->pr.a, pt->pr.b, &pt->pr.itermax, &pt->pr.tol, &pt->pr.chat, z, w, &it_end, &res, &info);

    pt->pr.err = res;
    pt->pr.iter = it_end;
  }

  else printf("Warning : Unknown solving method : %s\n", pt->pr.name);

  t2 = clock();


  printf("%.4lf seconds of processing\n", (t2 - t1) / (double)CLOCKS_PER_SEC);

  return info;
}
