
/*!\file solve_rd.c

  This subroutine allows the dual resolution of relay problems.

  Try \f$(z,w)\f$ such that:
\f$
\left\lbrace
\begin{array}{l}
M z- w=q\\
-z \in \partial\psi_{[-a,a]}(w)\\
\end{array}
\right.
\f$

 here M is an n by n  matrix, q an n-dimensional vector, z an n-dimensional  vector and w an n-dimensional vector.
 This system of equations and inequalities is solved thanks to @ref dr solvers.
 The routine's call is due to the function solve_rd.c.
*/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "SiconosNumerics.h"


/*!\fn   int solve_rd (double vec[],double *q,int *nn, methode *pt,double z[],double w[])

   solve_rd is a generic interface allowing the call of one of the DR solvers.

   \param double* : vec On enter double vector containing the components of the double matrix with a fortran90 allocation.
   \param double* : q On enter a pointer over doubles containing the components of the second member of the system.
   \param int* : nn On enter a pointer over integers, the dimension of the second member.
   \param methode* : pt On enter a pointer other a structure (::methode).
   \param double[] : z On return real vector, the solution of the problem.
   \param double[] : w On return real vector, the solution of the problem.

  \return On return integer, the termination reason (0 is successful otherwise 1).

   \author Nineb Sheherazade.
 */
int solve_rd(double *vec, double *q, int *nn, methode *pt, double z[], double w[])
{
  int info = -1, it_end, fail;
  char mot1[10] = "Gsnl", mot2[10] = "Gcp", mot3[10] = "Latin";
  double res;
  int n = *nn;

  if (strcmp(pt->rd.nom_method, mot1) == 0)
  {
    //    cfp_gsnl_(vec,q,& n,& pt->cfp.mu,& pt->cfp.itermax,& pt->cfp.tol,z,w,&it_end,&res,&info);
  }
  else if (strcmp(pt->rd.nom_method, mot2) == 0)
  {
    //    cfp_gcp_(vec,q,&n,& pt->cfp.mu,& pt->cfp.itermax,& pt->cfp.tol,z,w,&it_end,&res,&info);
  }
  else if (strcmp(pt->rd.nom_method, mot2) == 0)
  {
    //    rd_latin(vec,q,&n,& pt->rd.k_latin,pt->rd.a,& pt->rd.itermax,& pt->rd.tol,z,w,&it_end,&res,&info);
  }
  else printf("Warning : Unknown solving method : %s\n", pt->rd.nom_method);

  return info;
}
