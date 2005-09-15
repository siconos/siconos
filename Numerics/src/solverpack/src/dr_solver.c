/*!\file dr_solver.c
 *
 * This subroutine allows the dual resolution of relay problems.
 *
 * Try \f$(z,w)\f$ such that:
 *  \f$
 *   \left\lbrace
 *    \begin{array}{l}
 *      M z + q = w\\
 *     -w \in \partial\psi_{[-b, a]}(z)\\
 *    \end{array}
 *   \right.
 *  \f$
 *
 * where M is an (n x n)-matrix, q , w and z n-vectors.\n
 *
 * This system of equations and inequalities is solved thanks to @ref dr solvers.
 * The routine's call is due to the function pr_solver.c.
 * \fn int pr_solver( double *vec , double *q ,int *nn , method *pt , double *z , double *w )
 *
 * dr_solver is a generic interface allowing the call of one of the DR solvers.
 *
 * \param vec      On enter double vector containing the components of the double matrix with a fortran90 allocation.
 * \param q        On enter a pointer over doubles containing the components of the second member of the system.
 * \param nn       On enter a pointer over integers, the dimension of the second member.
 * \param pt       On enter a pointer other a structure (::method).
 * \param z        On return real vector, the solution of the problem.
 * \param w        On return real vector, the solution of the problem.
 *
 * \return  On return int, the termination reason (0 is successful otherwise 1).
 *
 * \author Nineb Sheherazade & Mathieu Renouf.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifndef MEXFLAG
#include "solverpack.h"
#endif

int dr_solver(double *vec , double *q , int *nn , method *pt , double *z , double *w)
{

  int info = -1, it_end;
  char drkey1[10] = "NLGS" , drkey2[10] = "CPG" , drkey3[10] = "Latin";
  double res;
  int n = *nn;

  clock_t t1, t2;

  t1 = clock();

  if (strcmp(pt->dr.name , drkey3) == 0)

    dr_latin(vec , q , &n , &pt->dr.k_latin , pt->dr.a , pt->dr.b , &pt->dr.itermax , &pt->dr.tol , z , w , &it_end , &res , &info);

  else if (strcmp(pt->dr.name , drkey1) == 0)

    dr_nlgs(vec , q , &n , pt->dr.a , pt->dr.b , &pt->dr.itermax , &pt->dr.tol , z , w , &it_end , &res , &info);

  else printf("Warning : Unknown solving method : %s\n", pt->dr.name);

  t2 = clock();
  printf("%.4lf seconds of processing\n", (t2 - t1) / (double)CLOCKS_PER_SEC);

  return info;

}
