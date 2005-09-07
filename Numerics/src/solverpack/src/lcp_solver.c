/*!\file lcp_solver.c
 *
 * This subroutine allows the resolution of LCP (Linear Complementary Problem).
 * Try \f$(z,w)\f$ such that:
 *
 * \f$
 *  \left\lbrace
 *   \begin{array}{l}
 *    0 \le z \perp Mz + q = w \ge 0\\
 *   \end{array}
 *  \right.
 * \f$
 *
 * M is an ( n x n ) matrix, q , w and z n-vector. This system of equalities and inequalities
 * is solved thanks to @ref lcp solvers. The routine's call is due to the function lcp_solver.c.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifndef MEXFLAG
#include "solverpack.h"
#endif

/*!\fn int lcp_solver( double *vec , double *q , int *nn , method *pt , double *z , double *w , int *it_end , double *res )
 *
 * lcp_solver is a generic interface allowing the call of one of the LCP solvers.
 *
 * \param double* vec  Unchanged parameter which contains the components of the LCP matrix with a Fortran storage.
 * \param double* q    Unchanged parameter which contains the components of the constant right hand side vector.
 * \param int* nn      Unchanged parameter which represents the dimension of the LCP problem.
 * \param method* pt  Unchanged parameter which represents the LCP structure.
 * \param double* z    Modified parameter which contains the initial value of the LCP and returns the solution of the problem.
 * \param double* w    Modified parameter which returns the complementary solution of the problem.
 * \param it_end       Modified parameter which returns the number of iterations performed by the algorithm.
 * \param res          Modified parameter which returns the final error value.
 *
 * \return integer     0 - successful / 1 - otherwise
 *
 * \author Nineb Sheherazade & Mathieu Renouf
 */

int lcp_solver(double *vec, double *q , int *n , method *pt , double *z , double *w , int *it_end , double *res)
{

  const char mot1[10] = "Lemke", mot2[10] = "NLGS", mot3[10] = "CPG";
  const char mot4[10] = "Latin", mot5[10] = "QP", mot6[10] = "NSQP";
  const char mot7[15] = "LexicoLemke";

  int info;

  *it_end = 0;
  *res    = 0.0;

  if (strcmp(pt->lcp.name , mot1) == 0)

    lcp_lemke(vec , q , n , &pt->lcp.itermax , z ,   /* in  */
              w , it_end , res , &info);            /* out */

  else if (strcmp(pt->lcp.name , mot2) == 0)

    lcp_nlgs(vec , q , n , &pt->lcp.itermax , &pt->lcp.tol , z , &pt->lcp.relax , &pt->lcp.iout ,  /* in  */
             w , it_end , res , &info);                                                           /* out */

  else if (strcmp(pt->lcp.name , mot3) == 0)

    lcp_cpg(vec , q , n , &pt->lcp.itermax , &pt->lcp.tol , z , &pt->lcp.iout ,  /* in  */
            w , it_end , res , &info);                                          /* out */

  else if (strcmp(pt->lcp.name , mot4) == 0)

    lcp_latin(vec , q , n , &pt->lcp.k_latin , &pt->lcp.itermax , &pt->lcp.tol , z ,   /* in  */
              w , it_end , res , &info);                                              /* out */

  else if (strcmp(pt->lcp.name , mot5) == 0)
  {

    // We assume that the LCP matrix M is symmetric

    lcp_qp(vec , q , n , &pt->lcp.tol , z ,  /* in  */
           w , &info);                      /* out */

  }
  else if (strcmp(pt->lcp.name , mot6) == 0)
  {

    // We assume that the LCP matrix M is not symmetric

    lcp_nsqp(vec , q , n , &pt->lcp.tol , z ,  /* in  */
             w , &info);                  /* out */
  }
  else if (strcmp(pt->lcp.name , mot7) == 0)

    lcp_lexicolemke(vec , q , n , &pt->lcp.itermax , z , &pt->lcp.iout ,  /* in  */
                    w , it_end , &info);                                 /* out */

  else printf("Warning : Unknown solver : %s\n", pt->lcp.name);

  return info;

}
