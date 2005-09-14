/*!\file dfc_2D_solver.c
 *
 * This subroutine allows the dual resolution of contact problems with friction\n
 *
 * Try \f$(z,w)\f$ such that:
 * \f$
 *  \left\lbrace
 *   \begin{array}{l}
 *    M z- w=q\\
 *    0 \le z_n \perp w_n \ge 0\\
 *    -z_t \in \partial\psi_{[-\mu w_n, \mu w_n]}(w_t)\\
 *   \end{array}
 *  \right.
 * \f$
 *
 * here M is an n by n  matrix, q an n-dimensional vector, z an n-dimensional  vector and w an n-dimensional vector.
 *
 * This system of equations and inequalities is solved thanks to @ref dfc solvers or thanks to @ref lcp routines after
 * a new formulation of this problem in the LCP form due to the dfc_2D2lcp.c and lcp2dfc_2D.c routines.
 * The routine's call is due to the function dfc_2D_solver.c.
 *
 * \fn  int dfc_2D_solver( double *vec , double *q , int *n ,method *pt , double *z , double *w )
 *
 *  dfc_2D_solver is a generic interface allowing the call of one of the DFC solvers.
 *
 * \param vec          Unchanged parameter which contains the components of the DFC_2D matrix with a Fortran storage.
 * \param q            Unchanged parameter which contains the components of the constant right hand side vector.
 * \param nn           Unchanged parameter which represents the dimension of the DFC_2D problem.
 * \param pt           Unchanged parameter which represents the DFC_2D structure.
 * \param z            Modified parameter which contains the initial value of the LCP and returns the solution of the problem.
 * \param w            Modified parameter which returns the complementary solution of the problem.
 *
 * \return integer     0 - successful\n
 *                     0 >  - otherwise (see specific solvers for more information about the log info)
 *
 *  \author Nineb Sheherazade.& Mathieu Renouf
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifndef MEXFLAG
#include "solverpack.h"
#endif

int dfc_2D_solver(double *K1, double *F1, int *n, method *pt, double *U2 , double *F2)
{

  const char mot1[10] = "Latin", mot2[10] = "Lemke", mot3[10] = "NLGS", mot4[10] = "CPG";
  double res;
  int info = -1, it_end;
  double info1 = -1, it_end1, itt;
  int aa, b = 3;
  int nn = *n, tail1, tail2, tail3;
  int *dim_nn, dim_MM, dim_q, tempo = 87;
  int ddim_i, ddim_c, ddim_tt, ddim_n;
  int *ddl_i, *ddl_n, *ddl_tt, *ddl_c, dim_i, dim_tt, dim_n, dim_c, i;
  double *MM, *q, *z, *w;
  clock_t t1, t2;

  int itmp;
  double rtmp;

  int iparamLCP[5];
  int dparamLCP[5];

  for (i = 0 ; i < 5 ; ++i) iparamLCP[i] = 0;
  for (i = 0 ; i < 5 ; ++i) dparamLCP[i] = 0.0;

  if (strcmp(pt->dfc_2D.name, mot1) == 0)
  {
    dim_MM = 4 * pt->dfc_2D.dim_tt * pt->dfc_2D.dim_tt;
    dim_q = 2 * pt->dfc_2D.dim_tt;
    MM =  malloc(dim_MM * sizeof(double));
    q =  malloc(dim_q * sizeof(double));



    aa = dfc_2D2lcp(&tempo, & pt->dfc_2D.mu, pt, K1,  pt->dfc_2D.ddl_i, &pt->dfc_2D.dim_i, pt->dfc_2D.ddl_n, &pt->dfc_2D.dim_n,  pt->dfc_2D.ddl_tt, &pt->dfc_2D.dim_tt, pt->dfc_2D.ddl_c, &pt->dfc_2D.dim_c,  pt->dfc_2D.J1, F1, &nn, MM, q);


    z = (double *) malloc((2 * pt->dfc_2D.dim_tt) * sizeof(double));
    w = (double *) malloc((2 * pt->dfc_2D.dim_tt) * sizeof(double));



    t1 = clock();

    dfc_2D_latin(MM, q, &tempo, & pt->dfc_2D.k_latin, & pt->dfc_2D.mu, & pt->dfc_2D.itermax, & pt->dfc_2D.tol, z, w, & it_end, &res, &info);

    t2 = clock();

    printf("%.4lf seconds of processing\n", (t2 - t1) / (double)CLOCKS_PER_SEC);

    aa = lcp2dfc_2D(&tempo , z, w, pt, K1, F1, &nn, pt->dfc_2D.J1, pt->dfc_2D.ddl_i, &pt->dfc_2D.dim_i, pt->dfc_2D.ddl_c, &pt->dfc_2D.dim_c, pt->dfc_2D.ddl_n, pt->dfc_2D.ddl_tt, &pt->dfc_2D.dim_tt, U2, F2);

    free(MM);
    free(q);
    free(z);
    free(w);

  }

  else if (strcmp(pt->dfc_2D.name, mot2) == 0)
  {

    dim_MM = 9 * pt->dfc_2D.dim_tt * pt->dfc_2D.dim_tt;
    dim_q = 3 * pt->dfc_2D.dim_tt;
    MM =  malloc(dim_MM * sizeof(double));
    q =  malloc(dim_q * sizeof(double));

    aa = dfc_2D2lcp(&tempo, & pt->dfc_2D.mu, pt, K1,  pt->dfc_2D.ddl_i, &pt->dfc_2D.dim_i, pt->dfc_2D.ddl_n, &pt->dfc_2D.dim_n,  pt->dfc_2D.ddl_tt, &pt->dfc_2D.dim_tt, pt->dfc_2D.ddl_c, &pt->dfc_2D.dim_c,  pt->dfc_2D.J1, F1, &nn, MM, q);

    z = (double *) malloc((3 * pt->dfc_2D.dim_tt) * sizeof(double));
    w = (double *) malloc((3 * pt->dfc_2D.dim_tt) * sizeof(double));

    t1 = clock();

    lcp_lemke(MM, q, &tempo, & pt->dfc_2D.itermax, z, w, &it_end, &res, &info);

    t2 = clock();

    printf("%.4lf seconds of processing\n", (t2 - t1) / (double)CLOCKS_PER_SEC);

    aa = lcp2dfc_2D(&tempo , z, w, pt, K1, F1, &nn, pt->dfc_2D.J1, pt->dfc_2D.ddl_i, &pt->dfc_2D.dim_i, pt->dfc_2D.ddl_c, &pt->dfc_2D.dim_c, pt->dfc_2D.ddl_n, pt->dfc_2D.ddl_tt, &pt->dfc_2D.dim_tt, U2, F2);


    free(MM);
    free(q);
    free(z);
    free(w);



  }
  else if (strcmp(pt->dfc_2D.name, mot3) == 0)
  {

    dim_MM = 9 * pt->dfc_2D.dim_tt * pt->dfc_2D.dim_tt;
    dim_q = 3 * pt->dfc_2D.dim_tt;
    MM =  malloc(dim_MM * sizeof(double));
    q =  malloc(dim_q * sizeof(double));

    aa = dfc_2D2lcp(&tempo, & pt->dfc_2D.mu, pt, K1,  pt->dfc_2D.ddl_i, &pt->dfc_2D.dim_i, pt->dfc_2D.ddl_n, &pt->dfc_2D.dim_n,  pt->dfc_2D.ddl_tt, &pt->dfc_2D.dim_tt, pt->dfc_2D.ddl_c, &pt->dfc_2D.dim_c,  pt->dfc_2D.J1, F1, &nn, MM, q);

    z = (double *) malloc((3 * pt->dfc_2D.dim_tt) * sizeof(double));
    w = (double *) malloc((3 * pt->dfc_2D.dim_tt) * sizeof(double));

    t1 = clock();
    itmp = 0;
    rtmp = 1.0;

    iparamLCP[0] = pt->dfc_2D.itermax;
    iparamLCP[1] = 0;
    dparamLCP[0] = pt->dfc_2D.tol;
    dparamLCP[1] = 1.0;

    lcp_nlgs(&tempo , MM , q , z , w , &info , iparamLCP , dparamLCP);

    it_end = iparamLCP[2];
    res    = dparamLCP[2];

    t2 = clock();

    printf("%.4lf seconds of processing\n", (t2 - t1) / (double)CLOCKS_PER_SEC);

    aa = lcp2dfc_2D(&tempo , z, w, pt, K1, F1, &nn, pt->dfc_2D.J1, pt->dfc_2D.ddl_i, &pt->dfc_2D.dim_i, pt->dfc_2D.ddl_c, &pt->dfc_2D.dim_c, pt->dfc_2D.ddl_n, pt->dfc_2D.ddl_tt, &pt->dfc_2D.dim_tt, U2, F2);

    free(MM);
    free(q);
    free(z);
    free(w);

  }
  else if (strcmp(pt->dfc_2D.name, mot4) == 0)
  {

    dim_MM = (3 * pt->dfc_2D.dim_tt) * (3 * pt->dfc_2D.dim_tt);
    dim_q = (3 * pt->dfc_2D.dim_tt);
    MM =  malloc(dim_MM * sizeof(double));
    q =  malloc(dim_q * sizeof(double));

    aa = dfc_2D2lcp(&tempo, & pt->dfc_2D.mu, pt, K1,  pt->dfc_2D.ddl_i, &pt->dfc_2D.dim_i, pt->dfc_2D.ddl_n, &pt->dfc_2D.dim_n,  pt->dfc_2D.ddl_tt, &pt->dfc_2D.dim_tt, pt->dfc_2D.ddl_c, &pt->dfc_2D.dim_c,  pt->dfc_2D.J1, F1, &nn, MM, q);

    z = (double *) malloc((3 * pt->dfc_2D.dim_tt) * sizeof(double));
    w = (double *) malloc((3 * pt->dfc_2D.dim_tt) * sizeof(double));

    t1 = clock();
    itmp = 0;
    rtmp = 1.0;

    iparamLCP[0] = pt->dfc_2D.itermax;
    iparamLCP[1] = 0;
    dparamLCP[0] = pt->dfc_2D.tol;

    lcp_cpg(&tempo , MM , q , z , w , &info , iparamLCP , dparamLCP);

    it_end = iparamLCP[2];
    res    = dparamLCP[2];

    t2 = clock();
    printf("%.4lf seconds of processing\n", (t2 - t1) / (double)CLOCKS_PER_SEC);

    aa = lcp2dfc_2D(&tempo , z, w, pt, K1, F1, &nn, pt->dfc_2D.J1, pt->dfc_2D.ddl_i, &pt->dfc_2D.dim_i, pt->dfc_2D.ddl_c, &pt->dfc_2D.dim_c, pt->dfc_2D.ddl_n, pt->dfc_2D.ddl_tt, &pt->dfc_2D.dim_tt, U2, F2);

    free(MM);
    free(q);
    free(z);
    free(w);

  }
  else printf("Warning : Unknown solving method : %s\n", pt->dfc_2D.name);



  return info;
}

