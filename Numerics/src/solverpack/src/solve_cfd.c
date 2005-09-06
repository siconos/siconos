/*!\file solve_cfd.c

  This subroutine allows the dual resolution of contact problems with friction.

  Try \f$(z,w)\f$ such that:
\f$
\left\lbrace
\begin{array}{l}
M z- w=q\\
0 \le z_n \perp w_n \ge 0\\
-z_t \in \partial\psi_{[-\mu w_n, \mu w_n]}(w_t)\\
\end{array}
\right.
\f$

 here M is an n by n  matrix, q an n-dimensional vector, z an n-dimensional  vector and w an n-dimensional vector.

 This system of equations and inequalities is solved thanks to @ref dfc solvers or thanks to @ref lcp routines after a new formulation of this problem in the LCP form due to the cfd_lcp.c and lcp_cfd.c routines.
 The routine's call is due to the function solve_cfd.c.
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifndef MEXFLAG
#include "solverpack.h"
#endif

/*!\fn  int solve_cfd (double *vec,double *q,int *n,methode *pt,double z[],double w[])

   solve_cfd is a generic interface allowing the call of one of the DFC solvers.
   \param vec On enter a double vector containing the components of the double matrix with a fortran90 allocation.
   \param q On enter a pointer over doubles containing the components of the double vector.
   \param n On enter a pointer over integers, the dimension of the second member.
   \param pt On enter a pointer over a structure (::methode).
   \param z On return double vector, the solution of the problem.
   \param w On return double vector, the solution of the problem.

   \return On return a pointer over integers, the termination reason (0 is successful otherwise 1).

   \author Nineb Sheherazade.
*/


int solve_cfd(double *K1, double *F1, int *n, methode *pt, double U2[], double F2[])
{

  const char mot1[10] = "Cfd_latin", mot2[10] = "Lemke", mot3[10] = "Gsnl", mot4[10] = "Gcp";
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


  if (strcmp(pt->cfd.nom_method, mot1) == 0)
  {
    dim_MM = 4 * pt->cfd.dim_tt * pt->cfd.dim_tt;
    dim_q = 2 * pt->cfd.dim_tt;
    MM =  malloc(dim_MM * sizeof(double));
    q =  malloc(dim_q * sizeof(double));



    aa = cfd_lcp(&tempo, & pt->cfd.mu, pt, K1,  pt->cfd.ddl_i, &pt->cfd.dim_i, pt->cfd.ddl_n, &pt->cfd.dim_n,  pt->cfd.ddl_tt, &pt->cfd.dim_tt, pt->cfd.ddl_c, &pt->cfd.dim_c,  pt->cfd.J1, F1, &nn, MM, q);


    z = (double *) malloc((2 * pt->cfd.dim_tt) * sizeof(double));
    w = (double *) malloc((2 * pt->cfd.dim_tt) * sizeof(double));



    t1 = clock();

    cfd_latin(MM, q, &tempo, & pt->cfd.k_latin, & pt->cfd.mu, & pt->cfd.itermax, & pt->cfd.tol, z, w, & it_end, &res, &info);

    t2 = clock();

    printf("%.4lf seconds of processing\n", (t2 - t1) / (double)CLOCKS_PER_SEC);

    aa = lcp_cfd(&tempo , z, w, pt, K1, F1, &nn, pt->cfd.J1, pt->cfd.ddl_i, &pt->cfd.dim_i, pt->cfd.ddl_c, &pt->cfd.dim_c, pt->cfd.ddl_n, pt->cfd.ddl_tt, &pt->cfd.dim_tt, U2, F2);

    free(MM);
    free(q);
    free(z);
    free(w);

  }

  else if (strcmp(pt->cfd.nom_method, mot2) == 0)
  {

    dim_MM = 9 * pt->cfd.dim_tt * pt->cfd.dim_tt;
    dim_q = 3 * pt->cfd.dim_tt;
    MM =  malloc(dim_MM * sizeof(double));
    q =  malloc(dim_q * sizeof(double));

    aa = cfd_lcp(&tempo, & pt->cfd.mu, pt, K1,  pt->cfd.ddl_i, &pt->cfd.dim_i, pt->cfd.ddl_n, &pt->cfd.dim_n,  pt->cfd.ddl_tt, &pt->cfd.dim_tt, pt->cfd.ddl_c, &pt->cfd.dim_c,  pt->cfd.J1, F1, &nn, MM, q);

    z = (double *) malloc((3 * pt->cfd.dim_tt) * sizeof(double));
    w = (double *) malloc((3 * pt->cfd.dim_tt) * sizeof(double));

    t1 = clock();

    lemke_lcp(MM, q, &tempo, & pt->cfd.itermax, z, w, &it_end, &res, &info);

    t2 = clock();

    printf("%.4lf seconds of processing\n", (t2 - t1) / (double)CLOCKS_PER_SEC);

    aa = lcp_cfd(&tempo , z, w, pt, K1, F1, &nn, pt->cfd.J1, pt->cfd.ddl_i, &pt->cfd.dim_i, pt->cfd.ddl_c, &pt->cfd.dim_c, pt->cfd.ddl_n, pt->cfd.ddl_tt, &pt->cfd.dim_tt, U2, F2);


    free(MM);
    free(q);
    free(z);
    free(w);



  }
  else if (strcmp(pt->cfd.nom_method, mot3) == 0)
  {

    dim_MM = 9 * pt->cfd.dim_tt * pt->cfd.dim_tt;
    dim_q = 3 * pt->cfd.dim_tt;
    MM =  malloc(dim_MM * sizeof(double));
    q =  malloc(dim_q * sizeof(double));

    aa = cfd_lcp(&tempo, & pt->cfd.mu, pt, K1,  pt->cfd.ddl_i, &pt->cfd.dim_i, pt->cfd.ddl_n, &pt->cfd.dim_n,  pt->cfd.ddl_tt, &pt->cfd.dim_tt, pt->cfd.ddl_c, &pt->cfd.dim_c,  pt->cfd.J1, F1, &nn, MM, q);

    z = (double *) malloc((3 * pt->cfd.dim_tt) * sizeof(double));
    w = (double *) malloc((3 * pt->cfd.dim_tt) * sizeof(double));

    t1 = clock();
    itmp = 0;
    rtmp = 1.0;
    gsnl_lcp(MM, q,  &tempo, & pt->cfd.itermax, & pt->cfd.tol, z, &rtmp , &itmp , w, &it_end, &res, &info);

    t2 = clock();

    printf("%.4lf seconds of processing\n", (t2 - t1) / (double)CLOCKS_PER_SEC);

    aa = lcp_cfd(&tempo , z, w, pt, K1, F1, &nn, pt->cfd.J1, pt->cfd.ddl_i, &pt->cfd.dim_i, pt->cfd.ddl_c, &pt->cfd.dim_c, pt->cfd.ddl_n, pt->cfd.ddl_tt, &pt->cfd.dim_tt, U2, F2);


    free(MM);
    free(q);
    free(z);
    free(w);

  }
  else if (strcmp(pt->cfd.nom_method, mot4) == 0)
  {

    dim_MM = (3 * pt->cfd.dim_tt) * (3 * pt->cfd.dim_tt);
    dim_q = (3 * pt->cfd.dim_tt);
    MM =  malloc(dim_MM * sizeof(double));
    q =  malloc(dim_q * sizeof(double));

    aa = cfd_lcp(&tempo, & pt->cfd.mu, pt, K1,  pt->cfd.ddl_i, &pt->cfd.dim_i, pt->cfd.ddl_n, &pt->cfd.dim_n,  pt->cfd.ddl_tt, &pt->cfd.dim_tt, pt->cfd.ddl_c, &pt->cfd.dim_c,  pt->cfd.J1, F1, &nn, MM, q);

    z = (double *) malloc((3 * pt->cfd.dim_tt) * sizeof(double));
    w = (double *) malloc((3 * pt->cfd.dim_tt) * sizeof(double));

    t1 = clock();

    gcp_lcp(MM, q, &tempo, & pt->cfd.itermax, & pt->cfd.tol, z, w, &it_end, &res, &info);

    t2 = clock();
    printf("%.4lf seconds of processing\n", (t2 - t1) / (double)CLOCKS_PER_SEC);

    aa = lcp_cfd(&tempo , z, w, pt, K1, F1, &nn, pt->cfd.J1, pt->cfd.ddl_i, &pt->cfd.dim_i, pt->cfd.ddl_c, &pt->cfd.dim_c, pt->cfd.ddl_n, pt->cfd.ddl_tt, &pt->cfd.dim_tt, U2, F2);

    free(MM);
    free(q);
    free(z);
    free(w);

  }
  else printf("Warning : Unknown solving method : %s\n", pt->cfd.nom_method);



  return info;
}

