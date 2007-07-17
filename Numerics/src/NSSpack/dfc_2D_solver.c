/* Siconos-Numerics version 2.1.1, Copyright INRIA 2005-2006.
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
/*!\file dfc_2D_solver.c

This subroutine allows the dual resolution of contact problems with friction in the 2D case (DFC_2D).\n

\fn int dfc_2D_solver( double *K1, double *F1, int *n, method *pt, double *U2 , double *F2 )


\param K1           On enter, the stiffness, a vector of double (in which the components
of the matrix have a Fortran storage),

\param F1           On enter, the right hand side, a vector of double,

\param n            On enter, the dimension of the DFC_2D problem, an integer,

\param pt           0n enter, the union (::method) containing the DFC_2D structure,
in this structure there is the following parameters:\n
- char   name:      the name of the solver we want to use (on enter),
- int    itermax:   the maximum number of iteration required (on enter)
- double tol:       the tolerance required (on enter)
- double  mu:       the friction coefficient (on enter)
- int    *ddl_n:    contact in normal direction dof (not prescribed) (on enter)
- int    *ddl_tt:   contact in tangential direction dof (not prescribed) (on enter)
- int    *ddl_d:    prescribed dof (on enter)
- int    dim_tt:    dimension of ddl_tt (= dimension of ddl_n) (on enter)
- int    dim_d:     dimension of ddl_d (on enter)
- double *J1:       gap in normal contact direction (on enter)
- int    chat:      an integer that can make on or off the chattering (0=off, >0=on)(on enter)
- int    iter:      the number of iteration carry out (on return)
- double err:       the residue (on return) \n\n

This problem can be solved thanks to @ref dfc_2D solvers or thanks to @ref lcp solvers after:\n
- either a condensation makes thanks to dfc_2D2cond_2D.c and cond_2D2dfc_2D.c,
- or a new formulation of this problem in the LCP form due to the dfc_2D2lcp.c and lcp2dfc_2D.c routines.

\param U2           On return, the solution of the problem, vector of double.

\param F2           On return, the solution of the problem, vector of double.

\return integer     0 - successful\n
0 >  - otherwise (see specific solvers for more information about the log info)

\author Nineb Sheherazade

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifndef MEXFLAG
#include "NSSpack.h"
#endif

int dfc_2D_solver(double *K1, double *F1, int *n, method *pt, double *U2 , double *F2)
{


  int           i, info, it_end;
  int           dim_q, dim_MM;
  int           iparamLCP[5];
  double        dparamLCP[5];

  double        res;
  double        *MM, *q, *z, *w;

  clock_t       t1, t2;

  const char    dfckey1[10] = "Cfd_latin", dfckey2[10] = "Lemke";
  const char    dfckey3[10] = "PGS",  dfckey4[10] = "CPG";
  const char    dfckey13[10] = "NLGS";


  info = -1;

  for (i = 0 ; i < 5 ; ++i) iparamLCP[i] = 0;
  for (i = 0 ; i < 5 ; ++i) dparamLCP[i] = 0.0;



  if (strcmp(pt->dfc_2D.name , dfckey1) == 0)
  {

    dim_q  = 2 * pt->dfc_2D.dim_tt;
    dim_MM = dim_q * dim_q;

    q      = (double *)malloc(dim_q * sizeof(double));
    z      = (double *)malloc(dim_q * sizeof(double));
    w      = (double *)malloc(dim_q * sizeof(double));
    MM     = (double *)malloc(dim_MM * sizeof(double));

    for (i = 0; i < dim_q; i++)
    {
      q[i] = 0.0;
      z[i] = 0.0;
      w[i] = 0.0;


    }


    dfc_2D2cond_2D(n , &pt->dfc_2D.mu , K1 , F1, pt->dfc_2D.ddl_n , pt->dfc_2D.ddl_tt , &pt->dfc_2D.dim_tt ,
                   pt->dfc_2D.ddl_d , &pt->dfc_2D.dim_d ,  pt->dfc_2D.J1 , MM , q);


    t1 = clock();


    dfc_2D_latin(MM , q , &dim_q , &pt->dfc_2D.k_latin , &pt->dfc_2D.mu , &pt->dfc_2D.itermax ,
                 & pt->dfc_2D.tol , &pt->dfc_2D.chat, z , w , &it_end, &res , &info);

    t2 = clock();

    printf("%.4lf seconds of processing\n", (t2 - t1) / (double)CLOCKS_PER_SEC);


    cond_2D2dfc_2D(n , z , w ,  K1 , F1 , pt->dfc_2D.J1 ,  pt->dfc_2D.ddl_n , pt->dfc_2D.ddl_tt , &pt->dfc_2D.dim_tt ,
                   pt->dfc_2D.ddl_d , &pt->dfc_2D.dim_d, U2 , F2);



    free(MM);
    free(q);
    free(z);
    free(w);

  }
  else if (strcmp(pt->dfc_2D.name , dfckey2) == 0)
  {

    dim_q  = 3 * pt->dfc_2D.dim_tt;
    dim_MM = dim_q * dim_q;

    q      = (double *)malloc(dim_q * sizeof(double));
    z      = (double *)malloc(dim_q * sizeof(double));
    w      = (double *)malloc(dim_q * sizeof(double));
    MM     = (double *)malloc(dim_MM * sizeof(double));

    for (i = 0; i < dim_q; i++)
    {
      q[i] = 0.0;
      z[i] = 0.0;
      w[i] = 0.0;


    }

    dfc_2D2lcp(n , &pt->dfc_2D.mu , K1 , F1, pt->dfc_2D.ddl_n , pt->dfc_2D.ddl_tt , &pt->dfc_2D.dim_tt ,
               pt->dfc_2D.ddl_d , &pt->dfc_2D.dim_d ,  pt->dfc_2D.J1 , MM , q);


    t1 = clock();


    iparamLCP[0] = pt->dfc_2D.itermax;
    iparamLCP[1] = 1;
    dparamLCP[0] = pt->dfc_2D.tol;
    dparamLCP[1] = 1.0;

    lcp_lexicolemke(&dim_q , MM , q , z , w , &info , iparamLCP , dparamLCP);

    t2 = clock();

    printf("%.4lf seconds of processing\n", (t2 - t1) / (double)CLOCKS_PER_SEC);

    lcp2dfc_2D(n , z , w ,  K1 , F1 , pt->dfc_2D.J1 ,  pt->dfc_2D.ddl_n , pt->dfc_2D.ddl_tt , &pt->dfc_2D.dim_tt ,
               pt->dfc_2D.ddl_d , &pt->dfc_2D.dim_d, U2 , F2);



    free(MM);
    free(q);
    free(z);
    free(w);



  }
  else if (strcmp(pt->dfc_2D.name , dfckey3) == 0)
  {

    dim_q  = 3 * pt->dfc_2D.dim_tt;
    dim_MM = dim_q * dim_q;

    q      = (double *)malloc(dim_q * sizeof(double));
    z      = (double *)malloc(dim_q * sizeof(double));
    w      = (double *)malloc(dim_q * sizeof(double));
    MM     = (double *)malloc(dim_MM * sizeof(double));


    for (i = 0; i < dim_q; i++)
    {
      q[i] = 0.0;
      z[i] = 0.0;
      w[i] = 0.0;


    }

    dfc_2D2lcp(n , &pt->dfc_2D.mu , K1 , F1, pt->dfc_2D.ddl_n , pt->dfc_2D.ddl_tt , &pt->dfc_2D.dim_tt ,
               pt->dfc_2D.ddl_d , &pt->dfc_2D.dim_d ,  pt->dfc_2D.J1 , MM , q);


    t1 = clock();


    iparamLCP[0] = pt->dfc_2D.itermax;
    iparamLCP[1] = 1;
    dparamLCP[0] = pt->dfc_2D.tol;
    dparamLCP[1] = 1.0;

    lcp_pgs(&dim_q , MM , q , z , w , &info , iparamLCP , dparamLCP);

    it_end = iparamLCP[2];
    res    = dparamLCP[2];

    t2     = clock();

    printf("%.4lf seconds of processing\n", (t2 - t1) / (double)CLOCKS_PER_SEC);

    lcp2dfc_2D(n , z , w ,  K1 , F1 , pt->dfc_2D.J1 ,  pt->dfc_2D.ddl_n , pt->dfc_2D.ddl_tt , &pt->dfc_2D.dim_tt ,
               pt->dfc_2D.ddl_d , &pt->dfc_2D.dim_d, U2 , F2);


    free(MM);
    free(q);
    free(z);
    free(w);
  }
  else if (strcmp(pt->dfc_2D.name , dfckey13) == 0)
  {
    printf("Warning: NLGS method is obsolete. Use PGS instead.\n");

  }
  else if (strcmp(pt->dfc_2D.name , dfckey4) == 0)
  {


    dim_q  = 3 * pt->dfc_2D.dim_tt;
    dim_MM = dim_q * dim_q;

    q      = (double *)malloc(dim_q * sizeof(double));
    z      = (double *)malloc(dim_q * sizeof(double));
    w      = (double *)malloc(dim_q * sizeof(double));
    MM     = (double *)malloc(dim_MM * sizeof(double));

    for (i = 0; i < dim_q; i++)
    {
      q[i] = 0.0;
      z[i] = 0.0;
      w[i] = 0.0;


    }


    dfc_2D2lcp(n , &pt->dfc_2D.mu , K1 , F1, pt->dfc_2D.ddl_n , pt->dfc_2D.ddl_tt , &pt->dfc_2D.dim_tt ,
               pt->dfc_2D.ddl_d , &pt->dfc_2D.dim_d ,  pt->dfc_2D.J1 , MM , q);


    t1 = clock();


    iparamLCP[0] = pt->dfc_2D.itermax;
    iparamLCP[1] = 1;
    dparamLCP[0] = pt->dfc_2D.tol;
    dparamLCP[1] = 1.0;


    lcp_cpg(&dim_q , MM , q , z , w , &info , iparamLCP , dparamLCP);

    it_end = iparamLCP[2];
    res    = dparamLCP[2];

    t2     = clock();

    printf("%.4lf seconds of processing\n", (t2 - t1) / (double)CLOCKS_PER_SEC);

    lcp2dfc_2D(n , z , w ,  K1 , F1 , pt->dfc_2D.J1 ,  pt->dfc_2D.ddl_n , pt->dfc_2D.ddl_tt , &pt->dfc_2D.dim_tt ,
               pt->dfc_2D.ddl_d , &pt->dfc_2D.dim_d, U2 , F2);


    free(MM);
    free(q);
    free(z);
    free(w);




  }
  else printf(" Warning !! Solver name unknown : %s\n", pt->dfc_2D.name);

  return info;
}

