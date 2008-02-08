/* Siconos-Numerics version 2.1.1, Copyright INRIA 2005-2007.
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

#include "FrictionContact2D_Solvers.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifndef MEXFLAG
#include "NonSmoothDrivers.h"
#endif

int dfc_2D_driver(double *K1, double *F1, int *n, method *pt, double *U2 , double *F2 , double *mu)
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


    dfc_2D2cond_2D(n , mu , K1 , F1, pt->dfc_2D.ddl_n , pt->dfc_2D.ddl_tt , &pt->dfc_2D.dim_tt ,
                   pt->dfc_2D.ddl_d , &pt->dfc_2D.dim_d ,  pt->dfc_2D.J1 , MM , q);


    t1 = clock();


    dfc_2D_latin(MM , q , &dim_q , &pt->dfc_2D.k_latin , mu , &pt->dfc_2D.itermax ,
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

    dfc_2D2lcp(n , mu , K1 , F1, pt->dfc_2D.ddl_n , pt->dfc_2D.ddl_tt , &pt->dfc_2D.dim_tt ,
               pt->dfc_2D.ddl_d , &pt->dfc_2D.dim_d ,  pt->dfc_2D.J1 , MM , q);


    t1 = clock();


    iparamLCP[0] = pt->dfc_2D.itermax;
    iparamLCP[1] = 1;
    dparamLCP[0] = pt->dfc_2D.tol;
    dparamLCP[1] = 1.0;

    /* \WARNING TMP Comment */
    //    lcp_lexicolemke( &dim_q , MM , q , z , w , &info , iparamLCP , dparamLCP );

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

    dfc_2D2lcp(n , mu , K1 , F1, pt->dfc_2D.ddl_n , pt->dfc_2D.ddl_tt , &pt->dfc_2D.dim_tt ,
               pt->dfc_2D.ddl_d , &pt->dfc_2D.dim_d ,  pt->dfc_2D.J1 , MM , q);


    t1 = clock();


    iparamLCP[0] = pt->dfc_2D.itermax;
    iparamLCP[1] = 1;
    dparamLCP[0] = pt->dfc_2D.tol;
    dparamLCP[1] = 1.0;

    /* \WARNING TMP Comment */
    /*     lcp_pgs( &dim_q , MM , q , z , w , &info , iparamLCP , dparamLCP ); */

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


    dfc_2D2lcp(n , mu , K1 , F1, pt->dfc_2D.ddl_n , pt->dfc_2D.ddl_tt , &pt->dfc_2D.dim_tt ,
               pt->dfc_2D.ddl_d , &pt->dfc_2D.dim_d ,  pt->dfc_2D.J1 , MM , q);


    t1 = clock();


    iparamLCP[0] = pt->dfc_2D.itermax;
    iparamLCP[1] = 1;
    dparamLCP[0] = pt->dfc_2D.tol;
    dparamLCP[1] = 1.0;

    /* \WARNING TMP COMMENT */
    //    lcp_cpg( &dim_q , MM , q , z , w , &info , iparamLCP , dparamLCP );

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

