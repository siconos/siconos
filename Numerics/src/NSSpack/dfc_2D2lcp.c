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

/*!\file dfc_2D2lcp.c

   This subroutine allows the formulation in the LCP (Linear  Complementary Problem) form
   of a 2D contact problem with friction.\n



   \fn void dfc_2D2lcp( int *dim_F1 , double *mu , double *K1 , double *F1, int *ddl_n ,
        int *ddl_tt ,int * dim_nc, int *ddl_d, int *dim_d , double *J1 , double *MM , double *q )




   \param dim_F1    On enter a pointer over integers, the dimension of the DFC_2D problem,

   \param mu      On enter a pointer over doubles, the friction coefficient,

   \param K1        On enter a pointer over doubles containing the components of the
                     rigidity matrix with a fortran90 storage,

   \param F1        On enter a pointer over doubles containing the right hand side,

   \param ddl_n     On enter a pointer over integers , the contact in normal direction dof
                     (not prescribed),

   \param ddl_tt    On enter a pointer over integers, the contact in tangential direction dof
                     (not prescribed)


   \param dim_nc    On enter a pointer over integers, the dimension of the vector ddl_tt.

   \param ddl_d     On enter a pointer over integers, the prescribed dof,

   \param dim_d     On enter a pointer over integers, the dimension of the vector ddl_d,

   \param J1        On enter a pointer over doubles, gap in normal contact direction.
   \n\n
   \param MM        On return a pointer over doubles containing the components of a double
                      matrix (3*dim_nc,3*dim_nc) with a fortran90 allocation.

   \param q         On return a pointer over doubles, a double vector (3*dim_nc).


   \author Nineb Sheherazade.

*/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef MEXFLAG
#include "NSSpack.h"
#endif
#include "LA.h"



void sortsn_(int *ddl_i, int *sort, int *n);


void dfc_2D2lcp(int *dim_F1 , double *mu , double *K1 , double *F1, int *ddl_n , int *ddl_tt , int * dim_nc, int *ddl_d, int *dim_d , double *J1 , double *MM , double *q)
{
  int            i, j, taille_i, taille_n, taille_tt, taille_c, taille_F1, kk, taille_sort;
  int            info2,  nnn, ind1, ind2, taille_d;
  int            *sort, *sort1, *sort2, *ddl_c, *vecF1, *vec_i, n3, *ddl_i, dim_nn, dim_i;
  integer incx = 1;

  double         alpha, invKii0;
  double         *M;
  double         *Kii, *Kin, *Kit, *Knn, *Knt, *Ktn, *Ktt;
  double         *Kni, *Kti, *Knn_bis, *Ktn_bis, *Knt_bis, *Ktt_bis;
  double         *Kii2, *invKii, *temp_ni, *temp_ti, *temp_nn, *Fi, *tempo_tn;
  double         *temp_tt, *temp_tn, *temp_nt, *R;
  double         *qn, *qt, *Jcn;
  double         *q1, *q0, *q2, *q3;

  /*  char           uplo='U';*/








  taille_n  = *dim_nc;
  taille_c  = 2 * taille_n;
  taille_tt = *dim_nc;
  taille_F1 = *dim_F1;
  taille_d  = *dim_d;
  dim_nn    = taille_n;


  ddl_c     = (int*) malloc(taille_c * sizeof(int));
  sort1     = (int*) malloc(taille_c * sizeof(int));
  sort      = (int*) malloc((taille_c + taille_d) * sizeof(int));
  sort2     = (int*) malloc((taille_c + taille_d) * sizeof(int));
  vecF1     = (int*) malloc(taille_F1 * sizeof(int));
  vec_i     = (int*) malloc(taille_F1 * sizeof(int));

  for (i = 0; i < taille_n; i++)
  {

    ddl_c[i]           = ddl_n[i];
    ddl_c[i + taille_n]  = ddl_tt[i];

  }


  sortsn_(ddl_c, sort1, &taille_c);



  for (i = 0; i < taille_c; i++)
  {

    sort[i] = ddl_c[i];

  }

  for (i = 0; i < taille_d; i++)
  {

    sort[i + taille_c] = ddl_d[i];

  }

  taille_sort = taille_d + taille_c;

  sortsn_(sort, sort2, &taille_sort);


  for (i = 0 ; i < taille_F1 ; i++)
  {

    vecF1[i] = i;

  }



  diffns(&taille_F1, vecF1, &taille_sort, sort, &taille_i, vec_i);


  ddl_i     = (int*) malloc(taille_i * sizeof(int));


  for (i = 0 ; i < taille_i ; i++)
  {

    ddl_i[i] = vec_i[i];

  }



  dim_i     = taille_i;


  q0        = (double*) malloc(dim_nn * sizeof(double));
  q1        = (double*) malloc(dim_nn * sizeof(double));
  q2        = (double*) malloc(dim_nn * sizeof(double));
  q3        = (double*) malloc(dim_nn * sizeof(double));

  R         = (double*) malloc(dim_i * dim_i * sizeof(double));
  temp_ni   = (double*) malloc(dim_nn * dim_i * sizeof(double));


  Kii       = (double*) malloc(taille_i * taille_i * sizeof(double));
  invKii    = (double*) malloc(taille_i * taille_i * sizeof(double));
  Kii2      = (double*) malloc(taille_i * taille_i * sizeof(double));
  Kin       = (double*) malloc(taille_i * taille_n * sizeof(double));
  Kit       = (double*) malloc(taille_i * taille_tt * sizeof(double));


  qn        = (double*) malloc(taille_n * sizeof(double));
  Knn       = (double*) malloc(taille_n * taille_n * sizeof(double));
  Knn_bis   = (double*) malloc(taille_n * taille_n * sizeof(double));
  Knt       = (double*) malloc(taille_n * taille_tt * sizeof(double));
  Knt_bis   = (double*) malloc(taille_n * taille_tt * sizeof(double));
  Kni       = (double*) malloc(taille_n * taille_i * sizeof(double));
  temp_nt   = (double*) malloc(taille_tt * taille_n * sizeof(double));
  temp_nn   = (double*) malloc(taille_n * taille_n * sizeof(double));




  qt        = (double*) malloc(taille_tt * sizeof(double));
  Ktn       = (double*) malloc(taille_tt * taille_n * sizeof(double));
  Ktn_bis   = (double*) malloc(taille_n * taille_tt * sizeof(double));
  Ktt       = (double*) malloc(taille_tt * taille_tt * sizeof(double));
  Ktt_bis   = (double*) malloc(taille_tt * taille_tt * sizeof(double));
  Kti       = (double*) malloc(taille_tt * taille_i * sizeof(double));
  temp_tt   = (double*) malloc(taille_tt * taille_tt * sizeof(double));
  temp_ti   = (double*) malloc(taille_i * taille_tt * sizeof(double));
  temp_tn   = (double*) malloc(taille_n * taille_tt * sizeof(double));
  tempo_tn  = (double*) malloc(taille_n * taille_tt * sizeof(double));

  M         = (double*) malloc(9 * taille_tt * taille_tt * sizeof(double));
  Fi        = (double*) malloc(taille_i * sizeof(double));
  Jcn       = (double*) malloc(taille_n * sizeof(double));




  for (i = 0; i < taille_i; i++)
  {
    Fi[i] = F1[ddl_i[i]];
    for (j = 0; j < taille_i; j++)
    {
      Kii[i + taille_i * j]  = 0.;
      Kii2[i + taille_i * j] = 0.;
    }
  }



  for (i = 0; i < taille_i; i++)
  {
    ind1 = ddl_i[i];
    for (j = 0; j < taille_i; j++)
    {
      ind2 = ddl_i[j];
      Kii[i + taille_i * j]  = K1[ind1 + taille_F1 * ind2];
      Kii2[i + taille_i * j] = Kii[i + taille_i * j];

    }

  }



  for (i = 0; i < taille_i * taille_n; i++)
    Kin[i] = 0.;

  for (i = 0; i < taille_i; i++)
  {
    for (j = 0; j < taille_n; j++)
    {
      Kin[i + taille_i * j] = K1[ddl_i[i] + taille_F1 * ddl_n[j]];
    }
  }



  for (i = 0; i < taille_i * taille_tt; i++)
    Kit[i] = 0.;

  for (i = 0; i < taille_i; i++)
    for (j = 0; j < taille_tt; j++)
    {
      Kit[i + taille_i * j] = K1[ddl_i[i] + taille_F1 * ddl_tt[j]];
    }


  for (i = 0; i < taille_n * taille_i; i++)
    Kni[i] = 0.;


  for (i = 0; i < taille_n; i++)
    for (j = 0; j < taille_i; j++)
    {
      Kni[i + taille_n * j] = K1[ddl_n[i] + taille_F1 * ddl_i[j]];
    }





  for (i = 0; i < taille_tt * taille_i; i++)
    Kti[i] = 0.;


  for (i = 0; i < taille_tt; i++)
    for (j = 0; j < taille_i; j++)
    {
      Kti[i + taille_tt * j] = K1[ddl_tt[i] + taille_F1 * ddl_i[j]];
    }

  for (i = 0; i < taille_n; i++)
    for (j = 0; j < taille_n; j++)
    {
      Knn[i + taille_n * j]     = 0.;
      temp_nn[i + taille_n * j] = 0.;

    }

  for (i = 0; i < taille_n; i++)
  {

    for (j = 0; j < taille_n; j++)
    {
      Knn[i + taille_n * j] = K1[ddl_n[i] + taille_F1 * ddl_n[j]];
    }

  }


  for (i = 0; i < taille_tt * taille_tt; i++)
    Ktt[i] = 0.;

  for (i = 0; i < taille_tt; i++)
    for (j = 0; j < taille_tt; j++)
    {
      Ktt[i + taille_tt * j] = K1[ddl_tt[i] + taille_F1 * ddl_tt[j]];
    }

  for (i = 0; i < taille_n; i++)
    for (j = 0; j < taille_tt; j++)
      Knt[i + taille_n * j] = 0.;



  for (i = 0;  i < taille_n; i++)
  {
    Jcn[i] = J1[ddl_n[i]];
    for (j = 0; j < taille_tt; j++)
      Knt[i + taille_n * j] = K1[ddl_n[i] + taille_F1 * ddl_tt[j]];
  }



  for (i = 0; i < taille_tt * taille_n; i++)
    Ktn[i] = 0.;

  for (i = 0; i < taille_tt; i++)
    for (j = 0; j < taille_n; j++)
      Ktn[i + taille_tt * j] = K1[ddl_tt[i] + taille_F1 * ddl_n[j]];



  for (i = 0; i < taille_i; i++)
    for (j = 0; j < taille_i; j++)
    {
      R[i + taille_i * j] = Kii2[i * taille_i + j];
    }




  /*                        Cholesky                                   */


  DPOTRF(LA_UP, taille_i, R , taille_i, info2);




  if (info2 != 0)
  {
    printf("Matter with Cholesky factorization \n");

    free(sort2);
    free(M);
    free(Kii);
    free(invKii);
    free(Kii2);
    free(Kin);
    free(Kit);
    free(Fi);
    free(qt);
    free(qn);
    free(Knn);
    free(Knn_bis);
    free(Knt);
    free(Knt_bis);
    free(Kni);
    free(temp_nt);
    free(temp_nn);
    free(Ktn);
    free(Ktn_bis);
    free(Ktt);
    free(Ktt_bis);
    free(Kti);
    free(temp_tt);
    free(temp_ti);
    free(temp_tn);
    free(tempo_tn);
    free(q0);
    free(q1);
    free(q2);
    free(q3);
    free(Jcn);
    free(R);
    free(temp_ni);


    free(sort);
    free(ddl_c);
    free(sort1);
    free(vecF1);
    free(vec_i);
    free(ddl_i);



  }


  DPOTRI(LA_UP, taille_i, R , taille_i, info2);

  if (info2 != 0)
  {
    printf("Matter with matrix inversion \n");

    free(sort2);
    free(M);
    free(Kii);
    free(invKii);
    free(Kii2);
    free(Kin);
    free(Kit);
    free(Fi);
    free(qt);
    free(qn);
    free(Knn);
    free(Knn_bis);
    free(Knt);
    free(Knt_bis);
    free(Kni);
    free(temp_nt);
    free(temp_nn);
    free(Ktn);
    free(Ktn_bis);
    free(Ktt);
    free(Ktt_bis);
    free(Kti);
    free(temp_tt);
    free(temp_ti);
    free(temp_tn);
    free(tempo_tn);
    free(q0);
    free(q1);
    free(q2);
    free(q3);
    free(Jcn);
    free(R);
    free(temp_ni);

    free(sort);
    free(ddl_c);
    free(sort1);
    free(vecF1);
    free(vec_i);
    free(ddl_i);


  }



  for (i = 0; i < taille_i; i++)
    for (j = 0; j < taille_i; j++)
    {
      invKii[j + taille_i * i] =  R[j + taille_i * i];

      invKii[i + taille_i * j] =  R[j + taille_i * i];

    }






  for (i = 0; i < taille_n; i++)
  {
    for (j = 0; j < taille_i; j++)
    {
      invKii0 = 0.0;
      for (kk = 0; kk < taille_i; kk++)
      {


        temp_ni[i + taille_n * j] = Kni[i + taille_n * kk] * invKii[kk + taille_i * j] + invKii0;

        invKii0 = temp_ni[i + taille_n * j];
      }
    }
  }




  for (i = 0; i < taille_n; i++)
  {
    for (j = 0; j < taille_n; j++)
    {
      invKii0 = 0.0;
      for (kk = 0; kk < taille_i; kk++)
      {
        temp_nn[i + taille_n * j] = temp_ni[i + taille_n * kk] * Kin[kk + taille_i * j] + invKii0;

        invKii0 = temp_nn[i + taille_n * j];
      }
    }
  }




  for (i = 0; i < taille_n; i++)
    for (j = 0; j < taille_n; j++)
    {
      Knn_bis[i + taille_n * j] = Knn[j + taille_n * i] - temp_nn[j + taille_n * i];
    }


  for (i = 0; i < taille_n; i++)
    for (j = 0; j < taille_n; j++)
      Knn[i + taille_n * j] = Knn_bis[i + taille_n * j];

  for (i = 0; i < taille_n; i++)
  {
    for (j = 0; j < taille_tt; j++)
    {
      invKii0 = 0.0;
      for (kk = 0; kk < taille_i; kk++)
      {
        temp_nt[i + taille_n * j] = temp_ni[i + taille_n * kk] * Kit[kk + taille_i * j] + invKii0;

        invKii0 = temp_nt[i + taille_n * j];
      }
    }
  }



  for (i = 0; i < taille_n; i++)
    for (j = 0; j < taille_tt; j++)
      Knt_bis[i + taille_n * j] = Knt[i + taille_n * j] - temp_nt[i + taille_n * j];





  for (i = 0; i < taille_n; i++)
    for (j = 0; j < taille_tt; j++)
      Knt[i + taille_n * j] = Knt_bis[i + taille_n * j];


  for (i = 0; i < taille_tt; i++)
  {
    for (j = 0; j < taille_i; j++)
    {
      invKii0 = 0.0;
      for (kk = 0; kk < taille_i; kk++)
      {
        temp_ti[i + taille_tt * j] = Kti[i + taille_tt * kk] * invKii[kk + taille_i * j] + invKii0;

        invKii0 = temp_ti[i + taille_tt * j];
      }
    }
  }


  for (i = 0; i < taille_tt; i++)
  {
    for (j = 0; j < taille_n; j++)
    {
      invKii0 = 0.0;
      for (kk = 0; kk < taille_i; kk++)
      {
        temp_tn[i + taille_tt * j] = temp_ti[i + taille_tt * kk] * Kin[kk + taille_i * j] + invKii0;

        invKii0 = temp_tn[i + taille_tt * j];
      }
    }
  }


  for (i = 0; i < taille_tt; i++)
    for (j = 0; j < taille_n; j++)
      Ktn_bis[i + taille_tt * j] = Ktn[j + taille_tt * i] - temp_tn[j + taille_tt * i];

  for (i = 0; i < taille_tt; i++)
    for (j = 0; j < taille_n; j++)
      Ktn[i + taille_tt * j] = Ktn_bis[j + taille_tt * i];

  for (i = 0; i < taille_tt; i++)
  {
    for (j = 0; j < taille_tt; j++)
    {
      invKii0 = 0.0;
      for (kk = 0; kk < taille_i; kk++)
      {
        temp_tt[i + taille_tt * j] = temp_ti[i + taille_tt * kk] * Kit[kk + taille_i * j] + invKii0;

        invKii0 = temp_tt[i + taille_tt * j];
      }
    }
  }



  for (i = 0; i < taille_tt; i++)
    for (j = 0; j < taille_tt; j++)
      Ktt_bis[i + taille_tt * j] = Ktt[j + taille_tt * i] - temp_tt[j + taille_tt * i];

  for (i = 0; i < taille_tt; i++)
    for (j = 0; j < taille_tt; j++)
      Ktt[i + taille_tt * j] = Ktt_bis[j + taille_tt * i];




  for (i = 0; i < taille_n; i++)
  {
    for (j = 0; j < taille_i; j++)
    {
      invKii0 = 0.0;
      for (kk = 0; kk < taille_i; kk++)
      {
        temp_ni[i + taille_n * j] = Kni[i + taille_n * kk] * invKii[kk + taille_i * j] + invKii0;

        invKii0 = temp_ni[i + taille_n * j];
      }
    }
  }





  for (i = 0; i < taille_n; i++)
  {
    invKii0 = 0.0;
    for (j = 0; j < taille_i; j++)
    {
      qn[i] = temp_ni[i + taille_n * j] * Fi[j] + invKii0;
      invKii0 = qn[i];
    }
  }






  for (i = 0; i < taille_tt; i++)
  {
    for (j = 0; j < taille_i; j++)
    {
      invKii0 = 0.0;
      for (kk = 0; kk < taille_i; kk++)
      {
        temp_ti[i + taille_tt * j] = Kti[i + taille_tt * kk] * invKii[kk + taille_i * j] + invKii0;

        invKii0 = temp_ti[i + taille_tt * j];
      }
    }
  }

  for (i = 0; i < taille_tt; i++)
  {
    invKii0 = 0.0;
    for (j = 0; j < taille_i; j++)
    {
      qt[i] = temp_ti[i + taille_tt * j] * Fi[j] + invKii0;

      invKii0 = qt[i];
    }
  }










  for (i = 0; i < taille_n; i++)
    for (j = 0; j < taille_n; j++)
    {
      M[i + 3 * taille_tt * j] = Knn[i + taille_n * j];

      M[i + 3 * taille_tt * (taille_n + j)] = -Knt[i + taille_n * j];

      M[i + 3 * taille_tt * (2 * taille_n + j)] = Knt[i + taille_n * j];

      M[taille_n + i + 3 * taille_tt * j] = -Ktn[i + taille_tt * j] + mu[i] * Knn[i + taille_n * j];

      M[taille_n + i + 3 * taille_tt * (taille_n + j)] = Ktt[i + taille_tt * j] - mu[i] * Knt[i + taille_n * j];

      M[taille_n + i + 3 * taille_tt * (2 * taille_n + j)] = mu[i] * Knt[i + taille_n * j] - Ktt[i + taille_tt * j];

      M[2 * taille_n + i + 3 * taille_tt * j] = Ktn[i + taille_tt * j] + mu[i] * Knn[i + taille_n * j];

      M[2 * taille_n + i + 3 * taille_tt * (taille_n + j)] = -Ktt[i + taille_tt * j] - mu[i] * Knt[i + taille_n * j];

      M[2 * taille_n + i + 3 * taille_tt * (2 * taille_n + j)] = mu[i] * Knt[i + taille_n * j] + Ktt[i + taille_tt * j];
    }



  for (i = 0; i < 3 * taille_n; i++)
    for (j = 0; j < 3 * taille_n; j++)
      MM[i * 3 * taille_n + j] = 0.;

  for (i = 0; i < 3 * taille_n; i++)
    for (j = 0; j < 3 * taille_n; j++)
      MM[j * 3 * taille_n + i] = M[i + 3 * taille_tt * j];

  n3 = 3 * taille_n;


  for (i = 0; i < taille_n; i++)
  {
    invKii0 = 0.0;
    for (kk = 0; kk < taille_n; kk++)
    {
      q2[i] = Knn[i + taille_n * kk] * Jcn[kk] + invKii0;
      invKii0 = q2[i];
    }
    q1[i] = -q2[i] + qn[i];

  }



  nnn = taille_n * taille_n;

  for (i = 0; i < taille_tt; i++)
  {
    for (j = 0; j < taille_tt; j++)
    {
      tempo_tn[i + taille_tt * j] = mu[i] * Knn[i + taille_n * j] - Ktn[i + taille_tt * j];

    }
  }


  for (i = 0; i < taille_n; i++)
  {

    invKii0 = 0.0;
    for (kk = 0; kk < taille_n; kk++)
    {
      q0[i] = tempo_tn[i + taille_tt * kk] * Jcn[kk] + invKii0;
      invKii0 = q0[i];
    }
    q0[i] = -q0[i] + mu[i] * qn[i];

  }



  for (i = 0; i < taille_n; i++)
  {
    q2[i] = q0[i] - qt[i];
  }


  for (i = 0; i < taille_tt; i++)
  {
    for (j = 0; j < taille_tt; j++)
    {
      tempo_tn[i + taille_tt * j] = mu[i] * Knn[i + taille_n * j] + Ktn[i + taille_tt * j];
    }
  }


  for (i = 0; i < taille_n; i++)
  {
    invKii0 = 0.0;

    for (kk = 0; kk < taille_n; kk++)
    {
      q3[i] = tempo_tn[i + taille_tt * kk] * Jcn[kk] + invKii0;

      invKii0 = q3[i];
    }

    q3[i] = -q3[i] + mu[i] * qn[i];

  }

  for (i = 0; i < taille_n; i++)
  {
    q3[i] = q3[i] + qt[i];
  }

  for (i = 0; i < taille_n; i++)
  {
    q[i] = -q1[i];

    q[taille_n + i] = -q2[i];

    q[2 * taille_n + i] = -q3[i];

  }



  alpha = -1.;
  DSCAL(n3, alpha, q, incx);



  free(M);
  free(Kii);
  free(invKii);
  free(Kii2);
  free(Kin);
  free(Kit);
  free(Fi);
  free(qt);
  free(qn);
  free(Knn);
  free(Knn_bis);
  free(Knt);
  free(Knt_bis);
  free(Kni);
  free(temp_nt);
  free(temp_nn);
  free(Ktn);
  free(Ktn_bis);
  free(Ktt);
  free(Ktt_bis);
  free(Kti);
  free(temp_tt);
  free(temp_ti);
  free(temp_tn);
  free(tempo_tn);
  free(q0);
  free(q1);
  free(q2);
  free(q3);
  free(Jcn);
  free(R);
  free(temp_ni);

  free(sort);
  free(sort2);
  free(ddl_c);
  free(sort1);
  free(vecF1);
  free(vec_i);
  free(ddl_i);

  /*        return ;*/
}
