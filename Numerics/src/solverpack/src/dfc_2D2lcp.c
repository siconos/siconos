/*!\file dfc_2D2lcp.c
 *
 * This file allows the formulation in the LCP (Linear  Complementary Problem) form of a contact problem with friction.
 *
 *
 *
 * \fn void dfc_2D2lcp( int *dim_n , double *mumu , methode *pt , double *K1 , int *ddl_i , int *dim_i ,
 *                      int *ddl_n , int *dim_nn , int *ddl_tt , int *dim_tt , int *ddl_c , int *dim_c ,
 *      double *J1 , double *F1 ,int *dim_F1 , double *MM , double *q )
 *
 *  dfc_2D2lcp subroutine allows the formulation in the LCP (Linear  Complementary Problem) form of a contact problem with friction.
 *
 * \param dim_n   Modified parameter which returns the dimension of the matrix after the reformulation of the problem.
 * \param mumu    Unchanged parameter which represents the friction coefficient.
 * \param pt      Unchanged parameter which represents union methode.
 * \param K1      Unchanged parameter which contains the components of the rigidity matrix with a standard fortran allocation.
 * \param ddl_i   Unchanged parameter which contains the components ddl imposed.
 * \param dim_i   Unchanged parameter which contains the dimension of the vector ddl_i.
 * \param ddl_n   Unchanged parameter which contains the components ddl .
 * \param dim_nn  Unchanged parameter which contains the dimension of the vector ddl_n.
 * \param ddl_tt  Unchanged parameter which contains the components ddl .
 * \param dim_tt  Unchanged parameter which contains the dimension of the vector ddl_tt.
 * \param ddl_c   Unchanged parameter which contains the components ddl .
 * \param dim_c   Unchanged parameter which contains the dimension of the vector ddl_c.
 * \param J1      Unchanged parameter which contains the components of the free motion.
 * \param F1      Unchanged parameter which contains the effort vector.
 * \param dim_F1  Unchanged parameter which contains the dimension of the vector F1.
 * \param MM      Modified parameter which returns the components of a double matrix (dim_n,dim_n) with a fortran allocation.
 * \param q       Modified parameter which returns the right-hand-side vector (dim_n).
 *
 * \author Nineb Sheherazade & Mathieu Renouf.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef MEXFLAG
#include "solverpack.h"
#endif


void dfc_2D2lcp(int *dim_n , double *mumu , method *pt , double *K1 , int *ddl_i , int *dim_i ,
                int * ddl_n , int *dim_nn, int *ddl_tt, int *dim_tt, int *ddl_c , int *dim_c ,
                double *J1 , double *F1 , int *dim_F1 , double *MM , double *q)
{

  int i, j, size_i, size_n, size_tt, size_c, size_F1, kk;
  int info2, nnn, ind1, ind2;
  double mu = *mumu, invKii0;
  double **M, **Mlat;
  double **Kin, **Kit, **Knn, **Knt, **Ktn, **Ktt, **R;
  double **Kic, **Kci, **Kcc, **Kni, **Kti, **Knn_bis, **Ktn_bis, **Knt_bis, **Ktt_bis;
  double **Kii2, **invKii, **temp_ni, **temp_ti, **temp_nn, *Fi, **tempo_tn;
  double **temp_tt, **temp_tn, **temp_nt;
  double *qn, *qt, *qbis, *Jcn , **temp_ci, **temp_ic, **temp_cc, *qi, *Jc;
  double *q1, *q0, *q2, *q3;
  double r3, rrr, rr, r1, r2;
  char uplo = 'U';

  const char mot1[10] = "Latin", mot2[10] = "NLGS", mot3[10] = "CPG",  mot4[10] = "Lemke";

  size_i  = *dim_i;
  size_n  = *dim_nn;
  size_c  = *dim_c;
  size_tt = *dim_tt;
  size_F1 = *dim_F1;

  q0 = (double*)malloc(size_n * sizeof(double));
  q1 = (double*)malloc(size_n * sizeof(double));
  q2 = (double*)malloc(size_n * sizeof(double));
  q3 = (double*)malloc(size_n * sizeof(double));


  R       = (double**)malloc(size_i * sizeof(double*));
  temp_ic = (double**)malloc(size_i * sizeof(double*));
  invKii  = (double**)malloc(size_i * sizeof(double*));
  Kii2    = (double**)malloc(size_i * sizeof(double*));
  Kin     = (double**)malloc(size_i * sizeof(double*));
  Kit     = (double**)malloc(size_i * sizeof(double*));
  Kic     = (double**)malloc(size_i * sizeof(double*));

  for (i = 0 ; i < size_i ; ++i)
  {
    R[i]       = (double*)malloc(size_i * sizeof(double));
    invKii[i]  = (double*)malloc(size_i * sizeof(double));
    Kii2[i]    = (double*)malloc(size_i * sizeof(double));
    Kin[i]     = (double*)malloc(size_n * sizeof(double));
    temp_ic[i] = (double*)malloc(size_c * sizeof(double));
    Kit[i]     = (double*)malloc(size_tt * sizeof(double));
    Kic[i]     = (double*)malloc(size_c * sizeof(double));
  }


  temp_ni = (double**)malloc(size_n * sizeof(double*));
  qn      = (double*)malloc(size_n * sizeof(double*));
  Knn     = (double**)malloc(size_n * sizeof(double*));
  Knn_bis = (double**)malloc(size_n * sizeof(double*));
  Knt     = (double**)malloc(size_n * sizeof(double*));
  Knt_bis = (double**)malloc(size_n * sizeof(double*));
  Kni     = (double**)malloc(size_n * sizeof(double*));
  temp_nt = (double**)malloc(size_n * sizeof(double*));
  temp_nn = (double**)malloc(size_n * sizeof(double*));

  for (i = 0 ; i < size_n; ++i)
  {
    temp_ni[i] = (double*)malloc(size_i * sizeof(double));
    Knn[i]     = (double*)malloc(size_n * sizeof(double));
    Knn_bis[i] = (double*)malloc(size_n * sizeof(double));
    Knt[i]     = (double*)malloc(size_tt * sizeof(double));
    Knt_bis[i] = (double*)malloc(size_tt * sizeof(double));
    Kni[i]     = (double*)malloc(size_i * sizeof(double));
    temp_nt[i] = (double*)malloc(size_tt * sizeof(double));
    temp_nn[i] = (double*)malloc(size_n * sizeof(double));
  }


  qt       = (double*) malloc(size_tt * sizeof(double*));
  Ktn      = (double**) malloc(size_tt * sizeof(double*));
  Ktn_bis  = (double**) malloc(size_tt * sizeof(double*));
  Ktt      = (double**) malloc(size_tt * sizeof(double*));
  Ktt_bis  = (double**) malloc(size_tt * sizeof(double*));
  Kti      = (double**) malloc(size_tt * sizeof(double*));
  temp_tt  = (double**) malloc(size_tt * sizeof(double*));
  temp_ti  = (double**) malloc(size_tt * sizeof(double*));
  temp_tn  = (double**) malloc(size_tt * sizeof(double*));
  tempo_tn = (double**) malloc(size_tt * sizeof(double*));

  for (i = 0 ; i < size_tt ; i++)
  {

    Ktn[i]      = (double*) malloc(size_n * sizeof(double));
    Ktn_bis[i]  = (double*) malloc(size_n * sizeof(double));
    Ktt[i]      = (double*) malloc(size_tt * sizeof(double));
    Ktt_bis[i]  = (double*) malloc(size_tt * sizeof(double));
    Kti[i]      = (double*) malloc(size_i * sizeof(double));
    temp_tt[i]  = (double*) malloc(size_tt * sizeof(double));
    temp_ti[i]  = (double*) malloc(size_i * sizeof(double));
    temp_tn[i]  = (double*) malloc(size_n * sizeof(double));
    tempo_tn[i] = (double*) malloc(size_n * sizeof(double));

  }

  Kcc     = (double**) malloc(size_c * sizeof(double*));
  Kci     = (double**) malloc(size_c * sizeof(double*));
  temp_ci = (double**) malloc(size_c * sizeof(double*));
  temp_cc = (double**) malloc(size_c * sizeof(double*));

  for (i = 0 ; i < size_c ; ++i)
  {

    Kcc[i]     = (double*) malloc(size_c * sizeof(double));
    Kci[i]     = (double*) malloc(size_i * sizeof(double));
    temp_ci[i] = (double*) malloc(size_i * sizeof(double));
    temp_cc[i] = (double*) malloc(size_c * sizeof(double));

  }

  Fi  = (double*) malloc(size_i * sizeof(double));
  Jc  = (double*) malloc(size_c * sizeof(double));
  Jcn = (double*) malloc(size_n * sizeof(double));

  for (i = 0 ; i < size_i ; ++i)
  {

    ind1 = ddl_i[i];
    Fi[i] = F1[ind1];

    for (j = 0 ; j < size_i  ; ++j) Kii2[j][i] = K1[ ind1 + size_F1 * ddl_i[j] ];
    for (j = 0 ; j < size_n  ; ++j) Kin[i][j]  = K1[ ind1 + size_F1 * ddl_n[j] ];
    for (j = 0 ; j < size_c  ; ++j) Kic[i][j]  = K1[ ind1 + size_F1 * ddl_c[j] ];
    for (j = 0 ; j < size_tt ; ++j) Kit[i][j]  = K1[ ind1 + size_F1 * ddl_tt[j] ];
  }

  for (i = 0 ; i < size_n ; ++i)
  {

    ind1   = ddl_n[i];
    Jcn[i] = J1[ind1];

    for (j = 0 ; j < size_i  ; ++j) Kni[i][j] = K1[ ind1 + size_F1 * ddl_i[j] ];
    for (j = 0 ; j < size_n  ; ++j) Knn[i][j] = K1[ ind1 + size_F1 * ddl_n[j] ];
    for (j = 0 ; j < size_tt ; ++j)  Knt[i][j] = K1[ ind1 + size_F1 * ddl_tt[j] ];
  }

  for (i = 0 ; i < size_tt ; ++i)
  {
    ind1 = ddl_tt[i];
    for (j = 0 ; j < size_i  ; ++j) Kti[i][j] = K1[ ind1 + size_F1 * ddl_i[j] ];
    for (j = 0 ; j < size_n  ; ++j) Ktn[i][j] = K1[ ind1 + size_F1 * ddl_n[j] ];
    for (j = 0 ; j < size_tt ; ++j) Ktt[i][j] = K1[ ind1 + size_F1 * ddl_tt[j] ];
  }

  for (i = 0 ; i < size_c ; ++i)
  {

    ind1  = ddl_c[i];
    Jc[i] = J1[ind1];

    for (j = 0 ; j < size_c ; ++j) Kcc[i][j] = K1[ ind1 + size_F1 * ddl_c[j] ];
    for (j = 0 ; j < size_i ; ++j) Kci[i][j] = K1[ ind1 + size_F1 * ddl_i[j] ];
  }

  for (i = 0 ; i < size_i ; ++i)
    for (j = 0 ; j < size_i ; ++j)
      R[i][j] = 0.;

  /*Cholesky*/

  R[0][0] = sqrt(Kii2[0][0]);


  for (i = 1; i < size_i; i++)
  {
    rr = 0.0;
    rrr = 0.0;
    for (j = 0; j <= i; j++)
    {
      r1 = 0.0;
      r2 = 0.0;
      for (kk = 0; kk <= j - 1; kk++)
      {
        rr = R[i][kk] * R[j][kk] + r1;
        r1 = rr;
      }

      if (fabs(R[j][j]) <= 1.e-10)
      {
        printf("nul pivot %d ,and R(%d,%d) %g/n", j, j, j, R[j][j]);
        break;
      }
      R[i][j] = (1 / R[j][j]) * (Kii2[i][j] - rr);

      r3 = 0.0;
      for (kk = 0; kk <= i - 1; kk++)
      {
        rrr = R[i][kk] * R[i][kk] + r3;
        r3 = rrr;
      }

      R[i][i] = sqrt(Kii2[i][i] - rrr);
    }
  }

  /* End of Cholesky */

  dpotri_(&uplo, &size_i, R , &size_i, &info2);

  for (i = 0; i < size_i; i++)
    for (j = 0; j < size_i; j++)
    {
      invKii[j][i] = R[i][j];
      invKii[i][j] = R[i][j];
    }

  for (i = 0; i < size_n; i++)
  {
    for (j = 0; j < size_i; j++)
    {

      invKii0 = 0.0;
      for (kk = 0; kk < size_i; kk++)
      {

        temp_ni[i][j] = Kni[i][kk] * invKii[kk][j] + invKii0;
        invKii0       = temp_ni[i][j];
      }
    }
  }

  for (i = 0; i < size_n; i++)
  {
    for (j = 0; j < size_n; j++)
    {
      invKii0 = 0.0;
      for (kk = 0; kk < size_i; kk++)
      {
        temp_nn[i][j] = temp_ni[i][kk] * Kin[kk][j] + invKii0;
        invKii0 = temp_nn[i][j];
      }
    }
  }

  for (i = 0; i < size_n; i++) /* attention il faut la transposee */
    for (j = 0; j < size_n; j++)
    {
      Knn_bis[i][j] = Knn[j][i] - temp_nn[j][i];
    }

  for (i = 0; i < size_n; i++)
    for (j = 0; j < size_n; j++)
      Knn[i][j] = Knn_bis[i][j];

  for (i = 0; i < size_n; i++)
  {
    for (j = 0; j < size_tt; j++)
    {
      invKii0 = 0.0;
      for (kk = 0; kk < size_i; kk++)
      {
        temp_nt[i][j] = temp_ni[i][kk] * Kit[kk][j] + invKii0;
        invKii0 = temp_nt[i][j];
      }
    }
  }

  for (i = 0; i < size_n; i++) /* attention il faut la transposee */
    for (j = 0; j < size_tt; j++)
      Knt_bis[i][j] = Knt[j][i] - temp_nt[j][i];

  for (i = 0; i < size_n; i++)
    for (j = 0; j < size_tt; j++)
      Knt[i][j] = Knt_bis[j][i];


  for (i = 0; i < size_tt; i++)
  {
    for (j = 0; j < size_i; j++)
    {
      invKii0 = 0.0;
      for (kk = 0; kk < size_i; kk++)
      {
        temp_ti[i][j] = Kti[i][kk] * invKii[kk][j] + invKii0;
        invKii0 = temp_ti[i][j];
      }
    }
  }

  for (i = 0; i < size_tt; i++)
  {
    for (j = 0; j < size_n; j++)
    {
      invKii0 = 0.0;
      for (kk = 0; kk < size_i; kk++)
      {
        temp_tn[i][j] = temp_ti[i][kk] * Kin[kk][j] + invKii0;
        invKii0 = temp_tn[i][j];
      }
    }
  }

  for (i = 0; i < size_tt; i++) /* attention il faut la transposee */
    for (j = 0; j < size_n; j++)
      Ktn_bis[i][j] = Ktn[j][i] - temp_tn[j][i];

  for (i = 0; i < size_tt; i++)
    for (j = 0; j < size_n; j++)
      Ktn[i][j] = Ktn_bis[j][i];

  for (i = 0; i < size_tt; i++)
  {
    for (j = 0; j < size_tt; j++)
    {
      invKii0 = 0.0;
      for (kk = 0; kk < size_i; kk++)
      {
        temp_tt[i][j] = temp_ti[i][kk] * Kit[kk][j] + invKii0;
        invKii0 = temp_tt[i][j];
      }
    }
  }

  for (i = 0; i < size_tt; i++) /* attention il faut la transposee */
    for (j = 0; j < size_tt; j++)
      Ktt_bis[i][j] = Ktt[j][i] - temp_tt[j][i];

  for (i = 0; i < size_tt; i++)
    for (j = 0; j < size_tt; j++)
      Ktt[i][j] = Ktt_bis[j][i];


  /*       qn        */

  for (i = 0; i < size_n; i++)
  {
    for (j = 0; j < size_i; j++)
    {
      invKii0 = 0.0;
      for (kk = 0; kk < size_i; kk++)
      {
        temp_ni[i][j] = Kni[i][kk] * invKii[kk][j] + invKii0;
        invKii0 = temp_ni[i][j];
      }
    }
  }

  /*       qn         */

  for (i = 0; i < size_n; i++)
  {
    invKii0 = 0.0;
    for (j = 0; j < size_i; j++)
    {
      qn[i] = temp_ni[i][j] * Fi[j] + invKii0;
      invKii0 = qn[i];
    }
  }

  /*       qt        */

  for (i = 0; i < size_tt; i++)
  {
    for (j = 0; j < size_i; j++)
    {
      invKii0 = 0.0;
      for (kk = 0; kk < size_i; kk++)
      {
        temp_ti[i][j] = Kti[i][kk] * invKii[kk][j] + invKii0;
        invKii0 = temp_ti[i][j];
      }
    }
  }

  for (i = 0; i < size_tt; i++)
  {
    invKii0 = 0.0;
    for (j = 0; j < size_i; j++)
    {
      qt[i] = temp_ti[i][j] * Fi[j] + invKii0;
      invKii0 = qt[i];
    }
  }

  /*       qt         */


  if (strcmp(pt->dfc_2D.name, mot2) == 0)
  {

    M = (double**) malloc(3 * size_tt/*size_F1*/*sizeof(double));
    for (i = 0; i < 3 * size_tt/*F1*/; i++)
      M[i] = (double*) malloc(3 * size_tt/*F1*/*sizeof(double));


    for (i = 0; i < size_n; i++)
      for (j = 0; j < size_n; j++)
      {
        M[i][j] = Knn[i][j];
        M[i][size_n + j] = -Knt[i][j];
        M[i][2 * size_n + j] = Knt[i][j];
        M[size_n + i][j] = -Ktn[i][j] + mu * Knn[i][j];
        M[size_n + i][size_n + j] = Ktt[i][j] - mu * Knt[i][j];
        M[size_n + i][2 * size_n + j] = mu * Knt[i][j] - Ktt[i][j];
        M[2 * size_n + i][j] = Ktn[i][j] + mu * Knn[i][j];
        M[2 * size_n + i][size_n + j] = -Ktt[i][j] - mu * Knt[i][j];
        M[2 * size_n + i][2 * size_n + j] = mu * Knt[i][j] + Ktt[i][j];
      }



    for (i = 0; i < 3 * size_n; i++)
      for (j = 0; j < 3 * size_n; j++)
        MM[i * 3 * size_n + j] = 0.;

    for (i = 0; i < 3 * size_n; i++)
      for (j = 0; j < 3 * size_n; j++)
        MM[i * 3 * size_n + j] = M[i][j]; /*/fortran compatibility*/

    *dim_n = 3 * size_n;


    for (i = 0; i < size_n; i++)
    {
      /*/intf(" ici 12 3\n"); */
      invKii0 = 0.0;
      for (kk = 0; kk < size_n; kk++)
      {
        q2[i] = Knn[i][kk] * Jcn[kk] + invKii0;
        invKii0 = q2[i];
      }
      q1[i] = -q2[i] + qn[i];

    }


    nnn = size_n * size_n;

    for (i = 0; i < size_tt; i++)
    {
      for (j = 0; j < size_tt; j++)
      {
        tempo_tn[i][j] = mu * Knn[i][j] - Ktn[i][j];

      }
    }


    for (i = 0; i < size_n; i++)
    {

      invKii0 = 0.0;
      for (kk = 0; kk < size_n; kk++)
      {
        q0[i] = tempo_tn[i][kk] * Jcn[kk] + invKii0;
        invKii0 = q0[i];
      }
      q0[i] = -q0[i] + mu * qn[i];

    }



    for (i = 0; i < size_n; i++)
    {
      q2[i] = q0[i] - qt[i];
    }


    for (i = 0; i < size_tt; i++)
    {
      for (j = 0; j < size_tt; j++)
      {
        tempo_tn[i][j] = mu * Knn[i][j] + Ktn[i][j];
      }
    }


    for (i = 0; i < size_n; i++)
    {
      invKii0 = 0.0;
      for (kk = 0; kk < size_n; kk++)
      {
        q3[i] = tempo_tn[i][kk] * Jcn[kk] + invKii0;
        invKii0 = q3[i];
      }
      q3[i] = -q3[i] + mu * qn[i];

    }

    for (i = 0; i < size_n; i++)
    {
      q3[i] = q3[i] + qt[i];
    }

    for (i = 0; i < size_n; i++)
    {
      q[i] = -q1[i];
      q[size_n + i] = -q2[i];
      q[2 * size_n + i] = -q3[i];
    }


    for (i = 0; i < 3 * size_n; i++)
      free(M[i]);

    free(M);

  }
  else if (strcmp(pt->dfc_2D.name, mot3) == 0)
  {

    M = (double**) malloc(3 * size_tt/*F1*/*sizeof(double));
    for (i = 0; i < 3 * size_tt/*F1*/; i++)
      M[i] = (double*) malloc(3 * size_tt/*F1*/*sizeof(double));

    for (i = 0; i < size_n; i++)
      for (j = 0; j < size_n; j++)
      {
        M[i][j] = Knn[i][j];
        M[i][size_n + j] = -Knt[i][j];
        M[i][2 * size_n + j] = Knt[i][j];
        M[size_n + i][j] = -Ktn[i][j] + mu * Knn[i][j];
        M[size_n + i][size_n + j] = Ktt[i][j] - mu * Knt[i][j];
        M[size_n + i][2 * size_n + j] = mu * Knt[i][j] - Ktt[i][j];
        M[2 * size_n + i][j] = Ktn[i][j] + mu * Knn[i][j];
        M[2 * size_n + i][size_n + j] = -Ktt[i][j] - mu * Knt[i][j];
        M[2 * size_n + i][2 * size_n + j] = mu * Knt[i][j] + Ktt[i][j];
      }


    for (i = 0; i < 3 * size_n; i++)
      for (j = 0; j < 3 * size_n; j++)
        MM[i * 3 * size_n + j] = 0.;

    for (i = 0; i < 3 * size_n; i++)
      for (j = 0; j < 3 * size_n; j++)
        MM[i * 3 * size_n + j] = M[i][j]; /*/fortran compatibility*/

    *dim_n = 3 * size_n;

    for (i = 0; i < size_n; i++)
    {
      invKii0 = 0.0;
      for (kk = 0; kk < size_n; kk++)
      {
        q2[i] = Knn[i][kk] * Jcn[kk] + invKii0;
        invKii0 = q2[i];
      }
      q1[i] = -q2[i] + qn[i];

    }




    for (i = 0; i < size_tt; i++)
    {
      for (j = 0; j < size_tt; j++)
      {
        tempo_tn[i][j] = mu * Knn[i][j] - Ktn[i][j];

      }
    }


    for (i = 0; i < size_n; i++)
    {

      invKii0 = 0.0;
      for (kk = 0; kk < size_n; kk++)
      {
        q0[i] = tempo_tn[i][kk] * Jcn[kk] + invKii0;
        invKii0 = q0[i];
      }
      q0[i] = -q0[i] + mu * qn[i];

    }


    for (i = 0; i < size_n; i++)
    {
      q2[i] = q0[i] - qt[i];
    }


    for (i = 0; i < size_tt; i++)
    {
      for (j = 0; j < size_tt; j++)
      {
        tempo_tn[i][j] = mu * Knn[i][j] + Ktn[i][j];
      }
    }


    for (i = 0; i < size_n; i++)
    {
      invKii0 = 0.0;
      for (kk = 0; kk < size_n; kk++)
      {
        q3[i] = tempo_tn[i][kk] * Jcn[kk] + invKii0;
        invKii0 = q3[i];
      }
      q3[i] = -q3[i] + mu * qn[i];

    }

    for (i = 0; i < size_n; i++)
    {
      q3[i] = q3[i] + qt[i];
    }

    for (i = 0; i < size_n; i++)
    {
      q[i] = -q1[i];
      q[size_n + i] = -q2[i];
      q[2 * size_n + i] = -q3[i];
    }


    for (i = 0; i < 3 * size_n; i++)
      free(M[i]);

    free(M);


  }
  else if (strcmp(pt->dfc_2D.name, mot4) == 0)
  {

    M = (double**) malloc(3 * size_tt/*F1*/*sizeof(double));
    for (i = 0; i < 3 * size_tt/*F1*/; i++)
      M[i] = (double*) malloc(3 * size_tt/*F1*/*sizeof(double));

    for (i = 0; i < size_n; i++)
      for (j = 0; j < size_n; j++)
      {
        M[i][j] = Knn[i][j];
        M[i][size_n + j] = -Knt[i][j];
        M[i][2 * size_n + j] = Knt[i][j];
        M[size_n + i][j] = -Ktn[i][j] + mu * Knn[i][j];
        M[size_n + i][size_n + j] = Ktt[i][j] - mu * Knt[i][j];
        M[size_n + i][2 * size_n + j] = mu * Knt[i][j] - Ktt[i][j];
        M[2 * size_n + i][j] = Ktn[i][j] + mu * Knn[i][j];
        M[2 * size_n + i][size_n + j] = -Ktt[i][j] - mu * Knt[i][j];
        M[2 * size_n + i][2 * size_n + j] = mu * Knt[i][j] + Ktt[i][j];
      }


    for (i = 0; i < 3 * size_n; i++)
      for (j = 0; j < 3 * size_n; j++)
        MM[i * 3 * size_n + j] = 0.;

    for (i = 0; i < 3 * size_n; i++)
      for (j = 0; j < 3 * size_n; j++)
        MM[i * 3 * size_n + j] = M[i][j]; /*/fortran compatibility*/

    *dim_n = 3 * size_n;

    for (i = 0; i < size_n; i++)
    {
      invKii0 = 0.0;
      for (kk = 0; kk < size_n; kk++)
      {
        q2[i] = Knn[i][kk] * Jcn[kk] + invKii0;
        invKii0 = q2[i];
      }
      q1[i] = -q2[i] + qn[i];

    }




    for (i = 0; i < size_tt; i++)
    {
      for (j = 0; j < size_tt; j++)
      {
        tempo_tn[i][j] = mu * Knn[i][j] - Ktn[i][j];

      }
    }


    for (i = 0; i < size_n; i++)
    {

      invKii0 = 0.0;
      for (kk = 0; kk < size_n; kk++)
      {
        q0[i] = tempo_tn[i][kk] * Jcn[kk] + invKii0;
        invKii0 = q0[i];
      }
      q0[i] = -q0[i] + mu * qn[i];

    }


    for (i = 0; i < size_n; i++)
    {
      q2[i] = q0[i] - qt[i];
    }


    for (i = 0; i < size_tt; i++)
    {
      for (j = 0; j < size_tt; j++)
      {
        tempo_tn[i][j] = mu * Knn[i][j] + Ktn[i][j];
      }
    }


    for (i = 0; i < size_n; i++)
    {
      invKii0 = 0.0;
      for (kk = 0; kk < size_n; kk++)
      {
        q3[i] = tempo_tn[i][kk] * Jcn[kk] + invKii0;
        invKii0 = q3[i];
      }
      q3[i] = -q3[i] + mu * qn[i];

    }

    for (i = 0; i < size_n; i++)
    {
      q3[i] = q3[i] + qt[i];
    }

    for (i = 0; i < size_n; i++)
    {
      q[i] = -q1[i];
      q[size_n + i] = -q2[i];
      q[2 * size_n + i] = -q3[i];
    }


    for (i = 0; i < 3 * size_n; i++)
      free(M[i]);

    free(M);


  }

  else if (strcmp(pt->dfc_2D.name, mot1) == 0)
  {

    printf("dans mot1");

    for (i = 0; i < size_i; i++)
    {
      for (j = 0; j < size_c; j++)
      {
        invKii0 = 0.0;
        for (kk = 0; kk < size_i; kk++)
        {

          temp_ic[i][j] = invKii[i][kk] * Kic[kk][j] + invKii0;
          invKii0 = temp_ic[i][j];
        }
      }
    }




    for (i = 0; i < size_c; i++)
    {
      for (j = 0; j < size_c; j++)
      {
        invKii0 = 0.0;
        for (kk = 0; kk < size_i; kk++)
        {
          temp_cc[i][j] = Kci[i][kk] * temp_ic[kk][j] + invKii0;
          invKii0 = temp_cc[i][j];
        }
      }
    }

    Mlat = (double**) malloc(size_c * sizeof(double));
    for (i = 0; i < size_c/*F1*/; i++)
      Mlat[i] = (double*) malloc(size_c * sizeof(double));


    for (i = 0; i < size_c; i++)
      for (j = 0; j < size_c; j++)
        Mlat[i][j] = Kcc[i][j] - temp_cc[i][j];

    for (i = 0; i < 2 * size_tt; i++)
      for (j = 0; j < 2 * size_tt; j++)
        MM[i * 2 * size_tt + j] = 0.;

    for (i = 0; i < 2 * size_tt; i++)
      for (j = 0; j < 2 * size_tt; j++)
        MM[i + 2 * size_tt * j] = Mlat[i][j]; /*/fortran compatibility*/

    *dim_n = 2 * size_tt;

    /* reste q latin*/

    qi = (double*) malloc(size_i * sizeof(double));

    for (i = 0; i < size_i; i++)
    {
      invKii0 = 0.0;
      for (kk = 0; kk < size_i; kk++)
      {
        qi[i] = invKii[i][kk] * Fi[kk] + invKii0;
        invKii0 = qi[i];
      }

    }


    for (i = 0; i < size_c; i++)
    {
      invKii0 = 0.0;
      for (kk = 0; kk < size_i; kk++)
      {
        q[i] = Kci[i][kk] * qi[kk] + invKii0;
        invKii0 = q[i];
      }

    }


    qbis = (double*) malloc(size_c * sizeof(double));

    for (i = 0; i < size_c; i++)
    {
      invKii0 = 0.0;
      for (kk = 0; kk < size_c; kk++)
      {
        qbis[i] = Mlat[i][kk] * Jc[kk] + invKii0;
        invKii0 = qbis[i];
      }

    }


    for (i = 0 ; i < size_c ; ++i) q[i] -= qbis[i];

    for (i = 0 ; i < size_c ; ++i) free(Mlat[i]);

    free(Mlat);
    free(qi);
    free(qbis);

    printf("dans mot1 fin;");
  }
  else printf("Warning : Unknown solving method : %s\n", pt->dfc_2D.name);

  for (i = 0; i < size_i; i++)
  {
    free(temp_ic[i]);
    free(invKii[i]);
    free(Kii2[i]);
    free(Kin[i]);
    free(Kit[i]);
    free(Kic[i]);
  }
  free(temp_ic);
  free(invKii);
  free(Kii2);
  free(Kin);
  free(Kit);
  free(Kic);
  free(Fi);


  for (i = 0; i < size_tt; i++)
  {

    free(Knn[i]);
    free(Knn_bis[i]);
    free(Knt[i]);
    free(Knt_bis[i]);
    free(Kni[i]);
    free(temp_nt[i]);
    free(temp_nn[i]);
    free(Ktn[i]);
    free(Ktn_bis[i]);
    free(Ktt[i]);
    free(Ktt_bis[i]);
    free(Kti[i]);
    free(temp_tt[i]);
    free(temp_ti[i]);
    free(temp_tn[i]);
    free(tempo_tn[i]);
  }

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

  for (i = 0; i < size_c; i++)
  {
    free(Kcc[i]);
    free(Kci[i]);
    free(temp_ci[i]);
    free(temp_cc[i]);
  }

  free(Kcc);
  free(Kci);
  free(temp_ci);
  free(temp_cc);

  free(Jcn);
  free(Jc);

  free(R);
  free(temp_ni);


  return;
}
