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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "NSSTools.h"
#include "LA.h"

void dfc_2D2cond_2D(int *dim_F1, double *mu, double *K1, double *F1, int *ddl_n, int * ddl_tt, int *dim_nc, int *ddl_d, int *dim_d, double * J1, double *MM, double *q)
{


  int            i, j, taille_i, taille_n, taille_tt, taille_c, taille_F1, kk, taille_d;
  int            info2,  ind1, ind2;
  int            *sort, *sort1, *sort2, *ddl_c, *vecF1, *vec_i, *ddl_i, dim_nn, dim_i;
  int            taille_sort;
  int incx = 1;

  double         alpha, invKii0;
  double         *Mlat;
  double         *Kii;
  double         *Kic, *Kci, *Kcc;
  double         *Kii2, *invKii, *Fi;
  double         *R;
  double         *qbis, *Jcn , *temp_ic, *temp_cc, *qi, *Jc;
  double         *q1, *q0, *q2, *q3;

  taille_n  = *dim_nc;
  taille_c  = 2 * taille_n;
  taille_tt = *dim_nc;
  taille_F1 = *dim_F1;
  taille_d  = *dim_d;


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
  dim_nn    = taille_n;


  q0        = (double*) malloc(dim_nn * sizeof(double));
  q1        = (double*) malloc(dim_nn * sizeof(double));
  q2        = (double*) malloc(dim_nn * sizeof(double));

  q3        = (double*) malloc(dim_nn * sizeof(double));
  R         = (double*) malloc(dim_i * dim_i * sizeof(double));
  temp_ic   = (double*) malloc(taille_c * taille_i * sizeof(double));
  Kii       = (double*) malloc(taille_i * taille_i * sizeof(double));
  invKii    = (double*) malloc(taille_i * taille_i * sizeof(double));

  Kii2      = (double*) malloc(taille_i * taille_i * sizeof(double));
  Kic       = (double*) malloc(taille_i * taille_c * sizeof(double));
  Kcc       = (double*) malloc(taille_c * taille_c * sizeof(double));
  Kci       = (double*) malloc(taille_c * taille_i * sizeof(double));
  temp_cc   = (double*) malloc(taille_c * taille_c * sizeof(double));




  Fi        = (double*) malloc(taille_i * sizeof(double));

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





  Jcn = (double*) malloc(taille_n * sizeof(double));

  for (i = 0;  i < taille_n; i++)
  {
    Jcn[i] = J1[ddl_n[i]];

  }



  for (i = 0; i < taille_c * taille_c; i++)
    Kcc[i] = 0.;

  for (i = 0; i < taille_c; i++)
    for (j = 0; j < taille_c; j++)
    {
      Kcc[i + taille_c * j] = K1[ddl_c[i] + taille_F1 * ddl_c[j]];
    }



  for (i = 0; i < taille_i * taille_c; i++)
    Kic[i] = 0.;

  for (i = 0; i < taille_i; i++)
    for (j = 0; j < taille_c; j++)
      Kic[i + taille_i * j] = K1[ddl_i[i] + taille_F1 * ddl_c[j]];




  for (i = 0; i < taille_c * taille_i; i++)
    Kci[i] = 0.;

  Jc = (double*) malloc(taille_c * sizeof(double));

  for (i = 0; i < taille_c; i++)
  {
    Jc [i] = J1[ddl_c[i]];
    for (j = 0; j < taille_i; j++)
      Kci[i + taille_c * j] = K1[ddl_c[i] + taille_F1 * ddl_i[j]];
  }




  for (i = 0; i < taille_i; i++)
    for (j = 0; j < taille_i; j++)
    {
      R[i + taille_i * j] = Kii2[i + taille_i * j];
    }




  /*                        Cholesky                                   */


  DPOTRF(LA_UP , taille_i, R , taille_i, info2);



  if (info2 != 0)
  {
    printf("Matter with Cholesky factorization \n");


    free(temp_ic);
    free(Kii);

    free(invKii);
    free(Kii2);
    free(Kic);
    free(Fi);
    free(q0);

    free(q1);
    free(q2);
    free(q3);
    free(Kcc);
    free(Kci);

    free(temp_cc);
    free(Jcn);
    free(Jc);
    free(R);
    free(sort);

    free(sort2);
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


    free(temp_ic);
    free(Kii);
    free(invKii);
    free(Kii2);
    free(Kic);
    free(Fi);

    free(q0);
    free(q1);
    free(q2);
    free(q3);
    free(Kcc);
    free(Kci);
    free(temp_cc);
    free(Jcn);
    free(Jc);
    free(R);
    free(sort);
    free(sort2);
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








  for (i = 0; i < taille_i; i++)
  {
    for (j = 0; j < taille_c; j++)
    {
      invKii0 = 0.0;

      for (kk = 0; kk < taille_i; kk++)
      {

        temp_ic[i + taille_i * j] = invKii[i + taille_i * kk] * Kic[kk + taille_i * j] + invKii0;

        invKii0 = temp_ic[i + taille_i * j];
      }
    }
  }




  for (i = 0; i < taille_c; i++)
  {
    for (j = 0; j < taille_c; j++)
    {
      invKii0 = 0.0;

      for (kk = 0; kk < taille_i; kk++)
      {
        temp_cc[i + taille_c * j] = Kci[i + taille_c * kk] * temp_ic[kk + taille_i * j] + invKii0;

        invKii0 = temp_cc[i + taille_c * j];
      }
    }
  }


  Mlat = (double*) malloc(taille_c * taille_c * sizeof(double));



  for (i = 0; i < taille_c; i++)
    for (j = 0; j < taille_c; j++)
      Mlat[i + taille_c * j] = Kcc[i + taille_c * j] - temp_cc[i + taille_c * j];


  for (i = 0; i < 2 * taille_tt; i++)
    for (j = 0; j < 2 * taille_tt; j++)
      MM[j * 2 * taille_tt + i] = 0.;


  for (i = 0; i < 2 * taille_tt; i++)
    for (j = 0; j < 2 * taille_tt; j++)
      MM[i + 2 * taille_tt * j] = Mlat[i + taille_c * j];




  qi = (double*) malloc(taille_i * sizeof(double));

  for (i = 0; i < taille_i; i++)
  {
    invKii0 = 0.0;

    for (kk = 0; kk < taille_i; kk++)
    {
      qi[i] = invKii[i + taille_i * kk] * Fi[kk] + invKii0;

      invKii0 = qi[i];
    }

  }


  for (i = 0; i < taille_c; i++)
  {
    invKii0 = 0.0;

    for (kk = 0; kk < taille_i; kk++)
    {
      q[i] = Kci[i + taille_c * kk] * qi[kk] + invKii0;

      invKii0 = q[i];
    }

  }


  qbis = (double*) malloc(taille_c * sizeof(double));

  for (i = 0; i < taille_c; i++)
  {
    invKii0 = 0.0;

    for (kk = 0; kk < taille_c; kk++)
    {
      qbis[i] = Mlat[i + taille_c * kk] * Jc[kk] + invKii0;

      invKii0 = qbis[i];
    }

  }


  for (i = 0; i < taille_c; i++)
  {
    q[i] = -q[i] + qbis[i];
  }



  alpha = -1;
  DSCAL(taille_c, alpha, q, incx);



  free(Mlat);
  free(qi);
  free(qbis);
  free(temp_ic);
  free(Kii);
  free(invKii);
  free(Kii2);
  free(Kic);
  free(Fi);

  free(q0);
  free(q1);
  free(q2);
  free(q3);
  free(Kcc);
  free(Kci);
  free(temp_cc);
  free(Jcn);
  free(Jc);
  free(R);
  free(sort);
  free(sort2);
  free(ddl_c);
  free(sort1);
  free(vecF1);
  free(vec_i);
  free(ddl_i);




  /*  return 1;*/
}
