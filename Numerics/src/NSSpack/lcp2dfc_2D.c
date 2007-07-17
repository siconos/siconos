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

/*!\file lcp2dfc_2D.c

   This file allows to give the solution of the 2D contact problem with friction given. \n

  \fn void lcp2dfc_2D (int *dim_F1, double *ztel, double *wtel, double *K1, double *F1, double *J1, int *ddl_n, int *ddl_tt, int * dim_tt, int *ddl_d, int *dim_d, double *U2,double *F2)

   lcp2dfc_2D.c subroutine allows to give the solution of the 2D contact problem with friction given.
    \sa dfc_2D2lcp subroutine.\n


   \param dim_F1    On enter a pointer over integers, the dimension of the DFC_2D problem,

   \param ztel      On enter a pointer over doubles, the solution given by a LCP solver.

   \param wtel      On enter a pointer over doubles, the solution given by a LCP solver.

   \param K1        On enter a pointer over doubles containing the components of the
                     rigidity matrix with a fortran90 storage,

   \param F1        On enter a pointer over doubles containing the right hand side,

   \param J1        On enter a pointer over doubles, gap in normal contact direction.

   \param ddl_n     On enter a pointer over integers , the contact in normal direction dof
                     (not prescribed),

   \param ddl_tt    On enter a pointer over integers, the contact in tangential direction dof
                     (not prescribed)


   \param dim_tt    On enter a pointer over integers, the dimension of the vector ddl_tt.

   \param ddl_d     On enter a pointer over integers, the prescribed dof,

   \param dim_d     On enter a pointer over integers, the dimension of the vector ddl_d,


   \n\n

   \param U2        On return a pointer over doubles, the solution of the contact friction problem U2(dim_F1).
   \param F2        On return a pointer over doubles, the solution of the contact friction problem F2(dim_F1).





   \author Nineb Sheherazade.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef MEXFLAG
#include "NSSpack.h"
#endif
#include "blaslapack.h"



void sortsn_(int *ddl_i, int *sort, int *n);

void lcp2dfc_2D(int *dim_F1, double *ztel, double *wtel, double *K1, double *F1, double *J1, int *ddl_n, int *ddl_tt, int * dim_tt, int *ddl_d, int *dim_d, double *U2, double *F2)
{





  int      info2,  *vectnt, ind1, ind2;
  integer incx = 1, incy = 1;
  int      dim_n, i, j, taille_n = *dim_tt, taille_tt = *dim_tt, taille_i;
  int      taille_c, *ddl_nr, *ddl_tr, nc, taille_nr, taille_tr, taille_F1 = *dim_F1;
  int      taille_d = *dim_d, taille_sort;
  int      *ddl_c, *sort1, *sort, *sort2, *vecF1, *vec_i, *ddl_i;


  double   *a1, *a2, *b1, *b2, *Ut, * Un, *Fn, *Ft, *Ui;
  double   *Jcn, *z, *w, *zz, *ww, invKii0;
  double   *temp_i, *temp_ibis;


  double   *Kii;
  double   *Kii2, *invKii;
  double   *Kit;
  double   *Kin;

  double   *Fi, *R;

  char     uplo = 'U';



  taille_c = 2 * taille_n;

  ddl_c     = (int*) malloc(taille_c * sizeof(int));
  sort1     = (int*) malloc(taille_c * sizeof(int));
  sort      = (int*) malloc((taille_c + taille_d) * sizeof(int));
  sort2     = (int*) malloc((taille_c + taille_d) * sizeof(int));
  vecF1     = (int*) malloc(taille_F1 * sizeof(int));
  vec_i     = (int*) malloc(taille_F1 * sizeof(int));




  for (i = 0; i < taille_n; i++)
  {

    ddl_c[i]          = ddl_n[i];
    ddl_c[i + taille_n] = ddl_tt[i];

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





  dim_n = 2 * taille_n;
  nc    = taille_n;


  Kii       = (double *) malloc(taille_i * taille_i * sizeof(double));
  Kii2      = (double *) malloc(taille_i * taille_i * sizeof(double));
  Kit       = (double *) malloc(taille_i * taille_tt * sizeof(double));
  Kin       = (double *) malloc(taille_i * taille_n * sizeof(double));
  Fi        = (double *) malloc(taille_i * sizeof(double));
  invKii    = (double *) malloc(taille_i * taille_i * sizeof(double));
  R         = (double *) malloc(taille_i * taille_i * sizeof(double));
  Jcn       = (double *) malloc(nc * sizeof(double));
  z         = (double *) malloc(dim_n * sizeof(double));
  w         = (double *) malloc(dim_n * sizeof(double));
  a1        = (double *) malloc(nc * sizeof(double));
  a2        = (double *) malloc(nc * sizeof(double));
  b1        = (double *) malloc(nc * sizeof(double));
  b2        = (double *) malloc(nc * sizeof(double));
  Ut        = (double *) malloc(nc * sizeof(double));
  Ft        = (double *) malloc(nc * sizeof(double));
  Un        = (double *) malloc(nc * sizeof(double));
  Fn        = (double *) malloc(nc * sizeof(double));
  temp_i    = (double *) malloc(taille_i * sizeof(double));
  temp_ibis = (double *) malloc(taille_i * sizeof(double));
  vectnt    = (int *)    malloc(dim_n * sizeof(int));
  ddl_nr    = (int *)    malloc(nc * sizeof(int));
  ddl_tr    = (int *)    malloc(nc * sizeof(int));
  zz        = (double *) malloc(dim_n * sizeof(double));
  ww        = (double *) malloc(dim_n * sizeof(double));
  Ui        = (double *) malloc(taille_i * sizeof(double));








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

  for (i = 0; i < taille_i; i++)
    for (j = 0; j < taille_n; j++)
      Kin[i + taille_i * j] = 0.;

  for (i = 0; i < taille_i; i++)
    for (j = 0; j < taille_tt; j++)
    {
      Kin[i + taille_i * j] = K1[ddl_i[i] + taille_F1 * ddl_n[j]];
    }

  for (i = 0; i < taille_i; i++)
    for (j = 0; j < taille_tt; j++)
      Kit[i + taille_i * j] = 0.;

  for (i = 0; i < taille_i; i++)
    for (j = 0; j < taille_tt; j++)
      Kit[i + taille_i * j] = K1[ddl_i[i] + taille_F1 * ddl_tt[j]];



  for (i = 0; i < taille_i; i++)
    for (j = 0; j < taille_i; j++)
    {
      R[i + taille_i * j] = Kii2[i + taille_i * j];
    }





  /*                   Cholesky            */


  dpotrf_(&uplo, (integer *)&taille_i, R , (integer *)&taille_i, (integer *)&info2);


  if (info2 != 0)
  {

    free(Jcn);
    free(z);
    free(w);
    free(zz);
    free(ww);
    free(a1);
    free(a2);
    free(b1);
    free(b2);
    free(Ut);
    free(Ft);
    free(Un);
    free(Fn);
    free(Ui);
    free(temp_i);
    free(temp_ibis);
    free(vectnt);
    free(ddl_nr);
    free(ddl_tr);
    free(invKii);
    free(R);
    free(Kii);
    free(Kii2);
    free(Kit);
    free(Kin);
    free(Fi);

    free(vec_i);
    free(ddl_c);
    free(sort1);
    free(sort2);
    free(sort);
    free(vecF1);
    free(ddl_i);

    /*return ;*/
  }



  dpotri_(&uplo, &taille_i, R , &taille_i, &info2);

  if (info2 != 0)
  {

    free(Jcn);
    free(z);
    free(w);
    free(zz);
    free(ww);
    free(a1);
    free(a2);
    free(b1);
    free(b2);
    free(Ut);
    free(Ft);
    free(Un);
    free(Fn);
    free(Ui);
    free(temp_i);
    free(temp_ibis);
    free(vectnt);
    free(ddl_nr);
    free(ddl_tr);
    free(invKii);
    free(R);
    free(Kii);
    free(Kii2);
    free(Kit);
    free(Kin);
    free(Fi);


    free(vec_i);
    free(ddl_c);
    free(sort1);
    free(sort2);
    free(sort);
    free(vecF1);
    free(ddl_i);

    /*return ;*/

  }



  for (i = 0; i < taille_i; i++)
    for (j = 0; j < taille_i; j++)
    {
      invKii[j + taille_i * i] =  R[j + taille_i * i];

      invKii[i + taille_i * j] =  R[j + taille_i * i];


    }






  for (i = 0; i < nc; i++)
  {
    Jcn[i] = J1[ddl_n[i]];
  }




  for (i = 0; i < nc; i++)
  {
    z[i]     = ztel[i] - Jcn[i];      /* Un*/
    Un[i]    = z[i];
    b1[i]    = ztel[nc + i];
    b2[i]    = ztel[2 * nc + i];
    z[nc + i]  = b2[i] - b1[i];       /* Ut */
    Ut[i]    = z[nc + i];
    w[i]     = wtel[i];               /* Fn */
    Fn[i]    = w[i];
    a1[i]    = wtel[nc + i];
    a2[i]    = wtel[2 * nc + i];
    w[nc + i]  = 0.5 * (a2[i] - a1[i]); /* Ft */
    Ft[i]    = w[nc + i];
  }





  for (i = 0; i < dim_n; i++)
  {
    vectnt[i] = i + 1;
  }

  for (i = 0; i < nc; i++)
  {
    ddl_nr[i] = vectnt[2 * i];

    if (i != 0) ddl_tr[i] = vectnt[2 * i - 1];

    else ddl_tr[i] = 0;

  }

  taille_nr = nc;
  taille_tr = nc;


  for (i = 0; i < taille_nr; i++)
  {
    zz[ddl_nr[i]] = z[i] + Jcn[i];
    zz[ddl_tr[i]] = z[taille_nr + i];
    ww[ddl_nr[i]] = w[i];
    ww[ddl_tr[i]] = w[taille_nr + i];
  }





  dcopy_((integer *)&taille_i, Fi, &incx, temp_i, &incy);


  for (i = 0; i < taille_i; i++)
  {
    invKii0 = 0.0;

    for (j = 0; j < taille_tt; j++)
    {
      temp_i[i] = -Kin[i + taille_i * j] * Un[j] + invKii0;
      invKii0   = temp_i[i];

    }

    temp_i[i] = temp_i[i] + Fi[i];
  }





  for (i = 0; i < taille_i; i++)
  {
    invKii0 = 0.0;

    for (j = 0; j < taille_tt; j++)
    {
      temp_ibis[i] = -Kit[i + taille_i * j] * Ut[j] + invKii0;
      invKii0      = temp_ibis[i];

    }

    temp_i[i] = temp_ibis[i] + temp_i[i];
  }



  for (i = 0; i < taille_i; i++)
  {
    invKii0 = 0.0;

    for (j = 0; j < taille_i; j++)
    {
      Ui[i]   = invKii[i + taille_i * j] * temp_i[j] + invKii0;
      invKii0 = Ui[i];
    }

  }


  for (i = 0; i < taille_i; i++)
  {
    U2[ddl_i[i]] = Ui[i];
  }

  for (i = 0; i < taille_tt; i++)
  {
    U2[ddl_n[i]] = Un[i];
  }

  for (i = 0; i < taille_tt; i++)
  {
    U2[ddl_tt[i]] = Ut[i];
  }

  for (i = 0; i < taille_i; i++)
  {
    F2[ddl_i[i]] = Fi[i];
  }

  for (i = 0; i < taille_tt; i++)
  {
    F2[ddl_n[i]] = Fn[i];
  }

  for (i = 0; i < taille_tt; i++)
  {
    F2[ddl_tt[i]] = Ft[i];
  }




  free(Jcn);
  free(z);
  free(w);
  free(zz);
  free(ww);
  free(a1);
  free(a2);
  free(b1);
  free(b2);
  free(Ut);
  free(Ft);
  free(Un);
  free(Fn);
  free(Ui);
  free(temp_i);
  free(temp_ibis);
  free(vectnt);
  free(ddl_nr);
  free(ddl_tr);
  free(invKii);
  free(R);
  free(Kii);
  free(Kii2);
  free(Kit);
  free(Kin);
  free(Fi);

  free(vec_i);
  free(ddl_c);
  free(sort1);
  free(sort2);
  free(sort);
  free(vecF1);
  free(ddl_i);



  /*  return ;   */
}
