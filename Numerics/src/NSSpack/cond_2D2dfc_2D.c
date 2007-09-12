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

/*!\file cond_2D2dfc_2D.c

   This file allows to give the solution of the 2D contact problem with friction given.\n

  \fn  void cond_2D2dfc_2D (int *dim_F1, double *ztel, double *wtel, double *K1, double *F1, double *J1, int *ddl_n, int *ddl_tt, int * dim_tt, int *ddl_d, int *dim_d, double *U2,double *F2)

   cond_2D2dfc_2D.c subroutine allows to give the solution of the 2D contact problem with friction given.
    \sa dfc_2D2cond_2D subroutine.\n




   \param dim_F1    On enter a pointer over integers, the dimension of the DFC_2D problem,

   \param ztel      On enter a pointer over doubles, the solution given by a dfc_2D solver.

   \param wtel      On enter a pointer over doubles, the solution given by a dfc_2D solver.

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


void cond_2D2dfc_2D(int *dim_F1, double *ztel, double *wtel, double *K1, double *F1, double *J1, int *ddl_n, int *ddl_tt, int * dim_tt, int *ddl_d, int *dim_d, double *U2, double *F2)
{




  int      info2, ind1, ind2;
  int      i, j, taille_n = *dim_tt, taille_i;
  int      taille_c = 2 * taille_n, taille_F1 = *dim_F1;
  int      taille_d = *dim_d, taille_sort;
  int      *ddl_c, *sort1, *sort, *sort2, *vecF1, *vec_i, *ddl_i;




  double   *Uc, *Jc, *Ui;
  double   invKii0;
  double   *temp_i;


  double   *Kii;
  double   *Kii2, *invKii;
  double   *Kic;

  double   *Fi, *R;

  char     uplo = 'U';



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


  Kii     = (double *) malloc(taille_i * taille_i * sizeof(double));
  Kii2    = (double *) malloc(taille_i * taille_i * sizeof(double));
  Kic     = (double *) malloc(taille_i * taille_c * sizeof(double));

  Fi      = (double *) malloc(taille_i * sizeof(double));
  invKii  = (double *) malloc(taille_i * taille_i * sizeof(double));
  R       = (double *) malloc(taille_i * taille_i * sizeof(double));
  Jc      = (double *) malloc(taille_c * sizeof(double));
  Uc      = (double *) malloc(taille_c * sizeof(double));

  Ui      = (double *) malloc(taille_i * sizeof(double));
  temp_i  = (double *) malloc(taille_i * sizeof(double));



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
    for (j = 0; j < taille_c; j++)
      Kic[i + taille_i * j] = 0.;


  for (i = 0; i < taille_i; i++)
    for (j = 0; j < taille_c; j++)
      Kic[i + taille_i * j] = K1[ddl_i[i] + taille_F1 * ddl_c[j]];


  for (i = 0; i < taille_i; i++)
    for (j = 0; j < taille_i; j++)
    {
      R[i + taille_i * j] = 0.;

    }




  for (i = 0; i < taille_i; i++)
    for (j = 0; j < taille_i; j++)
    {
      R[i + taille_i * j] = Kii2[i + taille_i * j];
    }





  /*                   Cholesky                 */


  dpotrf_(&uplo, (integer *)&taille_i, R , (integer *)&taille_i, (integer *)&info2);

  if (info2 != 0)
  {

    free(Jc);
    free(Uc);
    free(Ui);
    free(temp_i);
    free(invKii);

    free(R);
    free(Kii);
    free(Kii2);
    free(Kic);
    free(Fi);

    free(vec_i);
    free(ddl_c);
    free(sort1);
    free(sort2);
    free(sort);

    free(vecF1);
    free(ddl_i);


    return ;

  }



  dpotri_(&uplo, (integer*)&taille_i, R , (integer*)&taille_i, (integer*)&info2);

  if (info2 != 0)
  {

    free(Jc);
    free(Uc);
    free(Ui);
    free(temp_i);
    free(invKii);

    free(R);
    free(Kii);
    free(Kii2);
    free(Kic);
    free(Fi);


    free(vec_i);
    free(ddl_c);
    free(sort1);
    free(sort2);
    free(sort);

    free(vecF1);
    free(ddl_i);


    return ;

  }


  for (i = 0; i < taille_i; i++)
    for (j = 0; j < taille_i; j++)
    {
      invKii[j + taille_i * i] =  R[j + taille_i * i];

      invKii[i + taille_i * j] =  R[j + taille_i * i];


    }




  for (i = 0; i < taille_c; i++)
  {
    Jc[i] = J1[ddl_c[i]];
    Uc[i] = ztel[i] - Jc[i];
  }



  for (i = 0; i < taille_i; i++)
  {
    invKii0 = 0.0;
    for (j = 0; j < taille_c; j++)
    {
      temp_i[i] = Kic[i + taille_i * j] * Uc[j] + invKii0;
      invKii0   = temp_i[i];
    }
    temp_i[i] = Fi[i] - temp_i[i];
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

  for (i = 0; i < taille_c; i++)
  {
    U2[ddl_c[i]] = Uc[i];
  }

  for (i = 0; i < taille_i; i++)
  {
    F2[ddl_i[i]] = Fi[i];
  }

  for (i = 0; i < taille_c; i++)
  {
    F2[ddl_c[i]] = wtel[i];
  }



  free(Jc);
  free(Uc);
  free(Ui);
  free(temp_i);
  free(invKii);

  free(R);
  free(Kii);
  free(Kii2);
  free(Kic);
  free(Fi);

  free(vec_i);
  free(ddl_c);
  free(sort1);
  free(sort2);
  free(sort);

  free(vecF1);
  free(ddl_i);



  /*        return ;   */
}
