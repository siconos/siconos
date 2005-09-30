
/*!\file lcp_cfd.c

   This file allows to give the solution of the contact problem with friction given.

*/

/*!\fn  int lcp_cfd (int *dim_nn, double *ztel, double *wtel, methode *pt, double *K1, double *F1, int * dim_F1, double *J1, int *ddl_i, int * dim_i, int *ddl_c, int *dim_c,  int *ddl_n, int *ddl_tt, int *dim_tt, double *U2,double *F2)

   lcp_cfd subroutine allows to give the solution of the contact problem with friction given.
    \sa cfd_lcp subroutine.

   \param dim_nn On enter a pointer over integers, the size of the vector solution.
   \param ztel On enter a pointer over doubles, the solution given by an LCP ztel.
   \param wtel On enter a pointer over doubles, the solution given by an LCP wtel.
   \param pt On enter a pointer over the union methode.
   \param K1 On enter a pointer over doubles, the rigidity matrix.
   \param F1 On enter a pointer over doubles, the effort vector.
   \param dim_F1 On enter a pointer over integers, the dimension of F1.
   \param J1 On enter a pointer over doubles, the free motion vector.
   \param ddl_i On enter a pointer over integers, the ddl .
   \param dim_i On enter a pointer over integers, the ddl_i size.
   \param ddl_c On enter a pointer over integers, the ddl .
   \param dim_c On enter a pointer over integers, the ddl_c size.
   \param ddl_n On enter a pointer over integers, the ddl.
   \param ddl_tt On enter a pointer over integers, the ddl .
   \param dim_tt On enter a pointer over integers, the ddl_tt size.
   \param U2 On return a pointer over doubles, the solution of the contact friction problem U2(dim_nn).
   \param F2 On return a pointer over doubles, the solution of the contact friction problem F2(dim_nn).

   \author Nineb Sheherazade.
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef MEXFLAG
#include "solverpack.h"
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


  dpotrf_(&uplo, &taille_i, R , &taille_i, &info2);

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


    // return ;

  }



  dpotri_(&uplo, &taille_i, R , &taille_i, &info2);

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


    // return ;

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



  //        return ;
}
