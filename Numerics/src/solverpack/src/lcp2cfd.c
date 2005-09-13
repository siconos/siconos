
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




int lcp_cfd(int *dim_nn, double *ztel, double *wtel, method *pt, double *K1, double *F1, int * dim_F1, double *J1, int *ddl_i, int * dim_i, int *ddl_c, int *dim_c,  int *ddl_n, int *ddl_tt, int *dim_tt, double *U2, double *F2)
{

  FILE *f101, *f202, *f303, *f404, *f000;
  int dim_n, i, j, kk, taille_n = *dim_nn, taille_tt = *dim_tt, taille_i = *dim_i;
  int taille_c = *dim_c, *ddl_nr, *ddl_tr, nc, taille_nr, taille_tr, taille_F1 = *dim_F1;
  double *a1, *a2, *b1, *b2, *Ut, * Un, *Uc, *Fn, *Ft, *Jc, *Ui;
  double *Jcn, *z, *w, *zz, *ww, Mij, invKii0;
  double *temp_i, alpha, beta, r, rrr, r1, r2, r3, rr, *temp_ibis;
  char mot1[10] = "Cfd_latin", mot2[10] = "Gsnl", mot3[10] = "Gcp", mot4[10] = "Lemke", uplo = 'U', trans = 'T', val[14];
  int info1, info2, incx = 1, incy = 1, *vectnt, nl, ind1, ind2;
  int boucle, bouclei, iter1;
  char vall1[14], vall2[14];
  double(*Kc)[taille_F1], (*Kii)[taille_i]; /* Kc[taille_F1][taille_F1], Kii[taille_i][taille_i];*/
  double(*Kii2)[taille_i], **invKii;  /*Kii2[taille_i][taille_i], **invKii;*/
  double(*Kic)[taille_c], (*Kit)[taille_tt]; /*Kic[taille_i][taille_c], Kit[taille_i][taille_tt];*/
  double(*Kin)[taille_tt];
  double *Fi, **R, **RRT;
  double **invR, **invRT, invRTinvR0, invR0, invRT0;


  Kc = malloc(taille_F1 * taille_F1 * sizeof(double));
  Kii = malloc(taille_i * taille_i * sizeof(double));
  Kii2 = malloc(taille_i * taille_i * sizeof(double));
  Kic = malloc(taille_i * taille_c * sizeof(double));
  Kit = malloc(taille_i * taille_tt * sizeof(double));
  Kin = malloc(taille_i * taille_n * sizeof(double));
  Fi = (double *)malloc(taille_i * sizeof(double));


  for (i = 0; i < taille_F1; i++)
    for (j = 0; j < taille_F1; j++)
    {
      Kc[i][j] = K1[i + taille_F1 * j];
    }




  invKii = (double**) malloc(taille_i * sizeof(double));
  invR = (double**) malloc(taille_i * sizeof(double));
  invRT = (double**) malloc(taille_i * sizeof(double));
  R = (double**) malloc(taille_i * sizeof(double));
  RRT = (double**) malloc(taille_i * sizeof(double));
  for (i = 0; i < taille_i; i++)
  {
    invKii[i] = (double*) malloc(taille_i * sizeof(double));
    invR[i] = (double*) malloc(taille_i * sizeof(double));
    invRT[i] = (double*) malloc(taille_i * sizeof(double));
    R[i] = (double*) malloc(taille_i * sizeof(double));
    RRT[i] = (double*) malloc(taille_i * sizeof(double));
  }

  for (i = 0; i < taille_i; i++)
  {
    Fi[i] = F1[ddl_i[i]];
    for (j = 0; j < taille_i; j++)
    {
      Kii[i][j] = 0.;
      Kii2[i][j] = 0.;
    }
  }


  for (i = 0; i < taille_i; i++)
  {
    ind1 = ddl_i[i];
    for (j = 0; j < taille_i; j++)
    {
      ind2 = ddl_i[j];
      Kii[i][j] = Kc[ind1][ind2];/*/16/03/05*/
      Kii2[j][i] = Kii[i][j];
    }
  }

  for (i = 0; i < taille_i; i++)
    for (j = 0; j < taille_n; j++)
      Kin[i][j] = 0.;

  for (i = 0; i < taille_i; i++)
    for (j = 0; j < taille_tt/*n*/; j++)
    {
      Kin[i][j] = Kc[ddl_i[i]][ddl_n[j]];
    }

  for (i = 0; i < taille_i; i++)
    for (j = 0; j < taille_tt; j++)
      Kit[i][j] = 0.;

  for (i = 0; i < taille_i; i++)
    for (j = 0; j < taille_tt; j++)
      Kit[i][j] = Kc[ddl_i[i]][ddl_tt[j]];

  for (i = 0; i < taille_i; i++)
    for (j = 0; j < taille_c; j++)
      Kic[i][j] = 0.;

  for (i = 0; i < taille_i; i++)
    for (j = 0; j < taille_c; j++)
      Kic[i][j] = Kc[ddl_i[i]][ddl_c[j]];


  for (i = 0; i < taille_i; i++)
    for (j = 0; j < taille_i; j++)
    {
      R[i][j] = 0.;
      RRT[i][j] = 0.;
    }


  /*    // !!!!!!!!!!!!!!!!!!!!!Cholesky!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  R[0][0] = sqrt(Kii2[0][0]);


  for (i = 1; i < taille_i; i++)
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

  /*    // !!!!!end of cholesky!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  /*//  !determination of the R tranposeted*/

  for (i = 0; i < taille_i; i++)
    for (j = 0; j < taille_i; j++)
    {
      RRT[i][j] = R[j][i];
      invRT[i][j] = 0.;
      invR[i][j] = 0.;
    }

  /*// !!!!!!!!!inversion of the inf triangular matrix!!!!!!!!!!!!!*/

  for (i = 0; i < taille_i; i++)
  {
    for (j = 0; j < taille_i; j++)
    {
      if (i == j) invR[i][j] = 1 / R[i][j];
      else
      {
        invR0 = 0.;
        for (kk = j; kk <= (i - 1); kk++)
        {
          invR[i][j] = (-1 / R[i][i]) * R[i][kk] * invR[kk][j] + invR0;
          invR0 = invR[i][j];
        }
      }
    }
  }

  /*// !!!!!!!!!!!!!!!!!!!end of inversion!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  /*// !!!!!!!!!!!!!!!!!!!!!!inversion of the sup triangular matrix!!!!!!!*/
  for (i = 0; i < taille_i; i++)
    invRT[i][i] = 1 / RRT[i][i];

  for (i = taille_i - 2; i >= 0; i--)
  {
    for (j = taille_i - 1; j >= 0; j--)
    {
      invRT0 = 0.;
      for (kk = i + 1; kk <= j; kk++)
      {
        invRT[i][j] = (-1 / RRT[i][i]) * RRT[i][kk] * invRT[kk][j] + invRT0;
        invRT0 = invRT[i][j];
      }
    }
  }

  /*// !!!!!!!!!!!!!!!!!!!end of inversion!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/


  for (i = 0; i < taille_i; i++)
    for (j = 0; j < taille_i; j++)
    {
      invKii[i][j] = 0.;
    }


  for (i = 0; i < taille_i; i++)
    for (j = 0; j < taille_i; j++)
    {
      invRTinvR0 = 0.;
      for (kk = 0; kk < taille_i; kk++)
      {
        invKii[i][j] = invRT[i][kk] * invR[kk][j] + invRTinvR0;
        invRTinvR0 = invKii[i][j];
      }
    }


  if (strcmp(pt->cfd.nom_method, mot2) == 0)
  {



    if ((f000 = fopen("resultat_gsnl.dat", "r")) == NULL)
    {
      perror("fopen 1");
      exit(1);
    }

    f101 = fopen("resultat_gsnl_cfd.dat", "w+");

    dim_n = 2 * ((*dim_nn) / 3);
    nc = dim_n / 2;

    boucle = 0;
    while (boucle < 1500)
    {

      for (bouclei = 0; bouclei < *dim_nn; bouclei++)
      {
        fscanf(f000, "%d  %d  %s %s\n", &iter1, &i, &vall1, &vall2);
        ztel[i] = atof(vall1);
        wtel[i] = atof(vall2);
      }
      /*       fprintf(f101,"%d  %d  %14.7e\n",iter1-1,i,z[i]);*/
      //  }



      Jcn = (double*) malloc(nc * sizeof(double));
      z = (double*) malloc(dim_n/*taille_n*/*sizeof(double));
      w = (double*) malloc(dim_n/*taille_n*/*sizeof(double));

      a1 = (double*) malloc(nc * sizeof(double));
      a2 = (double*) malloc(nc * sizeof(double));
      b1 = (double*) malloc(nc * sizeof(double));
      b2 = (double*) malloc(nc * sizeof(double));
      Ut = (double*) malloc(nc * sizeof(double));
      Ft = (double*) malloc(nc * sizeof(double));
      Un = (double*) malloc(nc * sizeof(double));
      Fn = (double*) malloc(nc * sizeof(double));
      temp_i = (double*) malloc(taille_i * sizeof(double));
      temp_ibis = (double*) malloc(taille_i * sizeof(double));

      for (i = 0; i < nc; i++)
      {
        Jcn[i] = J1[ddl_n[i]];
      }


      for (i = 0; i < nc; i++)
      {
        z[i] = ztel[i] - Jcn[i]; /* Un*/
        Un[i] = z[i];
        b1[i] = ztel[nc + i];
        b2[i] = ztel[2 * nc + i];
        z[nc + i] = b2[i] - b1[i]; /* Ut */
        Ut[i] = z[nc + i];
        w[i] = wtel[i]; /* Fn */
        Fn[i] = w[i];
        a1[i] = wtel[nc + i];
        a2[i] = wtel[2 * nc + i];
        w[nc + i] = 0.5 * (a2[i] - a1[i]); /* Ft */
        Ft[i] = w[nc + i];
      }

      vectnt = (int*) malloc(dim_n * sizeof(int));
      ddl_nr = (int*) malloc(nc * sizeof(int));
      ddl_tr = (int*) malloc(nc * sizeof(int));

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

      taille_nr = nc;/*/sizeof(ddl_nr)/sizeof(ddl_nr[0]);*/
      taille_tr = nc;/*/sizeof(ddl_tr)/sizeof(ddl_tr[0]);*/

      zz = (double*) malloc(dim_n * sizeof(double));
      ww = (double*) malloc(dim_n/*taille_n*/*sizeof(double));
      Ui = (double*) malloc(taille_i * sizeof(double));


      for (i = 0; i < taille_nr; i++)
      {
        zz[ddl_nr[i]] = z[i] + Jcn[i]; /*z[i];*/
        zz[ddl_tr[i]] = z[taille_nr + i];
        ww[ddl_nr[i]] = w[i];
        ww[ddl_tr[i]] = w[taille_nr + i];
      }

      /*  f101=fopen("resultat_gsnl.dat","w+");*/

      for (i = 0; i < 2 * taille_nr; i++)
      {
        fprintf(f101, "%d  %d %14.7e\n", iter1, i, zz[i]); //
      }


      dcopy_(&taille_i, Fi, &incx, temp_i, &incy);


      for (i = 0; i < taille_i; i++)
      {
        invKii0 = 0.0;
        for (j = 0; j < taille_tt; j++)
        {
          temp_i[i] = -Kin[i][j] * Un[j] + invKii0;
          invKii0 = temp_i[i];

        }
        temp_i[i] = temp_i[i] + Fi[i];
      }





      for (i = 0; i < taille_i; i++)
      {
        invKii0 = 0.0;
        for (j = 0; j < taille_tt; j++)
        {
          temp_ibis[i] = -Kit[i][j] * Ut[j] + invKii0;
          invKii0 = temp_ibis[i];

        }
        temp_i[i] = temp_ibis[i] + temp_i[i];
      }



      for (i = 0; i < taille_i; i++)
      {
        invKii0 = 0.0;
        for (j = 0; j < taille_i; j++)
        {
          Ui[i] = invKii[i][j] * temp_i[j] + invKii0;
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

      boucle = boucle + 1;
    }

    fclose(f000);
    fclose(f101);


  }
  else if (strcmp(pt->cfd.nom_method, mot3) == 0)
  {

    dim_n = 2 * ((*dim_nn) / 3);
    nc = dim_n / 2;

    Jcn = (double*) malloc(nc * sizeof(double));
    z = (double*) malloc(dim_n/*taille_n*/*sizeof(double));
    w = (double*) malloc(dim_n/*taille_n*/*sizeof(double));
    a1 = (double*) malloc(nc * sizeof(double));
    a2 = (double*) malloc(nc * sizeof(double));
    b1 = (double*) malloc(nc * sizeof(double));
    b2 = (double*) malloc(nc * sizeof(double));
    Ut = (double*) malloc(nc * sizeof(double));
    Ft = (double*) malloc(nc * sizeof(double));
    Un = (double*) malloc(nc * sizeof(double));
    Fn = (double*) malloc(nc * sizeof(double));
    temp_i = (double*) malloc(taille_i * sizeof(double));
    temp_ibis = (double*) malloc(taille_i * sizeof(double));

    for (i = 0; i < nc; i++)
    {
      Jcn[i] = J1[ddl_n[i]];
    }


    for (i = 0; i < nc; i++)
    {
      z[i] = ztel[i] - Jcn[i]; /* Un */
      Un[i] = z[i];
      b1[i] = ztel[nc + i];
      b2[i] = ztel[2 * nc + i];
      z[nc + i] = b2[i] - b1[i]; /* Ut */
      Ut[i] = z[nc + i];
      w[i] = wtel[i]; /* Fn */
      Fn[i] = w[i];
      a1[i] = wtel[nc + i];
      a2[i] = wtel[2 * nc + i];
      w[nc + i] = 0.5 * (a2[i] - a1[i]); /* Ft */
      Ft[i] = w[nc + i];
    }

    vectnt = (int*) malloc(dim_n * sizeof(int));
    ddl_nr = (int*) malloc(nc * sizeof(int));
    ddl_tr = (int*) malloc(nc * sizeof(int));

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

    taille_nr = nc;/*/sizeof(ddl_nr)/sizeof(ddl_nr[0]);*/
    taille_tr = nc;/*/sizeof(ddl_tr)/sizeof(ddl_tr[0]);*/

    zz = (double*) malloc(dim_n * sizeof(double));
    ww = (double*) malloc(dim_n/*taille_n*/*sizeof(double));
    Ui = (double*) malloc(taille_i * sizeof(double));


    for (i = 0; i < taille_nr; i++)
    {
      zz[ddl_nr[i]] = z[i];
      zz[ddl_tr[i]] = z[taille_nr + i];
      ww[ddl_nr[i]] = w[i];
      ww[ddl_tr[i]] = w[taille_nr + i];
    }


    dcopy_(&taille_i, Fi, &incx, temp_i, &incy);

    for (i = 0; i < taille_i; i++)
    {
      invKii0 = 0.0;
      for (j = 0; j < taille_tt; j++)
      {
        temp_i[i] = -Kin[i][j] * Un[j] + invKii0;
        invKii0 = temp_i[i];
      }
      temp_i[i] = temp_i[i] + Fi[i];
    }



    for (i = 0; i < taille_i; i++)
    {
      invKii0 = 0.0;
      for (j = 0; j < taille_tt; j++)
      {
        temp_ibis[i] = -Kit[i][j] * Ut[j] + invKii0;
        invKii0 = temp_ibis[i];
      }
      temp_i[i] = temp_ibis[i] + temp_i[i];
    }



    for (i = 0; i < taille_i; i++)
    {
      invKii0 = 0.0;
      for (j = 0; j < taille_i; j++)
      {
        Ui[i] = invKii[i][j] * temp_i[j] + invKii0;
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


  }
  else if (strcmp(pt->cfd.nom_method, mot4) == 0)
  {

    dim_n = 2 * ((*dim_nn) / 3);
    nc = dim_n / 2;

    Jcn = (double*) malloc(nc * sizeof(double));
    z = (double*) malloc(dim_n/*taille_n*/*sizeof(double));
    w = (double*) malloc(dim_n/*taille_n*/*sizeof(double));
    a1 = (double*) malloc(nc * sizeof(double));
    a2 = (double*) malloc(nc * sizeof(double));
    b1 = (double*) malloc(nc * sizeof(double));
    b2 = (double*) malloc(nc * sizeof(double));
    Ut = (double*) malloc(nc * sizeof(double));
    Ft = (double*) malloc(nc * sizeof(double));
    Un = (double*) malloc(nc * sizeof(double));
    Fn = (double*) malloc(nc * sizeof(double));
    temp_i = (double*) malloc(taille_i * sizeof(double));
    temp_ibis = (double*) malloc(taille_i * sizeof(double));

    for (i = 0; i < nc; i++)
    {
      Jcn[i] = J1[ddl_n[i]];
    }


    for (i = 0; i < nc; i++)
    {
      z[i] = ztel[i] - Jcn[i]; /* Un */
      Un[i] = z[i];
      b1[i] = ztel[nc + i];
      b2[i] = ztel[2 * nc + i];
      z[nc + i] = b2[i] - b1[i]; /* Ut */
      Ut[i] = z[nc + i];
      w[i] = wtel[i]; /* Fn */
      Fn[i] = w[i];
      a1[i] = wtel[nc + i];
      a2[i] = wtel[2 * nc + i];
      w[nc + i] = 0.5 * (a2[i] - a1[i]); /* Ft */
      Ft[i] = w[nc + i];
    }

    vectnt = (int*) malloc(dim_n * sizeof(int));
    ddl_nr = (int*) malloc(nc * sizeof(int));
    ddl_tr = (int*) malloc(nc * sizeof(int));

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

    taille_nr = nc;/*/sizeof(ddl_nr)/sizeof(ddl_nr[0]);*/
    taille_tr = nc;/*/sizeof(ddl_tr)/sizeof(ddl_tr[0]);*/

    zz = (double*) malloc(dim_n * sizeof(double));
    ww = (double*) malloc(dim_n/*taille_n*/*sizeof(double));
    Ui = (double*) malloc(taille_i * sizeof(double));


    for (i = 0; i < taille_nr; i++)
    {
      zz[ddl_nr[i]] = z[i] + Jcn[i]; /*z(i)*/
      zz[ddl_tr[i]] = z[taille_nr + i];
      ww[ddl_nr[i]] = w[i];
      ww[ddl_tr[i]] = w[taille_nr + i];
    }


    f404 = fopen("resultat_lemke_cfd.dat", "w+");
    for (i = 0; i < 2 * taille_nr; i++)
    {
      fprintf(f404, "%d  %d %14.7e\n", i, i, zz[i]);
    }/**/
    fclose(f404);


    dcopy_(&taille_i, Fi, &incx, temp_i, &incy);

    for (i = 0; i < taille_i; i++)
    {
      invKii0 = 0.0;
      for (j = 0; j < taille_tt; j++)
      {
        temp_i[i] = -Kin[i][j] * Un[j] + invKii0;
        invKii0 = temp_i[i];
      }
      temp_i[i] = temp_i[i] + Fi[i];
    }



    for (i = 0; i < taille_i; i++)
    {
      invKii0 = 0.0;
      for (j = 0; j < taille_tt; j++)
      {
        temp_ibis[i] = -Kit[i][j] * Ut[j] + invKii0;
        invKii0 = temp_ibis[i];
      }
      temp_i[i] = temp_ibis[i] + temp_i[i];
    }



    for (i = 0; i < taille_i; i++)
    {
      invKii0 = 0.0;
      for (j = 0; j < taille_i; j++)
      {
        Ui[i] = invKii[i][j] * temp_i[j] + invKii0;
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


  }

  else if (strcmp(pt->cfd.nom_method, mot1) == 0)
  {

    Jc = (double*) malloc(taille_c * sizeof(double));
    Uc = (double*) malloc(taille_c * sizeof(double));
    Ui = (double*) malloc(taille_i * sizeof(double));
    temp_i = (double*) malloc(taille_i * sizeof(double));


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
        temp_i[i] = Kic[i][j] * Uc[j] + invKii0;
        invKii0 = temp_i[i];
      }
      temp_i[i] = Fi[i] - temp_i[i];
    }

    for (i = 0; i < taille_i; i++)
    {
      invKii0 = 0.0;
      for (j = 0; j < taille_i; j++)
      {
        Ui[i] = invKii[i][j] * temp_i[j] + invKii0;
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


  }
  else printf("Warning : Unknown solving method : %s\n", pt->cfd.nom_method);

  for (i = 0; i < taille_i; i++)
  {
    free(invKii[i]);
    free(invR[i]);
    free(invRT[i]);
    free(R[i]);
    free(RRT[i]);

  }
  free(invKii);
  free(invR);
  free(invRT);
  free(R);
  free(RRT);

  free(Kc);
  free(Kii);
  free(Kii2);
  free(Kic);
  free(Kit);
  free(Kin);
  free(Fi);


  return 3;
}
