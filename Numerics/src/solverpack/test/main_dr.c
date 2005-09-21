/*!
 *  This main file allows the primal resolution of contact problems with friction:
 *  try (z,w) such that:
 *
 *
 *
 *                    M z  - w = q                 (1)
 *                    0<= zn , 0<= wn , zn.wn= O   (2) (Signorini contact law)
 *                    -wt in diff(PSI)   (zt)      (3) (Coulomb friction law)
 *                                 [-mu*zn, mu*zn]
 *
 *
 *  here M is an n by n  matrix, q an n-dimensional vector, z an n-dimensional*  vector and w an n-dimensional vector.
 *
 *  This system of equations and inequalities is solved thanks to dr subroutines:
 *
 *         dr_nlgs( M,q,n,mu,itermax,tol,z,w,it_end,res,info)
 *        dr_latin( M,q,n,k_latin,mu,itermax,tol,z,w,it_end,res,info)
 *
 *  For more information about the methods see the chapter 4 of the Siconos manual theory.
 *
 *  The subroutine's call is due to the function dr_solver:
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "solverpack.h"

int main(void)
{

  FILE *f1, *f2, *f5, *f6;
  int i, j, nl, nc, nll, info, n = 40, dimM = n;
  double *q, *z, *w, *vec, *a, *b, *zt;
  double **M;
  double qi, Mij;
  char val[14], vall[50];

  method meth_dr;

  M = (double **)malloc(dimM * sizeof(double*));
  for (i = 0; i < n; i++) M[i] = (double*)malloc(dimM * sizeof(double));

  strcpy(meth_dr.dr.name , "Latin");

  meth_dr.dr.itermax = 1000;
  meth_dr.dr.tol = 0.0001;
  meth_dr.dr.k_latin = 0.007; //0.9;//0.00005;



  if ((f1 = fopen("DATA/M_relay_dr.dat", "r")) == NULL)
  {
    perror("fopen 1");
    exit(1);
  }




  vec = (double*)malloc(dimM * dimM * sizeof(double));


  while (!feof(f1))
  {
    fscanf(f1, "%d", &nl);
    fscanf(f1, "%d", &nc);
    fscanf(f1, "%s", val);
    Mij = atof(val);

    /////////////       on met la transpos       ////////////////
    //printf("nl=%i", nl);
    //printf("nc=%i", nc);


    M[nl - 1][nc - 1] = Mij;
    //*(*(M+nc-1)+nl-1)=Mij;

    //////////////         fin transpos         ////////////////////

  }

  //// valeurs du tableau dans vec (compatibilite allocation memoire f90)///
  for (i = 0; i < dimM; i++)
  {
    for (j = 0; j < dimM; j++)
    {
      vec[j * dimM + i] = M[i][j];

    }
  }
  ////////////////////////////////////////////////////////////////////////




  if ((f2 = fopen("DATA/q_relay_dr.dat", "r")) == NULL)
  {
    perror("fopen 2");
    exit(2);
  }


  if ((f5 = fopen("DATA/a_relay_dr.dat", "r")) == NULL)
  {
    perror("fopen 5");
    exit(5);
  }


  if ((f6 = fopen("DATA/b_relay_dr.dat", "r")) == NULL)
  {
    perror("fopen 6");
    exit(6);
  }



  q = malloc(dimM * sizeof(double));
  z = malloc(dimM * sizeof(double));
  zt = malloc(dimM * sizeof(double));
  w = malloc(dimM * sizeof(double));
  a = malloc(dimM * sizeof(double));
  b = malloc(dimM * sizeof(double));





  while (!feof(f2))
  {
    fscanf(f2, "%d", &nll);
    fscanf(f2, "%s", vall);
    qi = atof(vall);
    printf("nll=%i\n", nll);
    q[nll - 1] = qi;
  }


  while (!feof(f5))
  {
    fscanf(f5, "%d", &nll);
    fscanf(f5, "%s", vall);
    qi = atof(vall);
    fprintf(f5, "%d %.14e\n", nll, qi);
    a[nll - 1] = qi;
  }


  while (!feof(f6))
  {
    fscanf(f6, "%d", &nll);
    fscanf(f6, "%s", vall);
    qi = atof(vall);
    fprintf(f6, "%d %.14e\n", nll, qi);
    b[nll - 1] = qi;
  }



  meth_dr.dr.a = (double*)malloc(dimM * sizeof(double));
  meth_dr.dr.b = (double*)malloc(dimM * sizeof(double));

  for (i = 0; i <= n - 1; i++)
  {
    meth_dr.dr.a[i] = a[i];
    meth_dr.dr.b[i] = -b[i];
  }

  printf("\n\n we go in the function \n\n");

  info = dr_solver(vec, q, &n, &meth_dr, zt, w);



  printf("\n\n we go out the function and info is %d\n", info);

  fclose(f2);
  fclose(f1);
  fclose(f5);
  fclose(f6);

  for (i = 0; i < dimM; i++) free(M[i]);
  free(M);


  free(vec);
  free(q);
  free(z);
  free(zt);
  free(w);
  free(a);
  free(b);
  free(meth_dr.dr.a);
  free(meth_dr.dr.b);
  return 0;
}
