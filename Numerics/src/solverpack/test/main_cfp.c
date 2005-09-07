///////////////////////////////////////////////////////////////////////////
//  This main file allows the primal resolution of contact problems with friction:
//  try (z,w) such that:
//
//
//                    M z  - w = q                 (1)
//                    0<= zn , 0<= wn , zn.wn= O   (2) (Signorini contact law)
//                    -wt in diff(PSI)   (zt)      (3) (Coulomb friction law)
//                                 [-mu*zn, mu*zn]
//
//
//  here M is an n by n  matrix, q an n-dimensional vector, z an n-dimensional
//  vector and w an n-dimensional vector.
//
//  This system of equations and inequalities is solved thanks to cfp_subroutine:
//        cfp_gsnl (M,q,n,mu,itermax,tol,z,w,it_end,res,info)
//        cfp_gcp (M,q,n,mu,itermax,tol,z,w,it_end,res,info)
//        cfp_latin (M,q,n,k_latin,mu,itermax,tol,z,w,it_end,res,info)
//
//  where _ itermax is the maximum iterations required, it's an integer
//        _ mu is the friction coefficient, it's a float
//        _ res is the residue, it's a float
//        _ it_end is the number of iterations carried out, it's an integer
//        _ tol is the tolerance value desired, it's a float
//        _ k_latin is the parameter research of the latin, it's a float
//        _ z and w are the solutions of the problem
//        _ info shows the termination reason,0 is successful otherwise 1, it's an integer.
//
//
//    For more information about the methods see the chapter 4 of the Siconos manual theory.
//
//
//
//  The subroutine's call is due to the function solve_cfp:
//
//  int solve_cfp (double (*M)[maxcols],double *q,int n,methode *pt, double *z,double *w)
//
//  where M is an n by n matrix, q an n-dimensional vector, n is the row dimension
//  of M, and pt a pointer other a structure ( methode). z and w are n-dimensional
//  vectors solution.
//  methode is a variable with a structure type; this structure gives to the function
// solve_pcf, the name and the parameters (itermax, tol, k_latin,...) of the method we want to use.
//  This function return an interger:  0 successful return otherwise 1.
//
//
///////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "solverpack.h"



main()

{
  FILE *f1, *f2, *f3, *f4, *f5, *f6;
  int i, j, nl, nc, nll, it, n = 152, dimM = n;
  double *q, *z, *w, *vec;
  int info;
  double(*M)[n];
  double qi, Mij;
  char val[14], vall[14];
  static methode_cfp meth_cfp  = {"Gsnl", 109, -0.0001, 0.3, 0.7};

  if ((f1 = fopen("DATA/MM_gran_mu12.dat", "r")) == NULL)
  {
    perror("fopen 1");
    exit(1);
  }



  M = malloc(dimM * dimM * sizeof(double));
  vec = (double*)malloc(dimM * dimM * sizeof(double));


  while (!feof(f1))
  {
    fscanf(f1, "%d", &nl);
    fscanf(f1, "%d", &nc);
    fscanf(f1, "%s", val);
    Mij = atof(val);
    /////////////       on met la transpos       ////////////////
    *(*(M + nc - 1) + nl - 1) = Mij;
    //////////////         fin transpos         ////////////////////
  }

  //// valeurs du tableau dans vec (compatibilite allocation memoire f90)///
  for (i = 0; i < dimM; i++)
    for (j = 0; j < dimM; j++)
      vec[j * dimM + i] = M[i][j];
  ////////////////////////////////////////////////////////////////////////




  if ((f2 = fopen("DATA/qq_gran_mu12.dat", "r")) == NULL)
  {
    perror("fopen 2");
    exit(2);
  }



  q = malloc(dimM * sizeof(double));
  z = malloc(dimM * sizeof(double));
  w = malloc(dimM * sizeof(double));




  while (!feof(f2))
  {
    fscanf(f2, "%d", &nll);
    fscanf(f2, "%s", vall);
    qi = atof(vall);
    *(q + nll - 1) = qi;
  }


  printf("\n\n we go in the function Gsnl \n\n");

  info = solve_cfp(vec, q, &n, &meth_cfp, z, w);

  for (i = 0; i < n; i++)
    printf("z %g w %g\n", z[i], w[i]);


  printf("\n\n we go out the function and info is %d\n", info);

  fclose(f2);
  fclose(f1);

  free(q);
  free(z);
  free(w);
  free(M);
  free(vec);



  ///////////////////////////////////////
  // second test ////////////////////////
  ///////////////////////////////////////
  static methode_cfp meth_cfp2  = {"Latin", 1000, 0.0001, 0.3, 35};

  if ((f1 = fopen("DATA/MM_gran_mu12.dat", "r")) == NULL)
  {
    perror("fopen 1");
    exit(1);
  }



  M = malloc(dimM * dimM * sizeof(double));
  vec = (double*)malloc(dimM * dimM * sizeof(double));


  while (!feof(f1))
  {
    fscanf(f1, "%d", &nl);
    fscanf(f1, "%d", &nc);
    fscanf(f1, "%s", val);
    Mij = atof(val);
    /////////////       on met la transpos       ////////////////
    *(*(M + nc - 1) + nl - 1) = Mij;
    //////////////         fin transpos         ////////////////////
  }

  //// valeurs du tableau dans vec (compatibilite allocation memoire f90)///
  for (i = 0; i < dimM; i++)
    for (j = 0; j < dimM; j++)
      vec[j * dimM + i] = M[i][j];
  ////////////////////////////////////////////////////////////////////////




  if ((f2 = fopen("DATA/qq_gran_mu12.dat", "r")) == NULL)
  {
    perror("fopen 2");
    exit(2);
  }



  q = malloc(dimM * sizeof(double));
  z = malloc(dimM * sizeof(double));
  w = malloc(dimM * sizeof(double));




  while (!feof(f2))
  {
    fscanf(f2, "%d", &nll);
    fscanf(f2, "%s", vall);
    qi = atof(vall);
    *(q + nll - 1) = qi;
  }


  printf("\n\n we go in the function Latin\n\n");

  info = solve_cfp(vec, q, &n, &meth_cfp2, z, w);

  for (i = 0; i < n; i++)
    printf("z %g w %g\n", z[i], w[i]);


  printf("\n\n we go out the function and info is %d\n", info);

  fclose(f2);
  fclose(f1);

  free(q);
  free(z);
  free(w);
  free(M);
  free(vec);



  ///////////////////////////////////////
  //   third test ////////////////////////
  ///////////////////////////////////////
  static methode_cfp meth_cfp3  = {"Gcp", 1000, 0.0001, 0.3, 0.7};

  if ((f1 = fopen("DATA/MM_gran_mu12.dat", "r")) == NULL)
  {
    perror("fopen 1");
    exit(1);
  }



  M = malloc(dimM * dimM * sizeof(double));
  vec = (double*)malloc(dimM * dimM * sizeof(double));


  while (!feof(f1))
  {
    fscanf(f1, "%d", &nl);
    fscanf(f1, "%d", &nc);
    fscanf(f1, "%s", val);
    Mij = atof(val);
    /////////////       on met la transpos       ////////////////
    *(*(M + nc - 1) + nl - 1) = Mij;
    //////////////         fin transpos         ////////////////////
  }

  //// valeurs du tableau dans vec (compatibilite allocation memoire f90)///
  for (i = 0; i < dimM; i++)
    for (j = 0; j < dimM; j++)
      vec[j * dimM + i] = M[i][j];
  ////////////////////////////////////////////////////////////////////////




  if ((f2 = fopen("DATA/qq_gran_mu12.dat", "r")) == NULL)
  {
    perror("fopen 2");
    exit(2);
  }



  q = malloc(dimM * sizeof(double));
  z = malloc(dimM * sizeof(double));
  w = malloc(dimM * sizeof(double));




  while (!feof(f2))
  {
    fscanf(f2, "%d", &nll);
    fscanf(f2, "%s", vall);
    qi = atof(vall);
    *(q + nll - 1) = qi;
  }


  printf("\n\n we go in the function Gcp \n\n");

  info = solve_cfp(vec, q, &n, &meth_cfp3, z, w);

  for (i = 0; i < n; i++)
    printf("z %g w %g\n", z[i], w[i]);


  printf("\n\n we go out the function and info is %d\n", info);

  fclose(f2);
  fclose(f1);

  free(q);
  free(z);
  free(w);
  free(M);
  free(vec);
}


