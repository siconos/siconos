///////////////////////////////////////////////////////////////////////////
//  This main file allows the dual resolution of contact problems with
//  friction:
//  try (z,w) such that:
//
//
//                    M z  - w = q                 (1)
//                    0<= zn , 0<= wn , zn.wn= O   (2)
//                    -zt in diff(PSI)   (wt)      (3)
//                                 [-mu*wn, mu*wn]
//
//
//  here M is an n by n  matrix, q an n-dimensional vector, z an n-dimensional
//  vector and w an n-dimensional vector.
//
//  This system of equations and inequalities is solved thanks to cfd_subroutine:
//        cfd_latin (M,q,n,k_latin,mu,itermax,tol,z,w,it_end,res,info)
//        ....
//        ....
//  or thanks to LCP (Linear Complementary Problem) routines after a new
//  formulation of this problem in the shape of LCP due to the cfd_lcp and lcp_cfd routines:
//
//       cfd_lcp (nc,mu,M,q,Mtel,qtel)
//       lcp_cfd (nc,ztel,wtel,z,w)
//
//
//
//  where _ itermax is the maximum iterations required, it's an integer
//        _ mu is the friction coefficient, it's a float
//        _ res is the residue, it's a float
//        _ it_end is the number of iterations carried out, it's an integer
//        _ tol is the tolerance value desired, it's a float
//        _ k_latin is the parameter research of the latin, it's a double
//        _ z and w are the solutions of the problem.
//        _ nc: integer, nc=n/2 is the number of contacts.
//        _ Mtel(3*nc,3*nc): on return, real matrix after the new formulation.
//        _ qtel(3*nc): on return, real vector after the new formulation.
//        _ ztel(3*nc): real vector, given by a lcp routine.
//        _ wtel(3*nc): real vector, given by a lcp routine.
//        _ z(2*nc): on return, real vector, solution of the problem.
//        _ w(2*nc): on return, real vector, solution of the problem.
//        _ info : shows the termination reason, 0 is successful otherwise 1, it's an ineger.
//
//    For more information about the method(s) see the chapter 4 of the
//    Siconos manual theory.
//
//
//
//  The subroutine's call is due to the function solve_cfd:
//
//   int solve_cfd (double (*M)[maxcols],double *q,int n,methode *pt,double z[],double w[])
//
//  where M is an n by n matrix, q an n-dimensional vector, n is the row
//  dimension of M, and pt a pointer other a structure ( methode).
//  methode is a variable with a structure type; this structure gives to the
//  function solve_lcp, the name and the parameters (itermax, tol, k_latin,..)
//  of the method we want to use.
//  This function return an interger:  0 successful return otherwise 1.
//
//
///////////////////////////////////////////////////////////////////////////////


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "solverpack.h"



double ddot_(int *, double [], int *, double [], int*);
void dpotri_(char *, int *, double * , int *, int *);



main()

{
  FILE *f1, *f2, *f3, *f4, *f5, *f6, *f7;
  int i, j, nl, nc, nll, it, n = 596, dimM = n, dim_i = 529, dim_n = 25, dim_tt = 25, dim_c = 50, info = -7;
  int *ddl_i, *ddl_tt, *ddl_n, *ddl_c, m, *dim_nn, taille_i, taille_F1;
  double *q, *z, *w, *vec, *K1, *J1, *F1, *U2, *F2;
  double(*M)[n], *MM;
  double qi, Mij, mumu;
  char val[14], vall[14];
  char nom[64] = "Gsnl"; //Cfd_latin";//Gsnl";//Lemke";
  methode *pt;
  // static methode_cfd meth_cfd = {"Gnsl",1000,0.0000001,0.7,0.6,dim_i,dim_tt,dim_c,};

  methode meth_cfd;

  // meth_cfd.cfd.nom_method = nom;
  meth_cfd.cfd.itermax = 1000;
  meth_cfd.cfd.tol = 0.000001;
  meth_cfd.cfd.mu = 0.5;
  meth_cfd.cfd.k_latin = 0.6;
  meth_cfd.cfd.dim_i = dim_i;
  meth_cfd.cfd.dim_n = dim_n;
  meth_cfd.cfd.dim_tt = dim_tt;
  meth_cfd.cfd.dim_c = dim_c;



  if ((f1 = fopen("K1_mu.dat", "r")) == NULL)
  {
    perror("fopen 1");
    exit(1);
  }


  M = malloc(dimM * dimM * sizeof(double));
  vec = (double*)malloc(dimM * dimM * sizeof(double));
  K1 = (double*)malloc(dimM * dimM * sizeof(double));

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
  {
    for (j = 0; j < dimM; j++)
    {
      vec[j * dimM + i] = M[i][j];
      K1[j * dimM + i] = M[i][j];
    }
  }
  //// valeurs du tableau dans vec (compatibilite allocation memoire f90)///



  if ((f2 = fopen("J1_mu.dat", "r")) == NULL)
  {
    perror("fopen 2");
    exit(2);
  }




  meth_cfd.cfd.J1 = (double *) malloc(dimM * sizeof(double));
  F1 = (double *) malloc(dimM * sizeof(double));
  meth_cfd.cfd.ddl_i = (int*)malloc(dim_i * sizeof(int));
  meth_cfd.cfd.ddl_n = (int*)malloc(dim_n * sizeof(int));
  meth_cfd.cfd.ddl_tt = (int*)malloc(dim_tt * sizeof(int));
  meth_cfd.cfd.ddl_c = (int*)malloc(dim_c * sizeof(int));



  while (!feof(f2))
  {
    fscanf(f2, "%d", &nll);
    fscanf(f2, "%s", vall);
    qi = atof(vall);
    *(meth_cfd.cfd.J1 + nll - 1) = qi;
  }



  if ((f3 = fopen("F1_mu.dat", "r")) == NULL)
  {
    perror("fopen 3");
    exit(3);
  }


  while (!feof(f3))
  {
    fscanf(f3, "%d", &nll);
    fscanf(f3, "%s", vall);
    qi = atof(vall);
    *(F1 + nll - 1) = qi;
  }



  if ((f4 = fopen("ddl_i_mu.dat", "r")) == NULL)
  {
    perror("fopen 4");
    exit(4);
  }


  while (!feof(f4))
  {
    fscanf(f4, "%d", &nll);
    fscanf(f4, "%s", vall);
    qi = atof(vall);
    m = qi;
    *(meth_cfd.cfd.ddl_i + nll - 1) = m - 1;

  }



  if ((f5 = fopen("ddl_n_mu.dat", "r")) == NULL)
  {
    perror("fopen 5");
    exit(5);
  }


  while (!feof(f5))
  {
    fscanf(f5, "%d", &nll);
    fscanf(f5, "%s", vall);
    qi = atof(vall);
    m = qi;
    *(meth_cfd.cfd.ddl_n + nll - 1) = m - 1;
  }



  if ((f6 = fopen("ddl_t_mu.dat", "r")) == NULL)
  {
    perror("fopen 6");
    exit(6);
  }


  while (!feof(f6))
  {
    fscanf(f6, "%d", &nll);
    fscanf(f6, "%s", vall);
    qi = atof(vall);
    m = qi;
    *(meth_cfd.cfd.ddl_tt + nll - 1) = m - 1;

  }



  if ((f7 = fopen("ddl_c_mu.dat", "r")) == NULL)
  {
    perror("fopen 7");
    exit(7);
  }


  while (!feof(f7))
  {
    fscanf(f7, "%d", &nll);
    fscanf(f7, "%s", &vall);
    qi = atof(vall);
    m = qi;
    *(meth_cfd.cfd.ddl_c + nll - 1) = m - 1;
  }




  U2 = malloc(n * sizeof(double));
  F2 = malloc(n * sizeof(double));


  printf("\n\n  we go in the function name %s\n\n", nom);



  info = solve_cfd(K1, F1, &dimM, &meth_cfd, U2, F2);


  printf("\n\n we go out the function and info is %d\n", info);






  for (i = 0; i < 20; i++)
    printf("z %g w %g \n", U2[i], F2[i]);



  fclose(f1);
  fclose(f2);
  fclose(f3);
  fclose(f4);
  fclose(f5);
  fclose(f6);
  fclose(f7);


  free(vec);
  free(M);
  free(K1);
  free(F1);


  free(U2);
  free(F2);
  free(meth_cfd.cfd.ddl_i);
  free(meth_cfd.cfd.ddl_tt);
  free(meth_cfd.cfd.ddl_n);
  free(meth_cfd.cfd.ddl_c);
  free(meth_cfd.cfd.J1);


}








