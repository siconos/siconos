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
//  This system of equations and inequalities is solved thanks to d_subroutine_cf:
//        cfd_latin (M,q,n,k_latin,mu,itermax,tol,z,w,it_end,res,info)
//        ....
//        ....
//  or thanks to LCP (Linear Complementary Problem) routines after a new
//  formulation of this problem in the shape of LCP due to the d_lcp_cf and d_cf_lcp routines:
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
#include "SiconosNumerics.h"


main()

{
  FILE *f1, *f2, *f3, *f4, *f5, *f6;
  int i, j, nl, nc, nll, it, info, n = 50, dimM = n;
  double *q, *z, *w, *vec;
  double(*M)[n];
  double qi, Mij;
  char val[14], vall[14];


  //////////////////////////////////////////
  // first test ///////////////////////////
  //////////////////////////////////////////

  static methode_cfd meth_cfd  = {"Latin", 101, 0.0001, 0.5, 7};
  //  static methode_cfd meth_cfd  = {"gsnl_lcp",501, 0.0001,0.5,7};
  //  static methode_cfd meth_cfd  = {"gcp_lcp",501, 0.0001,0.5,7};

  if ((f1 = fopen("MM_mmc_mu.dat", "r")) == NULL)
  {
    perror("fopen 1");
    exit(1);
  }


  //  f3=fopen("aff_M.dat","w+");
  M = malloc(dimM * dimM * sizeof(double));
  vec = (double*)malloc(dimM * dimM * sizeof(double));

  while (!feof(f1))
  {
    fscanf(f1, "%d", &nl);
    fscanf(f1, "%d", &nc);
    fscanf(f1, "%s", val);
    Mij = atof(val);

    /////////////       on met la transpos       ////////////////
    /*   fprintf(f3,"%d %d %.14e\n",nc,nl,Mij);
    fscanf(f3,"%d %d %.14e\n", &nc, &nl, &Mij);*/
    //*(*(M+nc-1)+nl-1)=Mij;

    *(*(M + nl - 1) + nc - 1) = Mij;
    //////////////         fin transpos         ////////////////////

  }


  //// valeurs du tableau dans vec (compatibilite allocation memoire f90)///
  /* for (i=0;i<dimM;i++)
     {for (j=0;j<dimM;j++){
  vec[j*dimM+i]= M[i][j];
     }}
     */
  ////////////////////////////////////////////////////////////////////////


  if ((f2 = fopen("qq_mmc_mu.dat", "r")) == NULL)
  {
    perror("fopen 2");
    exit(2);
  }


  //  f4=fopen("aff_q.dat","w+");
  q = malloc(dimM * sizeof(double));
  z = malloc(dimM * sizeof(double));
  w = malloc(dimM * sizeof(double));

  while (!feof(f2))
  {
    fscanf(f2, "%d", &nll);
    fscanf(f2, "%s", vall);
    qi = atof(vall);
    //fprintf(f4,"%d %.14e\n",nll,qi);
    *(q + nll - 1) = qi;
  }


  printf("\n\n  we go in the function\n\n");

  info = solve_cfd(/*vec*/M, q, &n, &meth_cfd, z, w);

  printf("\n\n we go out the function and info is %d\n", info);

  fclose(f2);
  fclose(f1);

  free(w);
  free(q);
  free(z);
  free(vec);
  free(M);


  //////////////////////////////////////////
  // second test ///////////////////////////
  //////////////////////////////////////////
  static methode_cfd meth_cfd1  = {"Gsnl", 101, 0.0001, 0.5, 7};
  if ((f1 = fopen("MM_mmc_mu.dat", "r")) == NULL)
  {
    perror("fopen 1");
    exit(1);
  }


  //  f3=fopen("aff_M.dat","w+");
  M = malloc(dimM * dimM * sizeof(double));
  vec = (double*)malloc(dimM * dimM * sizeof(double));

  while (!feof(f1))
  {
    fscanf(f1, "%d", &nl);
    fscanf(f1, "%d", &nc);
    fscanf(f1, "%s", val);
    Mij = atof(val);

    /////////////       on met la transpos       ////////////////
    /*   fprintf(f3,"%d %d %.14e\n",nc,nl,Mij);
    fscanf(f3,"%d %d %.14e\n", &nc, &nl, &Mij);*/
    *(*(M + nc - 1) + nl - 1) = Mij;
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


  if ((f2 = fopen("qq_mmc_mu.dat", "r")) == NULL)
  {
    perror("fopen 2");
    exit(2);
  }


  //  f4=fopen("aff_q.dat","w+");
  q = malloc(dimM * sizeof(double));
  z = malloc(dimM * sizeof(double));
  w = malloc(dimM * sizeof(double));

  while (!feof(f2))
  {
    fscanf(f2, "%d", &nll);
    fscanf(f2, "%s", vall);
    qi = atof(vall);
    //fprintf(f4,"%d %.14e\n",nll,qi);
    *(q + nll - 1) = qi;
  }


  printf("\n\n  we go in the function\n\n");

  info = solve_cfd(vec, q, &n, &meth_cfd1, z, w);

  printf("\n\n we go out the function and info is %d\n", info);

  fclose(f2);
  fclose(f1);

  free(w);
  free(q);
  free(z);
  free(vec);
  free(M);






  //////////////////////////////////////////
  // third test ///////////////////////////
  //////////////////////////////////////////
  static methode_cfd meth_cfd2  = {"Gcp", 101, 0.0001, 0.5, 7};
  if ((f1 = fopen("MM_mmc_mu.dat", "r")) == NULL)
  {
    perror("fopen 1");
    exit(1);
  }


  //  f3=fopen("aff_M.dat","w+");
  M = malloc(dimM * dimM * sizeof(double));
  vec = (double*)malloc(dimM * dimM * sizeof(double));

  while (!feof(f1))
  {
    fscanf(f1, "%d", &nl);
    fscanf(f1, "%d", &nc);
    fscanf(f1, "%s", val);
    Mij = atof(val);

    /////////////       on met la transpos       ////////////////
    /*   fprintf(f3,"%d %d %.14e\n",nc,nl,Mij);
    fscanf(f3,"%d %d %.14e\n", &nc, &nl, &Mij);*/
    *(*(M + nc - 1) + nl - 1) = Mij;
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


  if ((f2 = fopen("qq_mmc_mu.dat", "r")) == NULL)
  {
    perror("fopen 2");
    exit(2);
  }


  //  f4=fopen("aff_q.dat","w+");
  q = malloc(dimM * sizeof(double));
  z = malloc(dimM * sizeof(double));
  w = malloc(dimM * sizeof(double));

  while (!feof(f2))
  {
    fscanf(f2, "%d", &nll);
    fscanf(f2, "%s", vall);
    qi = atof(vall);
    //fprintf(f4,"%d %.14e\n",nll,qi);
    *(q + nll - 1) = qi;
  }


  printf("\n\n  we go in the function\n\n");

  info = solve_cfd(vec, q, &n, &meth_cfd2, z, w);

  printf("\n\n we go out the function and info is %d\n", info);

  fclose(f2);
  fclose(f1);

  free(w);
  free(q);
  free(z);
  free(vec);
  free(M);
}













