///////////////////////////////////////////////////////////////////////////
//  This main file allows the resolution of LCP (Linear Complementary Problem):
//  try (z,w) such that:
//
//                        M z - w = q        (1)
//                        0<= z , 0<= w      (2)
//                        z.w= O             (3)
//
//  here M is an n by n  matrix, q an n-dimensional vector, z an n-dimensional
//  vector and w an n-dimensional vector.
//
//  This system of equations and inequalities is solved thanks to subroutine_lcp:
//        lemke_lcp (M,q,n,itermax,z,w,it_end,res,info)
//        gsnl_lcp(M,q,n,itermax,tol,z,w,it_end,res,info)
//        gcp_lcp(M,q,n,itermax,tol,z,w,it_end,res,info)
//        latin_lcp (M,q,n,k_latin,itermax,tol,z,w,it_end,res,info)
//        qp_lcp(M,q,n,res,z,w,fail) (not available for the moment)
//
//  where _ itermax is the maximum iterations required, it's an integer
//        _ res is the residue, it's a float
//        _ it_end is the number of iterations carried out, it's an integer
//        _ tol is the tolerance value desired, it's a float
//        _ k_latin is the parameter research of the latin, it's a float
//        _ fail shows the termination reason, it's an integer
//        _ info shows the trmination reason (0 successful otherwise 1), it's an integer.
//
//  The subroutine's call is due to the function solve_lcp:
//
//  int solve_lcp (float (*M)[maxcols],float * q, int n, methode *pt, float *z, float * w)
//
//  where M is an n by n matrix, q an n-dimensional vector, n is the row dimension of M, and pt a pointer other a structure (methode), z and w are n-dimensional vector, the solutions of the lcp.
//  methode is a variable with a structure type; this structure gives to the function solve_lcp, the name and the parameters (itermax, tol, k_latin) of the method we want to use.
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
  int i, j, nl, nc, nll, it, info, n = 31, dimM = n;
  double *q, *z, *w, *vec;
  double(*M)[n];
  double qi, Mij;
  char val[14], vall[14];

  static methode_lcp meth_lcp  = {"Gcp", 101, 0.0001, 0.6};

  if ((f1 = fopen("MM_mmc.dat", "r")) == NULL)
  {
    perror("fopen 1");
    exit(1);
  }

  if ((f2 = fopen("qq_mmc.dat", "r")) == NULL)
  {
    perror("fopen 2");
    exit(2);
  }

  //  f3=fopen("aff_M.dat","w+");
  M = malloc(dimM * dimM * sizeof(double));
  vec = (double*)malloc(dimM * dimM * sizeof(double));

  for (i = 0; i < dimM; i++)
    for (j = 0; j < dimM; j++)
      M[i][j] = 0.;

  while (!feof(f1))
  {
    fscanf(f1, "%d", &nl);
    fscanf(f1, "%d", &nc);
    fscanf(f1, "%s", val);
    Mij = atof(val);

    /////////////       on met la transpos       ////////////////
    /*fprintf(f3,"%d %d %.14e\n",nc,nl,Mij);
      fscanf(f3,"%d %d %.14e\n", &nc, &nl, &Mij);*/
    *(*(M + nc - 1) + nl - 1) = Mij;
    //////////////         fin transpos         ////////////////////
  }


  //// valeurs du tableau dans vec (compatibilite allocation memoire f90)///
  for (i = 0; i < dimM; i++)
    for (j = 0; j < dimM; j++)
      vec[j * dimM + i] = M[i][j];
  //       printf("vec(%d) = %.14e \n",i*dimM+j,M[i][j]);}
  ////////////////////////////////////////////////////////////////////////



  if ((f2 = fopen("qq_mmc.dat", "r")) == NULL)
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
    *(q + nll - 1) = -qi;
  }


  printf("\n we go in the function\n");

  info = solve_lcp(vec, q, &n, &meth_lcp, z, w);

  printf("\n we go out the function and info is %d\n", info);

  fclose(f2);
  fclose(f1);

  free(M);
  free(vec);
  free(q);
  free(z);
  free(w);

  /////////////////////////////////
  // second test //////////////////
  /////////////////////////////////
  static methode_lcp meth_lcp2  = {"Gsnl", 101, 0.0001, 0.6};

  if ((f1 = fopen("MM_mmc.dat", "r")) == NULL)
  {
    perror("fopen 1");
    exit(1);
  }

  if ((f2 = fopen("qq_mmc.dat", "r")) == NULL)
  {
    perror("fopen 2");
    exit(2);
  }

  //  f3=fopen("aff_M.dat","w+");
  M = malloc(dimM * dimM * sizeof(double));
  vec = (double*)malloc(dimM * dimM * sizeof(double));

  for (i = 0; i < dimM; i++)
  {
    for (j = 0; j < dimM; j++)
    {
      M[i][j] = 0.;
    }
  }

  while (!feof(f1))
  {
    fscanf(f1, "%d", &nl);
    fscanf(f1, "%d", &nc);
    fscanf(f1, "%s", val);
    Mij = atof(val);

    /////////////       on met la transpos       ////////////////
    /*fprintf(f3,"%d %d %.14e\n",nc,nl,Mij);
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
      //       printf("vec(%d) = %.14e \n",i*dimM+j,M[i][j]);}
    }
  }
  ////////////////////////////////////////////////////////////////////////



  if ((f2 = fopen("qq_mmc.dat", "r")) == NULL)
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
    *(q + nll - 1) = -qi;
  }


  printf("\n we go in the function\n");

  info = solve_lcp(vec, q, &n, &meth_lcp2, z, w);

  printf("\n we go out the function and info is %d\n", info);

  fclose(f2);
  fclose(f1);

  free(M);
  free(vec);
  free(q);
  free(z);
  free(w);



  /////////////////////////////////
  // third test //////////////////
  /////////////////////////////////
  static methode_lcp meth_lcp3  = {"Latin", 101, 0.0001, 0.6};

  if ((f1 = fopen("MM_mmc.dat", "r")) == NULL)
  {
    perror("fopen 1");
    exit(1);
  }

  if ((f2 = fopen("qq_mmc.dat", "r")) == NULL)
  {
    perror("fopen 2");
    exit(2);
  }

  //  f3=fopen("aff_M.dat","w+");
  M = malloc(dimM * dimM * sizeof(double));
  vec = (double*)malloc(dimM * dimM * sizeof(double));

  for (i = 0; i < dimM; i++)
  {
    for (j = 0; j < dimM; j++)
    {
      M[i][j] = 0.;
    }
  }

  while (!feof(f1))
  {
    fscanf(f1, "%d", &nl);
    fscanf(f1, "%d", &nc);
    fscanf(f1, "%s", val);
    Mij = atof(val);

    /////////////       on met la transpos       ////////////////
    /*fprintf(f3,"%d %d %.14e\n",nc,nl,Mij);
      fscanf(f3,"%d %d %.14e\n", &nc, &nl, &Mij);*/
    //*(*(M+nc-1)+nl-1)=Mij;

    // pas de transposée maintenant
    *(*(M + nl - 1) + nc - 1) = Mij;
    //////////////         fin transpos         ////////////////////

  }


  //// valeurs du tableau dans vec (compatibilite allocation memoire f90)///
  //    for (i=0;i<dimM;i++)
  //      {for (j=0;j<dimM;j++){
  //  vec[j*dimM+i]= M[i][j];
  //  //       printf("vec(%d) = %.14e \n",i*dimM+j,M[i][j]);}
  //      }}
  ////////////////////////////////////////////////////////////////////////



  if ((f2 = fopen("qq_mmc.dat", "r")) == NULL)
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
    *(q + nll - 1) = -qi;
  }


  printf("\n we go in the function\n");

  info = solve_lcp(/*vec*/M, q, &n, &meth_lcp3, z, w);

  printf("\n we go out the function and info is %d\n", info);

  fclose(f2);
  fclose(f1);

  free(M);
  free(vec);
  free(q);
  free(z);
  free(w);
}



