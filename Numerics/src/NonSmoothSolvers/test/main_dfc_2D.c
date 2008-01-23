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
//  This system of equations and inequalities is solved thanks to dfc_2D_subroutine:
//        dfc_2D_latin (M,q,n,k_latin,mu,itermax,tol,z,w,it_end,res,info)
//        ....
//        ....
//  or thanks to LCP (Linear Complementary Problem) routines after a new
//  formulation of this problem in the shape of LCP due to the dfc_2D2lcp and lcp2dfc_2D routines:
//
//       dfc_2D2lcp (nc,mu,M,q,Mtel,qtel)
//       lcp2dfc_2D (nc,ztel,wtel,z,w)
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
//  The subroutine's call is due to the function dfc_2D_driver:
//
//   int dfc_2D_driver (double (*M)[maxcols],double *q,int n,method *pt,double z[],double w[])
//
//  where M is an n by n matrix, q an n-dimensional vector, n is the row
//  dimension of M, and pt a pointer other a structure ( method).
//  method is a variable with a structure type; this structure gives to the
//  function solve_lcp, the name and the parameters (itermax, tol, k_latin,..)
//  of the method we want to use.
//  This function return an interger:  0 successful return otherwise 1.
//
//
///////////////////////////////////////////////////////////////////////////////


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "NonSmoothDrivers.h"

int main(void)
{

  FILE *f1, *f2, *f3, *f4, *f5, *f6;
  int i, nl, nc, nll, n = 596, dimM = n, dim_d = 17, dim_n = 25, dim_tt = 25, info = -7;
  int  m;
  double *vec, *K1, *F1, *U2, *F2;
  double *M;
  double qi, Mij;
  char val[50], vall[50];




  method meth_dfc_2D;

  // strcpy( meth_dfc_2D.dfc_2D.name, "NLGS");
  strcpy(meth_dfc_2D.dfc_2D.name, "Cfd_latin");
  meth_dfc_2D.dfc_2D.itermax = 1501;
  meth_dfc_2D.dfc_2D.tol = 0.000001;
  meth_dfc_2D.dfc_2D.mu = 0.5;
  meth_dfc_2D.dfc_2D.chat = 1;

  meth_dfc_2D.dfc_2D.k_latin = 0.6;
  meth_dfc_2D.dfc_2D.dim_d = dim_d;
  meth_dfc_2D.dfc_2D.dim_tt = dim_tt;




  if ((f1 = fopen("DATA/K1_mu.dat", "r")) == NULL)
  {
    perror("fopen 1");
    exit(1);
  }


  M = (double *)malloc(dimM * dimM * sizeof(double));
  vec = (double*)malloc(dimM * dimM * sizeof(double));
  K1 = (double*)malloc(dimM * dimM * sizeof(double));

  while (!feof(f1))
  {
    fscanf(f1, "%d", &nl);
    fscanf(f1, "%d", &nc);
    fscanf(f1, "%s", val);
    Mij = atof(val);



    K1 [(nc - 1)*dimM + nl - 1] = Mij;

  }



  if ((f2 = fopen("DATA/J1_mu.dat", "r")) == NULL)
  {
    perror("fopen 2");
    exit(2);
  }




  meth_dfc_2D.dfc_2D.J1 = (double *) malloc(dimM * sizeof(double));

  F1 = (double *) malloc(dimM * sizeof(double));

  meth_dfc_2D.dfc_2D.ddl_d = (int*)malloc(dim_d * sizeof(int));
  meth_dfc_2D.dfc_2D.ddl_n = (int*)malloc(dim_n * sizeof(int));
  meth_dfc_2D.dfc_2D.ddl_tt = (int*)malloc(dim_tt * sizeof(int));




  while (!feof(f2))
  {
    fscanf(f2, "%d", &nll);
    fscanf(f2, "%s", vall);
    qi = atof(vall);
    *(meth_dfc_2D.dfc_2D.J1 + nll - 1) = qi;
  }



  if ((f3 = fopen("DATA/F1_mu.dat", "r")) == NULL)
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



  if ((f4 = fopen("DATA/ddl_d_mu.dat", "r")) == NULL)
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
    *(meth_dfc_2D.dfc_2D.ddl_d + nll - 1) = m - 1;

  }



  if ((f5 = fopen("DATA/ddl_n_mu.dat", "r")) == NULL)
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
    *(meth_dfc_2D.dfc_2D.ddl_n + nll - 1) = m - 1;
  }



  if ((f6 = fopen("DATA/ddl_t_mu.dat", "r")) == NULL)
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
    *(meth_dfc_2D.dfc_2D.ddl_tt + nll - 1) = m - 1;

  }






  U2 = (double *)malloc(n * sizeof(double));
  F2 = (double*)malloc(n * sizeof(double));

  for (i = 0; i < n; i++)
  {
    U2[i] = 0.0;
    F2[i] = 0.0;

  }

  printf("\n\n  we go in the function  %s\n\n", meth_dfc_2D.dfc_2D.name);



  info = dfc_2D_driver(K1, F1, &dimM, &meth_dfc_2D, U2, F2);


  printf("\n\n we go out the function and info is %d\n", info);






  for (i = 0; i < 25; i++)
    printf("z %g w %g \n", U2[i], F2[i]);



  fclose(f1);
  fclose(f2);
  fclose(f3);
  fclose(f4);
  fclose(f5);
  fclose(f6);


  free(vec);
  free(M);
  free(K1);
  free(F1);


  free(U2);
  free(F2);
  free(meth_dfc_2D.dfc_2D.ddl_d);
  free(meth_dfc_2D.dfc_2D.ddl_tt);
  free(meth_dfc_2D.dfc_2D.ddl_n);

  free(meth_dfc_2D.dfc_2D.J1);

  return info;

}








