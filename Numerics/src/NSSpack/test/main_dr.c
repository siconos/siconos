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
// This main file allows the dual resolution of relay problems:
//  try (z,w) such that:
//
//
//
//                    M z  +q = w                      (1)
//                    -z in diff(PSI) (w)              (2)
//                               [-b,a]
//
//
//  here M is an n by n  matrix, q an n-dimensional vector, z,w,b,a an n-dimensional vectors.
//
//  This system of equations and inequalities is solved thanks to rp_subroutine:
//        dr_gsnl  ( M, q, n, a, b, itermax, tol, chat, z, w, it_end, res, info)
//        dr_latin ( M, q, n, k_latin,  a, b, itermax, tol, chat, z, w, it_end, res, info)
//
//
//  where _ itermax      is the maximum iterations required, it's an integer
//        _ res          is the residue, it's a float (positive float)
//        _ it_end       is the number of iterations carried out, it's an integer
//        _ tol          is the tolerance value desired, it's a positive float (if you make it non positive then you force the convergence)
//        _ chat         is the output log identifiant
//        _ k_latin      is the parameter research of the latin, it's a float (strictly positive)
//        _ a            is the upper bound, it's a vector of floats
//        _ b            is the down bound, it's a vector of floats
//        _ z and w      are the solutions of the problem
//        _ info         shows the termination reason,0 is successful (otherwise see the termination reason of the solver used), it's an integer.
//
//
//    For more information about the methods see the chapter 4 of the Siconos manual theory.
//
//
//
//  The subroutine's call is due to the function solve_rd:
//
//  int solve_dr ( double *M, double *q, int n, methode *pt, double *z, double *w )
//
//  where M       is an n by n matrix,
//        q       an n-dimensional vector,
//        n       is the row dimension of M,
//        pt      a pointer other a union ( methode see "NSSpack.h").
//        z and w are n-dimensional  vectors solution.
//
//        methode is a variable with a union type; in this union you find the structure (method_dr) that gives to the function
//  solve_dr, the name and the parameters (itermax, tol, k_latin, a, b, chat) of the method we want to use.
//  This function return an interger:  0 successful .
//
//
///////////////////////////////////////////////////////////////////////////////


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "NSSpack.h"
#include "LA.h"

#define CHAT


int test_mmc(void)
{
  FILE      *f1, *f2, *f3, *f4;
  int       i, j, nl, nc, nll, n = 40;
  int       info1 = -1, info2 = -1, incx = 1, incy = 1;
  int       nonsymmetric, dim;

  double    *q, *z1, *w1, *z2, *w2, *vec, *a, *b;
  double    *c, *d, *qt, *wt1, *wt2, *abso1, *abso2;
  double    mini, max1, max2, maxi_1, maxi_2, max11, max22;
  double    qi, Mij, alpha, beta, num, den, diff1, diff2;
  double    comp11, comp22, comp111, comp222, comp1111, comp2222;


  char      val[50], vall[50], NT = 'N';
  method    meth_dr1, meth_dr2;


  printf("\n* *** ******************** *** * \n");
  printf("* ***        TEST MMC      *** * \n");
  printf("* *** ******************** *** * \n");



  /*           Allocations                     */

  q     = (double*)malloc(n * sizeof(double));
  z1    = (double*)malloc(n * sizeof(double));
  w1    = (double*)malloc(n * sizeof(double));
  z2    = (double*)malloc(n * sizeof(double));
  w2    = (double*)malloc(n * sizeof(double));
  a     = (double*)malloc(n * sizeof(double));
  b     = (double*)malloc(n * sizeof(double));
  vec   = (double*)malloc(n * n * sizeof(double));
  c     = (double*)malloc(n * sizeof(double));
  d     = (double*)malloc(n * sizeof(double));
  qt    = (double*)malloc(n * sizeof(double));
  wt1   = (double*)malloc(n * sizeof(double));
  wt2   = (double*)malloc(n * sizeof(double));
  abso1 = (double*)malloc(n * sizeof(double));
  abso2 = (double*)malloc(n * sizeof(double));


  /*     Data loading of M , q , a and b           */

  if ((f1 = fopen("DATA/M_relay_rd2.dat", "r")) == NULL)
  {
    perror("fopen 1");
    exit(1);
  }


  if ((f2 = fopen("DATA/q_relay_rd2.dat", "r")) == NULL)
  {
    perror("fopen 2");
    exit(2);
  }


  if ((f3 = fopen("DATA/a_relay_rd2.dat", "r")) == NULL)
  {
    perror("fopen 5");
    exit(5);
  }


  if ((f4 = fopen("DATA/b_relay_rd2.dat", "r")) == NULL)
  {
    perror("fopen 6");
    exit(6);
  }



  while (!feof(f1))
  {
    fscanf(f1, "%d", &nl);
    fscanf(f1, "%d", &nc);
    fscanf(f1, "%s", val);
    Mij = atof(val);


    vec[(nc - 1)*n + nl - 1 ] = Mij;

  }




  while (!feof(f2))
  {
    fscanf(f2, "%d", &nll);
    fscanf(f2, "%s", vall);
    qi = atof(vall);
    q[nll - 1] = -qi;

  }


  while (!feof(f3))
  {
    fscanf(f3, "%d", &nll);
    fscanf(f3, "%s", vall);
    qi = atof(vall);
    a[nll - 1] = qi;
  }


  while (!feof(f4))
  {
    fscanf(f4, "%d", &nll);
    fscanf(f4, "%s", vall);
    qi = atof(vall);
    b[nll - 1] = qi;
  }


  fclose(f2);
  fclose(f4);
  fclose(f1);
  fclose(f3);



  /*           End od the data loading        */


  nonsymmetric = 0;


  /*           Is M symmetric ?                */

  for (i = 0 ; i < n ; ++i)
  {
    for (j = 0 ; j < i ; ++j)
    {
      if (abs(vec[i * n + j] - vec[j * n + i]) > 1e-10)
      {
        nonsymmetric = 1;
        break;
      }
    }
  }



#ifdef CHAT
  if (nonsymmetric) printf("\n !! WARNING !!\n M is a non symmetric matrix \n \n");
  else printf(" \n M is a symmetric matrix \n \n");
#endif





  /*           Initialization                 */



  strcpy(meth_dr1.dr.name, "NLGS");
  meth_dr1.dr.itermax  =  1000;
  meth_dr1.dr.tol      =  0.00001;
  meth_dr1.dr.chat     =  1;


  strcpy(meth_dr2.dr.name, "Latin");
  meth_dr2.dr.itermax  =  5000;
  meth_dr2.dr.tol      =  0.0005;
  meth_dr2.dr.k_latin  =  0.09;
  meth_dr2.dr.chat     =  1;




  meth_dr1.dr.a = (double*)malloc(n * sizeof(double));
  meth_dr1.dr.b = (double*)malloc(n * sizeof(double));
  meth_dr2.dr.a = (double*)malloc(n * sizeof(double));
  meth_dr2.dr.b = (double*)malloc(n * sizeof(double));

  for (i = 0; i <= n - 1; i++)
  {

    meth_dr1.dr.a[i] =  a[i];
    meth_dr1.dr.b[i] = -b[i];
    meth_dr2.dr.a[i] =  a[i];
    meth_dr2.dr.b[i] = -b[i];


  }




#ifdef CHAT
  printf("**** NLGS TEST ****\n \n");
#endif

  info1 = dr_solver(vec, q, &n, &meth_dr1, z1, w1);

#ifdef CHAT
  printf("\n**** LATIN TEST ***\n \n");
#endif

  info2 = dr_solver(vec, q, &n, &meth_dr2, z2, w2);


  if (info1 >= info2)
  {
    return info1;
  }
  else return info2;






#ifdef CHAT
  printf(" *** ************************************** ***\n");
  for (i = 0 ; i < n ; ++i)
    printf(" NLGS RESULT : %14.7e %14.7e   LATIN RESULT : % 14.7e  %14.7e \n" , z1[i] , w1 [i], z2[i], w2[i]);

  printf("\n -- SOLVEUR --- ITER --- ERR-----INT-------NUL-z-----POS-z-----NPOS-z----\n");





  for (i = 0; i < n ; i ++)
  {

    c[i] = (meth_dr1.dr.a[i] - meth_dr1.dr.b[i]) / 2;
    d[i] = (meth_dr1.dr.a[i] + meth_dr1.dr.b[i]) / 2;

  }


  DCOPY(n, w1, incx, wt1, incy);
  DCOPY(n, w2, incx, wt2, incy);
  DCOPY(n,  q, incx, qt,  incy);

  alpha = -1.;
  DAXPY(n, alpha, c, incx, wt1, incy);

  alpha = -1.;
  DAXPY(n , alpha , c , incx , wt2 , incy);

  alpha = -1.;
  DAXPY(n , alpha , c , incx , qt , incy);


  /*  TEST of behavior laws */


  /*                 Relay                */


  /*      wt in interval [-d,d]       */

  abs_part(wt1, abso1, &n);

  abs_part(d, abso2, &n);

  alpha = -1.0;
  DAXPY(n , alpha , abso1 , incx , abso2 , incy);

  min_part(abso2, &mini, &n);

  mini = - mini;

  dim = 1;
  pos_part(&mini, &mini, &dim) ;


  abs_part(wt1, abso2, &n);

  max_part(abso2 , &max11 , &n);


  if (max11 > 1e-10)
  {

    max11 = mini / max11;

  }
  else
  {

    abs_part(wt1, abso2, &n);
    max_part(abso2, &maxi_1, &n);

    abs_part(qt, abso1, &n);
    max_part(abso1, &maxi_2, &n);

    if (maxi_1 > maxi_2)
      max11 = maxi_1;
    else
      max11 = maxi_2;

    max11 = mini / max11;

  }



  /*       Test of |z| = 0        */



  abs_part(wt1, abso1, &n);

  abs_part(d, abso2, &n);

  alpha = -1.0;
  DAXPY(n , alpha , abso1 , incx , abso2 , incy);

  abs_part(abso2, abso1, &n);

  abs_part(z1, abso2, &n);

  comp11 = DDOT(n , abso1 , incx , abso2 , incy);


  abs_part(z1, abso2, &n);
  max_part(abso2 , &max2 , &n);

  abs_part(wt1, abso2, &n);
  max_part(abso2 , &max1 , &n);


  if (max1 > 1e-10)
  {

    comp11 = comp11 / (n * max1 * max2);

  }
  else
  {

    abs_part(wt1, abso2, &n);
    max_part(abso2, &maxi_1, &n);

    abs_part(qt, abso1, &n);
    max_part(abso1, &maxi_2, &n);

    if (maxi_1 > maxi_2)
      max1 = maxi_1;
    else
      max1 = maxi_2;

    comp11 = comp11 / (n * max2 * max1);

  }



  /*       Test of z > 0       */


  DCOPY(n, wt1, incx, abso1, incy);

  alpha = 1.0;
  DAXPY(n , alpha , d , incx , abso1 , incy);

  abs_part(abso1, abso2, &n);

  pos_part(z1, abso1, &n) ;

  comp111 = DDOT(n , abso1 , incx , abso2 , incy);

  abs_part(z1, abso2, &n);

  max_part(abso2 , &max2 , &n);

  abs_part(wt1, abso2, &n);

  max_part(abso2 , &max1 , &n);



  if (max1 > 1e-10)
  {

    comp111 = comp111 / (n * max1 * max2);

  }
  else
  {

    abs_part(wt1, abso2, &n);
    max_part(abso2, &maxi_1, &n);

    abs_part(qt, abso1, &n);
    max_part(abso1, &maxi_2, &n);

    if (maxi_1 > maxi_2)
      max1 = maxi_1;
    else
      max1 = maxi_2;

    comp111 = comp111 / (n * max2 * max1);

  }






  /*        Test of z < 0         */


  DCOPY(n, wt1, incx, abso1, incy);

  alpha = -1.0;
  DAXPY(n , alpha , d , incx , abso1 , incy);

  abs_part(abso1, abso2, &n);

  DCOPY(n, z1, incx, abso1, incy);

  alpha = -1.;
  DSCAL(n , alpha , abso1 , incx);

  pos_part(abso1, abso1, &n) ;

  comp1111 = DDOT(n , abso1 , incx , abso2 , incy);

  abs_part(z1, abso2, &n);

  max_part(abso2 , &max2 , &n);

  abs_part(wt1, abso2, &n);

  max_part(abso2 , &max1 , &n);



  if (max1 > 1e-10)
  {

    comp1111 = comp1111 / (n * max1 * max2);

  }
  else
  {

    abs_part(wt1, abso2, &n);
    max_part(abso2, &maxi_1, &n);

    abs_part(qt, abso1, &n);
    max_part(abso1, &maxi_2, &n);

    if (maxi_1 > maxi_2)
      max1 = maxi_1;
    else
      max1 = maxi_2;

    comp1111 = comp1111 / (n * max2 * max1);

  }



  /*        Equilibrium       */


  alpha = -1;
  DAXPY(n , alpha , q , incx , w1 , incy);

  beta  = 1;
  DGEMV(NT , n , n , beta , vec , n , z1 , incx , alpha , w1 , incy);

  num = DNRM2(n , w1 , incx);
  den = DNRM2(n , q ,  incx);

  diff1 = num / den ;


  /*       TEST of behavior laws     */


  /*                 Relay                */


  /*      wt in interval [-d,d]       */

  abs_part(wt2, abso1, &n);

  abs_part(d, abso2, &n);

  alpha = -1.0;
  DAXPY(n , alpha , abso1 , incx , abso2 , incy);

  min_part(abso2, mini, n);

  mini = - mini;

  dim = 1;
  pos_part(&mini, &mini, &dim) ;


  abs_part(wt2, abso2, &n);

  max_part(abso2 , &max11 , &n);


  if (max1 > 1e-10)
  {

    max22 = mini / max11;

  }
  else
  {

    abs_part(wt2, abso2, &n);
    max_part(abso2, &maxi_1, &n);

    abs_part(qt, abso1, &n);
    max_part(abso1, &maxi_2, &n);

    if (maxi_1 > maxi_2)
      max22 = maxi_1;
    else
      max22 = maxi_2;

    max22 = mini / max22;

  }



  /*       Test of |z| = 0        */



  abs_part(wt1, abso1, &n);

  abs_part(d, abso2, &n);

  alpha = -1.0;
  DAXPY(n , alpha , abso1 , incx , abso2 , incy);

  abs_part(abso2, abso1, &n);

  abs_part(z2, abso2, &n);

  comp22 = DDOT(n , abso1 , incx , abso2 , incy);


  abs_part(z2, abso2, &n);

  max_part(abso2 , &max2 , &n);

  abs_part(wt2, abso2, &n);

  max_part(abso2 , &max1 , &n);


  if (max1 > 1e-10)
  {

    comp22 = comp22 / (n * max1 * max2);

  }
  else
  {

    abs_part(wt2, abso2, &n);
    max_part(abso2, &maxi_1, &n);

    abs_part(qt, abso1, &n);
    max_part(abso1, &maxi_2, &n);

    if (maxi_1 > maxi_2)
      max1 = maxi_1;
    else
      max1 = maxi_2;

    comp22 = comp22 / (n * max2 * max1);

  }



  /*       Test of z > 0       */


  DCOPY(n, wt2, incx, abso1, incy);

  alpha = 1.0;
  DAXPY(n , alpha , d , incx , abso1 , incy);

  abs_part(abso1, abso2, &n);

  pos_part(z2, abso1, &n) ;

  comp222 = DDOT(n , abso1 , incx , abso2 , incy);

  abs_part(z2, abso2, &n);

  max_part(abso2 , &max2 , &n);

  abs_part(wt2, abso2, &n);

  max_part(abso2 , &max1 , &n);



  if (max1 > 1e-10)
  {

    comp222 = comp222 / (n * max1 * max2);

  }
  else
  {

    abs_part(wt2, abso2, &n);
    max_part(abso2, &maxi_1, &n);

    abs_part(qt, abso1, &n);
    max_part(abso1, &maxi_2, &n);

    if (maxi_1 > maxi_2)
      max1 = maxi_1;
    else
      max1 = maxi_2;

    comp222 = comp222 / (n * max2 * max1);

  }






  /*        Test of z < 0         */


  DCOPY(n, wt2, incx, abso1, incy);

  alpha = -1.0;
  DAXPY(n , alpha , d , incx , abso1 , incy);

  abs_part(abso1, abso2, &n);

  DCOPY(n, z2, incx, abso1, incy);

  alpha = -1.;
  DSCAL(n , alpha , abso1 , incx);

  pos_part(abso1, abso1, &n) ;

  comp2222 = DDOT(n , abso1 , incx , abso2 , incy);

  abs_part(z2, abso2, &n);

  max_part(abso2 , &max2 , &n);

  abs_part(wt2, abso2, &n);

  max_part(abso2 , &max1 , &n);



  if (max1 > 1e-10)
  {

    comp2222 = comp2222 / (n * max1 * max2);

  }
  else
  {

    abs_part(wt2, abso2, &n);
    max_part(abso2, &maxi_1, &n);

    abs_part(qt, abso1, &n);
    max_part(abso1, &maxi_2, &n);

    if (maxi_1 > maxi_2)
      max1 = maxi_1;
    else
      max1 = maxi_2;

    comp2222 = comp2222 / (n * max2 * max1);

  }




  /*        Equilibrium       */

  alpha = -1;
  DAXPY(n , alpha , q , incx , w2 , incy);

  beta  = 1;
  DGEMV(LA_NOTRANS , n , n , beta , vec , n , z2 , incx , alpha , w2 , incy);

  num = DNRM2(n , w2 , incx);
  den = DNRM2(n , q , incx);

  diff2 = num / den ;


  /*       PRINT RESULT       */

  printf("\n  NLGS RESULT:%5d|%10.4e|%10.4e|%10.4e|%10.4e|%10.4e|\n" , meth_dr1.dr.iter, diff1, max11, comp11, comp111, comp1111);
  printf(" LATIN RESULT:%5d|%10.4e|%10.4e|%10.4e|%10.4e|%10.4e|\n" , meth_dr2.dr.iter , diff2, max22, comp22, comp222, comp2222);





  printf(" *** ************************************** ***\n");
#endif


  free(vec);
  free(q);
  free(z1);
  free(w1);
  free(z2);
  free(w2);
  free(a);
  free(b);
  free(c);
  free(d);
  free(qt);
  free(wt1);
  free(wt2);
  free(abso1);
  free(abso2);
  free(meth_dr1.dr.a);
  free(meth_dr1.dr.b);
  free(meth_dr2.dr.a);
  free(meth_dr2.dr.b);
}


int main(void)
{
  return test_mmc();

}


