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
// This main file allows the primal resolution of relay problems:
//  try (z,w) such that:
//
//
//
//                    M z  + q = w                     (1)
//                    -w in diff(PSI) (z)              (2)
//                               [-b,a]
//
//
//  here M is an n by n  matrix, q an n-dimensional vector, z,w,b,a an n-dimensional vectors.
//
//  This system of equations and inequalities is solved thanks to rp_subroutine:
//        pr_gsnl  ( M, q, n, a, b, itermax, tol, chat, z, w, it_end, res, info)
//        pr_latin ( M, q, n, k_latin, a, b, itermax, tol, chat, z, w, it_end, res, info)
//
//
//  where _ itermax    is the maximum iterations required, it's an integer
//        _ res        is the residue, it's a float (positive float)
//        _ it_end     is the number of iterations carried out, it's an integer
//        _ tol        is the tolerance value desired, it's a positive float (if you make it non positive then you force the convergence)
//        _ k_latin    is the parameter research of the latin, it's a float (strictly positive)
//        _ a          is the upper bound, it's a vector of floats
//        _ b          is the down bound, it's a vector of floats
//        _ chat       is the output log identifiant
//        _ z and w    are the solutions of the problem
//        _ info       shows the termination reason,0 is successful (otherwise see the termination reason of the solver used), it's an integer.
//
//
//    For more information about the methods see the chapter 4 of the Siconos manual theory.
//
//
//
//  The subroutine's call is due to the function solve_rp:
//
//  int pr_solver ( double *M, double *q, int n, methode *pt, double *z, double *w )
//
//  where M        is an n by n matrix,
//        q        an n-dimensional vector,
//        n        is the row dimension of M,
//        pt       a pointer other a union ( methode see "NSSpack.h").
//        z and w are n-dimensional vectors solution.
//
//  methode is a variable with a union type; in this union you find the structure (method_rp) that gives to the function
// rp_solver, the name and the parameters ( itermax, tol, k_latin, a, b, chat) of the method we want to use.
//  This function return an interger:  0 for successful
//
//
///////////////////////////////////////////////////////////////////////////////


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "NSSpack.h"
#include "blaslapack.h"

#define CHAT


int test_mmc(void)
{
  FILE    *f1, *f2, *f3, *f4;
  int     i, j, nl, nc, nll, n = 15;
  int     info1 = -1, info2 = -1, incx = 1, incy = 1, dim;
  int     nonsymmetric;

  double  *q, *z1, *w1, *z2, *w2, *vec, *a, *b, *qt;
  double  *zt1, *zt2, *abso1, *abso2, *d, *c, *abso3, mini, max1, max2, maxi_1, maxi_2;
  double  qi, Mij, diff1, diff2, num, den, alpha, beta, max11, max22;
  double  comp11, comp22, comp111, comp222, comp1111, comp2222;

  char    val[50], vall[50], NT = 'T';
  method  meth_pr1, meth_pr2;


  printf("\n* *** ******************** *** * \n");
  printf("* ***        TEST MMC      *** * \n");
  printf("* *** ******************** *** * \n");



  /*           Allocations                     */

  q       = (double*) malloc(n * sizeof(double));
  qt      = (double*) malloc(n * sizeof(double));
  z1      = (double*) malloc(n * sizeof(double));
  w1      = (double*) malloc(n * sizeof(double));
  z2      = (double*) malloc(n * sizeof(double));
  w2      = (double*) malloc(n * sizeof(double));
  a       = (double*) malloc(n * sizeof(double));
  b       = (double*) malloc(n * sizeof(double));
  vec     = (double*) malloc(n * n * sizeof(double));
  zt1     = (double*) malloc(n * sizeof(double));
  zt2     = (double*) malloc(n * sizeof(double));
  abso1   = (double*) malloc(n * sizeof(double));
  abso2   = (double*) malloc(n * sizeof(double));
  d       = (double*) malloc(n * sizeof(double));
  abso3   = (double*) malloc(n * sizeof(double));
  c       = (double*) malloc(n * sizeof(double));


  /*     Data loading of M , q , a and b           */

  if ((f1 = fopen("DATA/M_relay1.dat", "r")) == NULL)
  {
    perror("fopen 1");
    exit(1);
  }


  if ((f2 = fopen("DATA/q_relay1.dat", "r")) == NULL)
  {
    perror("fopen 2");
    exit(2);
  }


  if ((f3 = fopen("DATA/a_relay1.dat", "r")) == NULL)
  {
    perror("fopen 5");
    exit(5);
  }


  if ((f4 = fopen("DATA/b_relay1.dat", "r")) == NULL)
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
      if (abs(vec[i * n + j] - vec[j * n + i]) > 1e-16)
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

  strcpy(meth_pr1.pr.name, "NLGS");
  meth_pr1.pr.itermax  =  50000;
  meth_pr1.pr.tol      =  0.000001;
  meth_pr1.pr.chat     =  1;


  strcpy(meth_pr2.pr.name, "Latin");
  meth_pr2.pr.itermax  =  50000;
  meth_pr2.pr.tol      =  0.0000001;
  meth_pr2.pr.k_latin  =  0.003;//0.003;
  meth_pr2.pr.chat     =  1;




  meth_pr1.pr.a = (double*)malloc(n * sizeof(double));
  meth_pr1.pr.b = (double*)malloc(n * sizeof(double));
  meth_pr2.pr.a = (double*)malloc(n * sizeof(double));
  meth_pr2.pr.b = (double*)malloc(n * sizeof(double));

  for (i = 0; i <= n - 1; i++)
  {

    meth_pr1.pr.a[i] =  a[i];
    meth_pr1.pr.b[i] = -b[i];
    meth_pr2.pr.a[i] =  a[i];
    meth_pr2.pr.b[i] = -b[i];


  }




#ifdef CHAT
  printf("**** NLGS TEST ****\n \n");
#endif

  info1 = pr_solver(vec, q, &n, &meth_pr1, z1, w1);

#ifdef CHAT
  printf("\n**** LATIN TEST ***\n \n");
#endif

  info2 = pr_solver(vec, q, &n, &meth_pr2, z2, w2);


  if (info1 >= info2)
  {
    return info1;
  }
  else return info2;



#ifdef CHAT
  printf(" *** ************************************** ***\n");
  for (i = 0 ; i < n ; ++i)
    printf(" NLGS RESULT : %10.4e  %10.4e| LATIN RESULT : % 10.4e  %10.4e \n" , z1[i], w1[i], z2[i], w2[i]);

  printf(" -- SOLVEUR --- ITER --- ERR-------INT-------NUL-w-----POS-w-----NPOS-w-----\n");




  /*  TEST of behavior laws */



  for (i = 0; i < n ; i ++)
  {

    c[i] = (meth_pr1.pr.a[i] - meth_pr1.pr.b[i]) / 2;
    d[i] = (meth_pr1.pr.a[i] + meth_pr1.pr.b[i]) / 2;

  }


  dcopy_(&n, z1, &incx, zt1, &incy);
  dcopy_(&n, z2, &incx, zt2, &incy);

  alpha = -1.;
  daxpy_(&n , &alpha , c , &incx , zt1 , &incy);

  alpha = -1.;
  daxpy_(&n , &alpha , c , &incx , zt2 , &incy);

  alpha = 1.;
  beta  = 0.;
  dgemv_(&NT , &n , &n , &beta , vec , &n , q , &incx , &alpha , qt , &incy);



  /*                 Relay                */


  /*      zt in interval [-d,d]       */

  abs_part(zt1, abso1, &n);

  abs_part(d, abso2, &n);

  alpha = -1.0;
  daxpy_(&n , &alpha , abso1 , &incx , abso2 , &incy);

  min_part(abso2, &mini, &n);

  mini = - mini;

  dim = 1;
  pos_part(&mini, &mini, &dim) ;


  abs_part(zt1, abso2, &n);

  max_part(abso2 , &max11 , &n);

  max11 = mini / max11;


  /*       Test of |w| = 0        */



  abs_part(zt1, abso1, &n);

  abs_part(d, abso2, &n);

  alpha = -1.0;
  daxpy_(&n , &alpha , abso1 , &incx , abso2 , &incy);

  abs_part(abso2, abso1, &n);

  abs_part(w1, abso2, &n);

  comp11 = ddot_(&n , abso1 , &incx , abso2 , &incy);


  abs_part(zt1, abso2, &n);

  max_part(abso2 , &max2 , &n);

  abs_part(w1, abso2, &n);

  max_part(abso2 , &max1 , &n);


  if (max1 > 1e-10)
  {

    comp11 = comp11 / (n * max1 * max2);

  }
  else
  {

    abs_part(w1, abso2, &n);
    max_part(abso2, &maxi_1, &n);

    abs_part(q, abso1, &n);
    max_part(abso1, &maxi_2, &n);

    if (maxi_1 > maxi_2)
      max1 = maxi_1;
    else
      max1 = maxi_2;

    comp11 = comp11 / (n * max2 * max1);

  }



  /*       Test of wt > 0       */


  dcopy_(&n, zt1, &incx, abso1, &incy);

  alpha = 1.0;
  daxpy_(&n , &alpha , d , &incx , abso1 , &incy);

  abs_part(abso1, abso2, &n);

  pos_part(w1, abso1, &n) ;

  comp111 = ddot_(&n , abso1 , &incx , abso2 , &incy);

  abs_part(zt1, abso2, &n);

  max_part(abso2 , &max2 , &n);

  abs_part(w1, abso2, &n);

  max_part(abso2 , &max1 , &n);



  if (max1 > 1e-10)
  {

    comp111 = comp111 / (n * max1 * max2);

  }
  else
  {

    abs_part(w1, abso2, &n);
    max_part(abso2, &maxi_1, &n);

    abs_part(q, abso1, &n);
    max_part(abso1, &maxi_2, &n);

    if (maxi_1 > maxi_2)
      max1 = maxi_1;
    else
      max1 = maxi_2;

    comp111 = comp111 / (n * max2 * max1);

  }






  /*        Test of wt < 0         */


  dcopy_(&n, zt1, &incx, abso1, &incy);

  alpha = -1.0;
  daxpy_(&n , &alpha , d , &incx , abso1 , &incy);

  abs_part(abso1, abso2, &n);

  dcopy_(&n, w1, &incx, abso1, &incy);

  alpha = -1.;
  dscal_(&n , &alpha , abso1 , &incx);

  pos_part(abso1, abso1, &n) ;

  comp1111 = ddot_(&n , abso1 , &incx , abso2 , &incy);

  abs_part(zt1, abso2, &n);

  max_part(abso2 , &max2 , &n);

  abs_part(w1, abso2, &n);

  max_part(abso2 , &max1 , &n);



  if (max1 > 1e-10)
  {

    comp1111 = comp1111 / (n * max1 * max2);

  }
  else
  {

    abs_part(w1, abso2, &n);
    max_part(abso2, &maxi_1, &n);

    abs_part(q, abso1, &n);
    max_part(abso1, &maxi_2, &n);

    if (maxi_1 > maxi_2)
      max1 = maxi_1;
    else
      max1 = maxi_2;

    comp1111 = comp1111 / (n * max2 * max1);

  }



  /*                    Equilibrium                  */



  alpha = -1;
  daxpy_(&n , &alpha , q , &incx , w1 , &incy);

  beta  = 1;
  dgemv_(&NT , &n , &n , &beta , vec , &n , z1 , &incx , &alpha , w1 , &incy);

  num = dnrm2_(&n , w1 , &incx);
  den = dnrm2_(&n , q , &incx);

  diff1 = num / den ;



  /*                 Relay        */


  /*      zt in interval [-d,d]       */


  abs_part(zt2, abso1, &n);

  abs_part(d, abso2, &n);

  alpha = -1.0;
  daxpy_(&n , &alpha , abso1 , &incx , abso2 , &incy);

  min_part(abso2, &mini, &n);

  mini = - mini;

  dim = 1;
  pos_part(&mini, &mini, &dim) ;


  abs_part(zt2, abso2, &n);

  max_part(abso2 , &max22 , &n);

  max22 = mini / max22;


  /*        Test of |w| = 0        */



  abs_part(zt2, abso1, &n);

  abs_part(d, abso2, &n);

  alpha = -1.0;
  daxpy_(&n , &alpha , abso1 , &incx , abso2 , &incy);

  abs_part(abso2, abso1, &n);

  abs_part(w2, abso2, &n);

  comp22 = ddot_(&n , abso1 , &incx , abso2 , &incy);


  abs_part(zt2, abso2, &n);

  max_part(abso2 , &max2 , &n);

  abs_part(w2, abso2, &n);

  max_part(abso2 , &max1 , &n);



  if (max1 > 1e-10)
  {

    comp22 = comp22 / (n * max1 * max2);

  }
  else
  {

    abs_part(w2, abso2, &n);
    max_part(abso2, &maxi_1, &n);

    abs_part(q, abso1, &n);
    max_part(abso1, &maxi_2, &n);

    if (maxi_1 > maxi_2)
      max1 = maxi_1;
    else
      max1 = maxi_2;

    comp22 = comp22 / (n * max2 * max1);

  }



  /*       Test of wt > 0       */



  dcopy_(&n, zt2, &incx, abso1, &incy);

  alpha = 1.0;
  daxpy_(&n , &alpha , d , &incx , abso1 , &incy);

  abs_part(abso1, abso2, &n);

  pos_part(w2, abso1, &n) ;

  comp222 = ddot_(&n , abso1 , &incx , abso2 , &incy);

  abs_part(zt2, abso2, &n);

  max_part(abso2 , &max2 , &n);

  abs_part(w2, abso2, &n);

  max_part(abso2 , &max1 , &n);



  if (max1 > 1e-10)
  {

    comp222 = comp222 / (n * max1 * max2);

  }
  else
  {

    abs_part(w2, abso2, &n);
    max_part(abso2, &maxi_1, &n);

    abs_part(q, abso1, &n);
    max_part(abso1, &maxi_2, &n);

    if (maxi_1 > maxi_2)
      max1 = maxi_1;
    else
      max1 = maxi_2;

    comp222 = comp222 / (n * max2 * max1);

  }

  /*        Test of wt < 0       */


  dcopy_(&n, zt2, &incx, abso1, &incy);

  alpha = -1.0;
  daxpy_(&n , &alpha , d , &incx , abso1 , &incy);

  abs_part(abso1, abso2, &n);

  dcopy_(&n, w2, &incx, abso1, &incy);

  alpha = -1.;
  dscal_(&n , &alpha , abso1 , &incx);

  pos_part(abso1, abso1, &n) ;

  comp2222 = ddot_(&n , abso1 , &incx , abso2 , &incy);

  abs_part(zt2, abso2, &n);

  max_part(abso2 , &max2 , &n);

  abs_part(w2, abso2, &n);

  max_part(abso2 , &max1 , &n);



  if (max1 > 1e-10)
  {

    comp2222 = comp2222 / (n * max1 * max2);

  }
  else
  {

    abs_part(w2, abso2, &n);
    max_part(abso2, &maxi_1, &n);

    abs_part(q, abso1, &n);
    max_part(abso1, &maxi_2, &n);

    if (maxi_1 > maxi_2)
      max1 = maxi_1;
    else
      max1 = maxi_2;

    comp2222 = comp2222 / (n * max2 * max1);

  }


  /*                    Equilibrium                  */

  alpha = -1;
  daxpy_(&n , &alpha , q , &incx , w2 , &incy);

  beta  = 1;
  dgemv_(&NT , &n , &n , &beta , vec , &n , z2 , &incx , &alpha , w2 , &incy);

  num = dnrm2_(&n , w2 , &incx);
  den = dnrm2_(&n , q , &incx);

  diff2 = num / den ;



  printf("\n  NLGS RESULT: %5d|%10.4e|%10.4e|%10.4e|%10.4e|%10.4e|\n" , meth_pr1.pr.iter ,  diff1, max11, comp11, comp111, comp1111);
  printf(" LATIN RESULT: %5d|%10.4e|%10.4e|%10.4e|%10.4e|%10.4e|\n" , meth_pr2.pr.iter ,  diff2, max22, comp22, comp222, comp2222);





  printf(" *** ************************************** ***\n");
#endif


  free(vec);
  free(q);
  free(qt);
  free(z1);
  free(w1);
  free(z2);
  free(w2);
  free(a);
  free(b);
  free(zt1);
  free(zt2);
  free(abso1);
  free(abso2);
  free(d);
  free(abso3);
  free(c);
  free(meth_pr1.pr.a);
  free(meth_pr1.pr.b);
  free(meth_pr2.pr.a);
  free(meth_pr2.pr.b);
}


int main(void)
{
  return test_mmc();
}


