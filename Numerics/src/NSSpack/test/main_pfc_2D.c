/* Siconos-Numerics version 2.1.1, Copyright INRIA 2005-2006.
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
/*!
 *
 *  This main file allows the primal resolution of contact problems with friction:
 *  try (z,w) such that:
 *
 *
 *                    M z  + q = w                 (1)
 *                    0<= zn , 0<= wn , zn.wn= O   (2) (Signorini contact law)
 *                    -wt in diff(PSI)   (zt)      (3) (Coulomb friction law)
 *                                 [-mu*zn, mu*zn]
 *
 *
 *  here M is an n by n  matrix, q an n-dimensional vector, z an n-dimensional
 *  vector and w an n-dimensional vector.
 *
 *  This system of equations and inequalities is solved thanks to pfc_2D_subroutine:
 *        pfc_2D_nlgs (M,q,n,mu,itermax,tol,z,w,it_end,res,info)
 *        pfc_2D_cpg (M,q,n,mu,itermax,tol,z,w,it_end,res,info)
 *        pfc_2D_latin (M,q,n,k_latin,mu,itermax,tol,z,w,it_end,res,info)
 *
 *  where _ itermax is the maximum iterations required, it's an integer
 *        _ mu is the friction coefficient, it's a float
 *        _ res is the residue, it's a float
 *        _ it_end is the number of iterations carried out, it's an integer
 *        _ tol is the tolerance value desired, it's a float
 *        _ k_latin is the parameter research of the latin, it's a float
 *        _ z and w are the solutions of the problem
 *        _ info shows the termination reason,0 is successful otherwise 1, it's an integer.
 *
 *
 * \author Shéhérazade Nineb.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "NSSpack.h"
#include "blaslapack.h"

#define BAVARD

void pfc_2D_series(int n , double *vec , double *q)
{

  int         i, j;
  int         nonsymmetric;
  int         info[3];
  int         nc;
  int         incx = 1, incy = 1, dim;

  double      den, num, alpha, beta , diff;
  double      *z1, *w1, *z2, *w2, *z3, *w3;
  double      *wn1, *zn1, *wt1, *zt1, *wn2;
  double      *zn2, *wt2, *zt2, *wn3, *zn3, *wt3, *zt3;
  double      *muzn1, *muzn2, *muzn3;
  double      mini, *abso, maxi, comp1, comp2, comp3, *absow, *absoz;
  double      max11, max22, max33, maxi_1, maxi_2;
  double      max1, max2, *abso1, *abso2;
  double      r1, r2, comp11, comp22, comp33, comp111, comp222, comp333;
  double      comp1111, comp2222, comp3333;

  char        NT = 'N';

  method      meth_pfc_2D1, meth_pfc_2D2, meth_pfc_2D3;



  printf(" START SERIES \n");


  for (i = 0 ; i < 3 ; ++i) info[i] = -1;



  strcpy(meth_pfc_2D1.pfc_2D.name, "NLGS");
  meth_pfc_2D1.pfc_2D.itermax  =  10000;
  meth_pfc_2D1.pfc_2D.tol      =  0.0000001;
  meth_pfc_2D1.pfc_2D.chat     =  1;
  meth_pfc_2D1.pfc_2D.mu       =  0.3;


  strcpy(meth_pfc_2D2.pfc_2D.name, "CPG");
  meth_pfc_2D2.pfc_2D.itermax  =  7000;
  meth_pfc_2D2.pfc_2D.tol      =  0.0000001;
  meth_pfc_2D2.pfc_2D.chat     =  1;
  meth_pfc_2D2.pfc_2D.mu       =  0.3;


  strcpy(meth_pfc_2D3.pfc_2D.name, "Latin");
  meth_pfc_2D3.pfc_2D.itermax  =  15000;
  meth_pfc_2D3.pfc_2D.tol      =  0.0000001;
  meth_pfc_2D3.pfc_2D.chat     =  1;
  meth_pfc_2D3.pfc_2D.mu       =  0.3;
  meth_pfc_2D3.pfc_2D.k_latin  =  5.5;




  z1     = malloc(n * sizeof(double));
  w1     = malloc(n * sizeof(double));
  z2     = malloc(n * sizeof(double));
  w2     = malloc(n * sizeof(double));
  z3     = malloc(n * sizeof(double));
  w3     = malloc(n * sizeof(double));
  absow  = malloc(n * sizeof(double));
  absoz  = malloc(n * sizeof(double));



  wn1    = malloc(n / 2 * sizeof(double));
  zn1    = malloc(n / 2 * sizeof(double));
  wt1    = malloc(n / 2 * sizeof(double));
  zt1    = malloc(n / 2 * sizeof(double));
  wn2    = malloc(n / 2 * sizeof(double));
  zn2    = malloc(n / 2 * sizeof(double));
  wt2    = malloc(n / 2 * sizeof(double));
  zt2    = malloc(n / 2 * sizeof(double));
  wn3    = malloc(n / 2 * sizeof(double));
  zn3    = malloc(n / 2 * sizeof(double));
  wt3    = malloc(n / 2 * sizeof(double));
  zt3    = malloc(n / 2 * sizeof(double));


  muzn1  = malloc(n / 2 * sizeof(double));
  muzn2  = malloc(n / 2 * sizeof(double));
  muzn3  = malloc(n / 2 * sizeof(double));

  abso   = malloc(n / 2 * sizeof(double));
  abso1  = malloc(n / 2 * sizeof(double));
  abso2  = malloc(n / 2 * sizeof(double));


  nonsymmetric = 0;


  /*         Is M symmetric ?         */

  for (i = 0 ; i < n ; ++i)
  {
    for (j = 0 ; j < i ; ++j)
    {
      if (abs(vec[i * n + j] - vec[j * n + i]) > 1e-12)
      {
        nonsymmetric = 1;
        break;
      }
    }
  }



#ifdef BAVARD
  if (nonsymmetric) printf("\n !! WARNING !!\n M is a non symmetric matrix \n");
  else printf(" M is a symmetric matrix \n");
#endif


#ifdef BAVARD
  printf("**** NLGS TEST ****\n");
#endif
  for (i = 0 ; i < n ; ++i)
  {
    z1[i] = 0.0;
    w1[i] = 0.0;
  }

  info[0] = pfc_2D_solver(vec , q , &n , &meth_pfc_2D1 , z1 , w1);


#ifdef BAVARD
  printf("**** CPG TEST *****\n");
#endif
  for (i = 0 ; i < n ; ++i)
  {
    z2[i] = 0.0;
    w2[i] = 0.;
  }

  info[1] = pfc_2D_solver(vec , q , &n , &meth_pfc_2D2 , z2 , w2);


#ifdef BAVARD
  printf("**** Latin TEST ***\n");
#endif
  for (i = 0 ; i < n ; ++i)
  {
    z3[i] = 0.0;
    w3[i] = 0.0;
  }

  info[2] = pfc_2D_solver(vec , q , &n , &meth_pfc_2D3 , z3 , w3);


#ifdef BAVARD
  printf(" *** ************************************** ***\n");

  for (i = 0 ; i < n ; i++)
    printf("\n   NLGS RESULT : %10.4g  %10.4g |CPG : %10.4g  %10.4g |LATIN : %10.4g  %10.4g" , z1[i], w1[i], z2[i], w2[i], z3[i], w3[i]);



  printf("\n\n");


  /*  TEST of behavior laws  and  equilibrium */


  nc = n / 2;

  for (i = 0; i < nc ; i ++)
  {

    wn1[i] = w1[2 * i];
    wt1[i] = w1[2 * i + 1];

    zn1[i] = z1[2 * i];
    zt1[i] = z1[2 * i + 1];


    wn2[i] = w2[2 * i];
    wt2[i] = w2[2 * i + 1];

    zn2[i] = z2[2 * i];
    zt2[i] = z2[2 * i + 1];

    wn3[i] = w3[2 * i];
    wt3[i] = w3[2 * i + 1];

    zn3[i] = z3[2 * i];
    zt3[i] = z3[2 * i + 1];

  }


  printf(" -- SOLVEUR -- ITER----- ERR ----POS--zn--------wn-------COMP------INT------NUL-wt----POS-wt----NPOS-wt--\n");



  /*        Complementary of normal part     */


  /*         zn            */

  min_part(zn1, &mini, &nc);

  mini = - mini;

  dim = 1;
  pos_part(&mini, &mini, &dim) ;

  abs_part(zn1, abso, &nc);

  max_part(abso , &maxi , &nc);

  r1 = mini / maxi;

  /*       wn                */

  min_part(wn1, &mini, &nc);

  mini = - mini;

  pos_part(&mini, &mini, &dim) ;

  abs_part(wn1, abso, &nc);

  max_part(abso , &maxi , &nc);


  if (maxi > 1e-10)
  {

    r2 = mini / maxi;

  }
  else
  {

    abs_part(w1, absow, &n);
    max_part(absow, &maxi_1, &n);

    abs_part(q, absow, &n);
    max_part(absow, &maxi_2, &n);

    if (maxi_1 > maxi_2)
      maxi = maxi_1;
    else
      maxi = maxi_2;

    r2 = mini / maxi;
  }



  /*        zn^t wn                  */

  abs_part(wn1, abso1, &nc);

  abs_part(zn1, abso2, &nc);


  comp1 = ddot_(&nc , abso1 , &incx , abso2 , &incy);

  max_part(abso1 , &max1 , &nc);

  max_part(abso2 , &max2 , &nc);



  if (max1 > 1e-10)
  {

    comp1 = comp1 / (nc * max2 * max1);


  }
  else
  {

    abs_part(w1, absow, &n);
    max_part(absow, &maxi_1, &n);

    abs_part(q, absow, &n);
    max_part(absow, &maxi_2, &n);

    if (maxi_1 > maxi_2)
      max1 = maxi_1;
    else
      max1 = maxi_2;

    comp1 = comp1 / (nc * max2 * max1);

  }



  /*                  Friction                */



  /*           Test in interval             */


  abs_part(zt1, abso1, &nc);

  dcopy_(&nc, zn1, &incx, muzn1, &incy);

  alpha = 0.3;
  dscal_(&nc , &alpha , muzn1 , &incx);

  abs_part(muzn1, abso2, &nc);


  alpha = -1.0;
  daxpy_(&nc , &alpha , abso1 , &incx , abso2 , &incy);

  min_part(abso2, &mini, &nc);

  mini = - mini;

  pos_part(&mini, &mini, &dim) ;


  abs_part(z1, absoz, &n);

  max_part(absoz , &max11 , &n);

  max11 = mini / max11;


  /*            Test of |wt| = 0               */



  abs_part(zt1, abso1, &nc);

  dcopy_(&nc, zn1, &incx, muzn1, &incy);

  alpha = 0.3;
  dscal_(&nc , &alpha , muzn1 , &incx);

  abs_part(muzn1, abso2, &nc);


  alpha = -1.0;
  daxpy_(&nc , &alpha , abso1 , &incx , abso2 , &incy);

  abs_part(abso2, abso1, &nc);

  abs_part(wt1, abso2, &nc);



  comp11 = ddot_(&nc , abso1 , &incx , abso2 , &incy);


  abs_part(z1, absoz, &n);

  max_part(absoz , &max2 , &n);




  abs_part(w1, absow, &n);

  max_part(absow , &max1 , &n);



  if (max1 > 1e-10)
  {


    comp11 = comp11 / (nc * max1 * max2);

  }
  else
  {

    abs_part(w1, absow, &n);
    max_part(absow, &maxi_1, &n);

    abs_part(q, absow, &n);
    max_part(absow, &maxi_2, &n);

    if (maxi_1 > maxi_2)
      max1 = maxi_1;
    else
      max1 = maxi_2;

    comp11 = comp11 / (nc * max1 * max2);

  }




  /*             Test of wt > 0           */

  dcopy_(&nc, zn1, &incx, muzn1, &incy);

  alpha = 0.3;
  dscal_(&nc , &alpha , muzn1 , &incx);

  alpha = -1.0;
  daxpy_(&nc , &alpha , zt1 , &incx , muzn1 , &incy);

  abs_part(muzn1, abso1, &nc);

  pos_part(wt1, abso2, &nc) ;

  comp111 = ddot_(&nc , abso1 , &incx , abso2 , &incy);

  abs_part(z1, absoz, &n);

  max_part(absoz , &max2 , &n);

  abs_part(w1, absow, &n);

  max_part(absow , &max1 , &n);


  if (max1 > 1e-10)
  {


    comp111 = comp111 / (nc * max1 * max2);

  }
  else
  {

    abs_part(w1, absow, &n);
    max_part(absow, &maxi_1, &n);

    abs_part(q, absow, &n);
    max_part(absow, &maxi_2, &n);

    if (maxi_1 > maxi_2)
      max1 = maxi_1;
    else
      max1 = maxi_2;

    comp111 = comp111 / (nc * max1 * max2);

  }



  /*             Test of wt < 0               */



  dcopy_(&nc, zn1, &incx, muzn1, &incy);

  alpha = 0.3;
  dscal_(&nc , &alpha , muzn1 , &incx);

  alpha = 1.0;
  daxpy_(&nc , &alpha , zt1 , &incx , muzn1 , &incy);

  abs_part(muzn1, abso1, &nc);

  dcopy_(&nc, wt1, &incx, muzn1, &incy);

  alpha = -1.;
  dscal_(&nc , &alpha , muzn1 , &incx);

  pos_part(muzn1, abso2, &nc) ;

  comp1111 = ddot_(&nc , abso1 , &incx , abso2 , &incy);

  abs_part(z1, absoz, &n);

  max_part(absoz , &max2 , &n);

  abs_part(w1, absow, &n);

  max_part(absow , &max1 , &n);



  if (max1 > 1e-10)
  {


    comp1111 = comp1111 / (nc * max1 * max2);

  }
  else
  {

    abs_part(w1, absow, &n);
    max_part(absow, &maxi_1, &n);

    abs_part(q, absow, &n);
    max_part(absow, &maxi_2, &n);

    if (maxi_1 > maxi_2)
      max1 = maxi_1;
    else
      max1 = maxi_2;

    comp1111 = comp1111 / (nc * max1 * max2);

  }


  /*                 Equilibrium                   */

  alpha = -1;
  daxpy_(&n , &alpha , q , &incx , w1 , &incy);

  beta  = 1;
  dgemv_(&NT , &n , &n , &beta , vec , &n , z1 , &incx , &alpha , w1 , &incy);

  num = dnrm2_(&n , w1 , &incx);
  den = dnrm2_(&n , q , &incx);

  diff = num / den ;



  printf("\n  NLGS (LOG:%1d)|%5d|%7.4e|%7.4e|%7.4e|%7.4e|%7.4e|%7.4e|%7.4e|%7.4e|" , info[0] , meth_pfc_2D1.pfc_2D.iter , diff , r1, r2, comp1, max11, comp11, comp111, comp1111);


  /*                Complementary of normal part        */


  /*              zn             */

  min_part(zn2, &mini, &nc);

  mini = - mini;

  pos_part(&mini, &mini, &dim) ;

  abs_part(zn2, abso, &nc);

  max_part(abso , &maxi , &nc);

  r1 = mini / maxi;


  /*              wn                 */


  min_part(wn2, &mini, &nc);

  mini = - mini;

  pos_part(&mini, &mini, &dim) ;

  abs_part(wn2, abso, &nc);

  max_part(abso , &maxi , &nc);


  if (maxi > 1e-10)
  {

    r2 = mini / maxi;

  }
  else
  {

    abs_part(w2, absow, &n);
    max_part(absow, &maxi_1, &n);

    abs_part(q, absow, &n);
    max_part(absow, &maxi_2, &n);

    if (maxi_1 > maxi_2)
      maxi = maxi_1;
    else
      maxi = maxi_2;

    r2 = mini / maxi;
  }


  /*          zn^t wn              */


  abs_part(wn2, abso1, &nc);

  abs_part(zn2, abso2, &nc);

  comp2 = ddot_(&nc , abso1 , &incx , abso2 , &incy);

  max_part(abso1 , &max1 , &nc);

  max_part(abso2 , &max2 , &nc);


  if (max1 > 1e-10)
  {

    comp2 = comp2 / (nc * max2 * max1);


  }
  else
  {

    abs_part(w2, absow, &n);
    max_part(absow, &maxi_1, &n);

    abs_part(q, absow, &n);
    max_part(absow, &maxi_2, &n);

    if (maxi_1 > maxi_2)
      max1 = maxi_1;
    else
      max1 = maxi_2;

    comp2 = comp2 / (nc * max2 * max1);

  }


  /*                  Friction                    */



  /*           Test in interval             */

  abs_part(zt2, abso1, &nc);

  dcopy_(&nc, zn2, &incx, muzn2, &incy);

  alpha = 0.3;
  dscal_(&nc , &alpha , muzn2 , &incx);

  abs_part(muzn2, abso2, &nc);


  alpha = -1.0;
  daxpy_(&nc , &alpha , abso1 , &incx , abso2 , &incy);

  min_part(abso2, &mini, &nc);

  mini = - mini;

  pos_part(&mini, &mini, &dim) ;


  abs_part(z2, absow, &n);

  max_part(absow , &max22 , &n);

  max22 = mini / max22;

  /*                  Test of |wt| = 0                 */


  abs_part(zt2, abso1, &nc);

  dcopy_(&nc, zn2, &incx, muzn2, &incy);

  alpha = 0.3;
  dscal_(&nc , &alpha , muzn2 , &incx);

  abs_part(muzn2, abso2, &nc);


  alpha = -1.0;
  daxpy_(&nc , &alpha , abso1 , &incx , abso2 , &incy);

  abs_part(abso2, abso1, &nc);

  abs_part(wt2, abso2, &nc);



  comp22 = ddot_(&nc , abso1 , &incx , abso2 , &incy);


  abs_part(z2, absoz, &n);

  max_part(absoz , &max2 , &n);

  abs_part(w2, absow, &n);

  max_part(absow , &max1 , &n);


  if (max1 > 1e-10)
  {


    comp22 = comp22 / (nc * max1 * max2);

  }
  else
  {

    abs_part(w2, absow, &n);
    max_part(absow, &maxi_1, &n);

    abs_part(q, absow, &n);
    max_part(absow, &maxi_2, &n);

    if (maxi_1 > maxi_2)
      max1 = maxi_1;
    else
      max1 = maxi_2;

    comp22 = comp22 / (nc * max1 * max2);

  }



  /*             Test of wt > 0           */


  dcopy_(&nc, zn2, &incx, muzn2, &incy);

  alpha = 0.3;
  dscal_(&nc , &alpha , muzn2 , &incx);

  alpha = -1.0;
  daxpy_(&nc , &alpha , zt2 , &incx , muzn2 , &incy);

  abs_part(muzn2, abso1, &nc);

  pos_part(wt2, abso2, &nc) ;

  comp222 = ddot_(&nc , abso1 , &incx , abso2 , &incy);

  abs_part(z2, absoz, &n);

  max_part(absoz , &max2 , &n);

  abs_part(w2, absow, &n);

  max_part(absow , &max1 , &n);



  if (max1 > 1e-10)
  {


    comp222 = comp222 / (nc * max1 * max2);

  }
  else
  {

    abs_part(w2, absow, &n);
    max_part(absow, &maxi_1, &n);

    abs_part(q, absow, &n);
    max_part(absow, &maxi_2, &n);

    if (maxi_1 > maxi_2)
      max1 = maxi_1;
    else
      max1 = maxi_2;

    comp222 = comp222 / (nc * max1 * max2);

  }


  /*           Test of wt < 0             */


  dcopy_(&nc, zn2, &incx, muzn2, &incy);

  alpha = 0.3;
  dscal_(&nc , &alpha , muzn2 , &incx);

  alpha = 1.0;
  daxpy_(&nc , &alpha , zt2 , &incx , muzn2 , &incy);

  abs_part(muzn2, abso1, &nc);

  dcopy_(&nc, wt2, &incx, muzn2, &incy);

  alpha = -1.;
  dscal_(&nc , &alpha , muzn2 , &incx);

  pos_part(muzn2, abso2, &nc) ;

  comp2222 = ddot_(&nc , abso1 , &incx , abso2 , &incy);

  abs_part(z2, absoz, &n);

  max_part(absoz , &max2 , &n);

  abs_part(w2, absow, &n);

  max_part(absow , &max1 , &n);


  if (max1 > 1e-10)
  {


    comp2222 = comp2222 / (nc * max1 * max2);

  }
  else
  {

    abs_part(w2, absow, &n);
    max_part(absow, &maxi_1, &n);

    abs_part(q, absow, &n);
    max_part(absow, &maxi_2, &n);

    if (maxi_1 > maxi_2)
      max1 = maxi_1;
    else
      max1 = maxi_2;

    comp2222 = comp2222 / (nc * max1 * max2);

  }



  /*           Equilibrium                   */

  alpha = -1;
  daxpy_(&n , &alpha , q , &incx , w2 , &incy);

  beta  = 1;
  dgemv_(&NT , &n , &n , &beta , vec , &n , z2 , &incx , &alpha , w2 , &incy);

  num = dnrm2_(&n , w2 , &incx);
  den = dnrm2_(&n , q , &incx);

  diff = num / den ;



  printf("\n   CPG (LOG:%1d)|%5d|%7.4e|%7.4e|%7.4e|%7.4e|%7.4e|%7.4e|%7.4e|%7.4e|" , info[1] , meth_pfc_2D2.pfc_2D.iter , diff , r1, r2, comp2, max22, comp22, comp222, comp2222);


  /*        Complementary of normal part            */

  min_part(zn3, &mini, &nc);

  mini = - mini;

  pos_part(&mini, &mini, &dim) ;

  abs_part(zn3, abso, &nc);

  max_part(abso , &maxi , &nc);

  r1 = mini / maxi;


  min_part(wn3, &mini, &nc);

  mini = - mini;

  pos_part(&mini, &mini, &dim) ;

  abs_part(wn3, abso, &nc);

  max_part(abso , &maxi , &nc);


  if (maxi > 1e-10)
  {

    r2 = mini / maxi;

  }
  else
  {

    abs_part(w3, absow, &n);
    max_part(absow, &maxi_1, &n);

    abs_part(q, absow, &n);
    max_part(absow, &maxi_2, &n);

    if (maxi_1 > maxi_2)
      maxi = maxi_1;
    else
      maxi = maxi_2;

    r2 = mini / maxi;
  }



  abs_part(wn3, abso1, &nc);

  abs_part(zn3, abso2, &nc);

  comp3 = ddot_(&nc , abso1 , &incx , abso2 , &incy);

  max_part(abso1 , &max1 , &nc);

  max_part(abso2 , &max2 , &nc);



  if (max1 > 1e-10)
  {

    comp3 = comp3 / (nc * max2 * max1);


  }
  else
  {

    abs_part(w3, absow, &n);
    max_part(absow, &maxi_1, &n);

    abs_part(q, absow, &n);
    max_part(absow, &maxi_2, &n);

    if (maxi_1 > maxi_2)
      max1 = maxi_1;
    else
      max1 = maxi_2;

    comp3 = comp3 / (nc * max2 * max1);

  }




  /*                   Friction                 */


  /*         Test in interval           */

  abs_part(zt3, abso1, &nc);

  dcopy_(&nc, zn3, &incx, muzn3, &incy);

  alpha = 0.3;
  dscal_(&nc , &alpha , muzn3 , &incx);

  abs_part(muzn3, abso2, &nc);


  alpha = -1.0;
  daxpy_(&nc , &alpha , abso1 , &incx , abso2 , &incy);

  min_part(abso2, &mini, &nc);

  mini = - mini;

  pos_part(&mini, &mini, &dim) ;


  abs_part(z3, absoz, &n);

  max_part(absoz , &max33 , &n);

  max33 = mini / max33;



  /*            Test of |wt| = 0               */

  abs_part(zt3, abso1, &nc);

  dcopy_(&nc, zn3, &incx, muzn3, &incy);

  alpha = 0.3;
  dscal_(&nc , &alpha , muzn3 , &incx);

  abs_part(muzn3, abso2, &nc);


  alpha = -1.0;
  daxpy_(&nc , &alpha , abso1 , &incx , abso2 , &incy);

  abs_part(abso2, abso1, &nc);

  abs_part(wt3, abso2, &nc);

  comp33 = ddot_(&nc , abso1 , &incx , abso2 , &incy);


  abs_part(z3, absoz, &n);

  max_part(absoz , &max2 , &n);

  abs_part(w3, absow, &n);

  max_part(absow , &max1 , &n);



  if (max1 > 1e-10)
  {


    comp33 = comp33 / (nc * max1 * max2);

  }
  else
  {

    abs_part(w3, absow, &n);
    max_part(absow, &maxi_1, &n);

    abs_part(q, absow, &n);
    max_part(absow, &maxi_2, &n);

    if (maxi_1 > maxi_2)
      max1 = maxi_1;
    else
      max1 = maxi_2;

    comp33 = comp33 / (nc * max1 * max2);

  }



  /*          Test of wt > 0            */


  dcopy_(&nc, zn3, &incx, muzn3, &incy);

  alpha = 0.3;
  dscal_(&nc , &alpha , muzn3 , &incx);

  alpha = -1.0;
  daxpy_(&nc , &alpha , zt3 , &incx , muzn3 , &incy);

  abs_part(muzn3, abso1, &nc);

  pos_part(wt3, abso2, &nc) ;

  comp333 = ddot_(&nc , abso1 , &incx , abso2 , &incy);

  abs_part(z3, absoz, &n);

  max_part(absoz , &max2 , &n);

  abs_part(w3, absow, &n);

  max_part(absow , &max1 , &n);


  if (max1 > 1e-10)
  {


    comp333 = comp333 / (nc * max1 * max2);

  }
  else
  {

    abs_part(w3, absow, &n);
    max_part(absow, &maxi_1, &n);

    abs_part(q, absow, &n);
    max_part(absow, &maxi_2, &n);

    if (maxi_1 > maxi_2)
      max1 = maxi_1;
    else
      max1 = maxi_2;

    comp333 = comp333 / (nc * max1 * max2);

  }





  /*            Test of wt < 0           */



  dcopy_(&nc, zn3, &incx, muzn3, &incy);

  alpha = 0.3;
  dscal_(&nc , &alpha , muzn3 , &incx);

  alpha = 1.0;
  daxpy_(&nc , &alpha , zt3 , &incx , muzn3 , &incy);

  abs_part(muzn3, abso1, &nc);

  dcopy_(&nc, wt3, &incx, muzn3, &incy);

  alpha = -1.;
  dscal_(&nc , &alpha , muzn3 , &incx);

  pos_part(muzn3, abso2, &nc) ;

  comp3333 = ddot_(&nc , abso1 , &incx , abso2 , &incy);

  abs_part(z3, absoz, &n);

  max_part(absoz , &max2 , &n);

  abs_part(w3, absow, &n);

  max_part(absow , &max1 , &n);


  if (max1 > 1e-10)
  {


    comp3333 = comp3333 / (nc * max1 * max2);

  }
  else
  {

    abs_part(w3, absow, &n);
    max_part(absow, &maxi_1, &n);

    abs_part(q, absow, &n);
    max_part(absow, &maxi_2, &n);

    if (maxi_1 > maxi_2)
      max1 = maxi_1;
    else
      max1 = maxi_2;

    comp3333 = comp3333 / (nc * max1 * max2);

  }


  /*              Equilibrium                   */

  alpha = -1;
  daxpy_(&n , &alpha , q , &incx , w3 , &incy);

  beta  = 1;
  dgemv_(&NT , &n , &n , &beta , vec , &n , z3 , &incx , &alpha , w3 , &incy);


  num = dnrm2_(&n , w3 , &incx);
  den = dnrm2_(&n , q , &incx);


  diff = num / den ;


  printf("\n LATIN (LOG:%1d)|%5d|%7.4e|%7.4e|%7.4e|%7.4e|%7.4e|%7.4e|%7.4e|%7.4e| \n \n" , info[2] , meth_pfc_2D3.pfc_2D.iter , diff , r1, r2, comp3, max33, comp33, comp333, comp3333);



#endif

  free(z1);
  free(w1);
  free(z2);
  free(w2);
  free(z3);
  free(w3);
  free(zn1);
  free(wn1);
  free(zt1);
  free(wt1);
  free(zn2);
  free(wn2);
  free(zt2);
  free(wt2);
  free(zn3);
  free(wn3);
  free(zt3);
  free(wt3);
  free(abso);
  free(abso1);
  free(abso2);

  free(absow);
  free(absoz);

  free(muzn1);
  free(muzn2);
  free(muzn3);

}

/* ****************************************************** */
/*          READ MATRICES OF MATRIX DIRECTORY             */
/* ****************************************************** */

void test_matrix(void)
{

  FILE *LCPfile;

  int i, j, itest, NBTEST;
  int isol;
  int dim , dim2;

  double *q, *sol;
  double *vecM;

  char val[20];

  NBTEST = 2;

  for (itest = 0 ; itest < NBTEST ; ++itest)
  {

    switch (itest)
    {
    case 0:
      printf("\n\n 2x2 LCP \n");
      if ((LCPfile = fopen("MATRIX/deudeu.dat", "r")) == NULL)
      {
        perror("fopen LCPfile: deudeu.dat");
        exit(1);
      }
      break;
    case 1:
      printf("\n\n DIAGONAL LCP \n");
      if ((LCPfile = fopen("MATRIX/trivial2.dat", "r")) == NULL)
      {
        perror("fopen LCPfile: trivial2.dat");
        exit(1);
      }
      break;
    }

    fscanf(LCPfile , "%d" , &dim);

    dim2 = dim * dim;
    isol = 1;

    vecM = (double*)malloc(dim2 * sizeof(double));
    q    = (double*)malloc(dim * sizeof(double));
    sol  = (double*)malloc(dim * sizeof(double));

    for (i = 0 ; i < dim ; ++i)  q[i]    = 0.0;
    for (i = 0 ; i < dim ; ++i)  sol[i]  = 0.0;
    for (i = 0 ; i < dim2 ; ++i) vecM[i] = 0.0;

    for (i = 0 ; i < dim ; ++i)
    {
      for (j = 0 ; j < dim ; ++j)
      {
        fscanf(LCPfile, "%s", val);
        vecM[ dim * j + i ] = atof(val);
      }
    }

    for (i = 0 ; i < dim ; ++i)
    {
      fscanf(LCPfile , "%s" , val);
      q[i] = atof(val);
    }

    fscanf(LCPfile , "%s" , val);

    if (!feof(LCPfile))
    {

      sol[0] = atof(val);

      for (i = 1 ; i < dim ; ++i)
      {
        fscanf(LCPfile , "%s" , val);
        sol[i] = atof(val);
      }
    }
    else
    {
      isol = 0;
      for (i = 0 ; i < dim ; ++i) sol[i] = 0.0;
    }

    fclose(LCPfile);

    pfc_2D_series(dim , vecM , q);

    free(sol);
    free(vecM);
    free(q);
  }
}

/* ****************************************************** *
 *          READ MATRICES OF DATA DIRECTORY               *
 * ****************************************************** */

void test_data(void)
{

  FILE *f1, *f2;

  int nl, nc, n;

  double qi, Mij;
  double *q, *vec;

  char val[20];

  /* WARNING STATIC SIZE */

  n = 314;

  printf("\n\n GRANUL TEST \n");
  if ((f1 = fopen("DATA/M_bille84mu.dat", "r")) == NULL)
  {
    perror("fopen 1");
    exit(1);
  }

  vec = (double*)malloc(n * n * sizeof(double));
  q   = (double*)malloc(n * sizeof(double));

  while (!feof(f1))
  {

    fscanf(f1, "%d", &nl);
    fscanf(f1, "%d", &nc);
    fscanf(f1, "%s", val);
    Mij = atof(val);
    vec[(nc - 1)*n + (nl - 1) ] = Mij;

  }

  if ((f2 = fopen("DATA/q_bille84mu.dat", "r")) == NULL)
  {
    perror("fopen 2");
    exit(2);
  }

  while (!feof(f2))
  {
    fscanf(f2, "%d", &nl);
    fscanf(f2, "%s", val);
    qi = atof(val);
    q[nl - 1] = -qi;
  }

  fclose(f2);
  fclose(f1);

  pfc_2D_series(n , vec , q);

  free(vec);
  free(q);

}

int main(void)
{

  /****************************************************************/
#ifdef BAVARD
  printf("\n ********** BENCHMARK FOR PFC 2D SOLVERS ********** \n\n");
#endif
  /****************************************************************/

  test_data();
  test_matrix();

}
