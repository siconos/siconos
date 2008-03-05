/* Siconos-Numerics version 3.0.0, Copyright INRIA 2005-2008.
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "LA.h"
#include <time.h>
#include "NonSmoothDrivers.h"
#include "pfc_3D_Fischer_Burmeister.h"

/* Compute function G */
void Compute_G_FB(int m, double *G, double *y , double *C, double *x , double *b, double *Ip, double *IP, double *I3, double coef21 , double r, double coef)
{

  double *A, *B;
  double coef2;
  int i, j, mm;
  mm = m * m;

  A     = (double*)malloc(mm * sizeof(double));
  B     = (double*)malloc(mm * sizeof(double));

  for (i = 0 ; i < m ; ++i)
    for (j = 0 ; j < m ; ++j)
      A[j * m + i] = B[j * m + i] = 0.;


  coef2 = coef * coef;

  /* Ip = matrix_inv2(Ip`)*/
  /* IP = matrix_inv2(Ip`)*mup */
  /* I3 = matrix_inv2(Ip)*e3 */

  A[0 * m + 0] =  C[0 * 3 + 0] + (C[1 * 3 + 0] * IP[0] + C[2 * 3 + 0] * IP[1]);
  A[1 * m + 0] = -(C[1 * 3 + 0] * Ip[0 * 2 + 0] + C[2 * 3 + 0] * Ip[0 * 2 + 1]);
  A[2 * m + 0] = -(C[1 * 3 + 0] * Ip[1 * 2 + 0] + C[2 * 3 + 0] * Ip[1 * 2 + 1]);

  A[0 * m + 1] =  -(C[0 * 3 + 1] * Ip[0 * 2 + 0] + C[0 * 3 + 2] * Ip[0 * 2 + 1] + (C[1 * 3 + 1] * IP[0] + C[2 * 3 + 1] * IP[1]) * Ip[0 * 2 + 0] + (C[1 * 3 + 2] * IP[0] + C[2 * 3 + 2] * IP[1]) * Ip[0 * 2 + 1]);
  A[0 * m + 2] =  -(C[0 * 3 + 1] * Ip[1 * 2 + 0] + C[0 * 3 + 2] * Ip[1 * 2 + 1] + (C[1 * 3 + 1] * IP[0] + C[2 * 3 + 1] * IP[1]) * Ip[1 * 2 + 0] + (C[1 * 3 + 2] * IP[0] + C[2 * 3 + 2] * IP[1]) * Ip[1 * 2 + 1]);

  A[1 * m + 1] = Ip[0 * 2 + 0] * (C[1 * 3 + 1] * Ip[0 * 2 + 0] + C[2 * 3 + 1] * Ip[0 * 2 + 1]) + Ip[0 * 2 + 1] * (C[1 * 3 + 2] * Ip[0 * 2 + 0] + C[2 * 3 + 2] * Ip[0 * 2 + 1]);
  A[2 * m + 1] = Ip[0 * 2 + 0] * (C[1 * 3 + 1] * Ip[1 * 2 + 0] + C[2 * 3 + 1] * Ip[1 * 2 + 1]) + Ip[0 * 2 + 1] * (C[1 * 3 + 2] * Ip[1 * 2 + 0] + C[2 * 3 + 2] * Ip[1 * 2 + 1]);

  A[1 * m + 2] = Ip[1 * 2 + 0] * (C[1 * 3 + 1] * Ip[0 * 2 + 0] + C[2 * 3 + 1] * Ip[0 * 2 + 1]) + Ip[1 * 2 + 1] * (C[1 * 3 + 2] * Ip[0 * 2 + 0] + C[2 * 3 + 2] * Ip[0 * 2 + 1]);
  A[2 * m + 2] = Ip[1 * 2 + 0] * (C[1 * 3 + 1] * Ip[1 * 2 + 0] + C[2 * 3 + 1] * Ip[1 * 2 + 1]) + Ip[1 * 2 + 1] * (C[1 * 3 + 2] * Ip[1 * 2 + 0] + C[2 * 3 + 2] * Ip[1 * 2 + 1]);

  A[3 * m + 1] = -I3[0];
  A[3 * m + 2] = -I3[1];

  A[0 * m + 3] = coef * (1 - I3[0] - I3[1]);
  A[1 * m + 3] = I3[0];
  A[2 * m + 3] = I3[1];

  B[0 * m + 1] = -2 * ((Ip[0 * 2 + 0] * Ip[0 * 2 + 0] + Ip[0 * 2 + 1] * Ip[0 * 2 + 1]) + (Ip[0 * 2 + 0] * Ip[1 * 2 + 0] + Ip[0 * 2 + 1] * Ip[1 * 2 + 1])) * coef * y[4];
  B[0 * m + 2] = -2 * ((Ip[1 * 2 + 0] * Ip[0 * 2 + 0] + Ip[1 * 2 + 1] * Ip[0 * 2 + 1]) + (Ip[1 * 2 + 0] * Ip[1 * 2 + 0] + Ip[1 * 2 + 1] * Ip[1 * 2 + 1])) * coef * y[4];

  B[1 * m + 1] =  2 * (Ip[0 * 2 + 0] * Ip[0 * 2 + 0] + Ip[0 * 2 + 1] * Ip[0 * 2 + 1]) * y[4];
  B[1 * m + 2] =  2 * (Ip[1 * 2 + 0] * Ip[0 * 2 + 0] + Ip[1 * 2 + 1] * Ip[0 * 2 + 1]) * y[4];
  B[2 * m + 1] =  2 * (Ip[1 * 2 + 0] * Ip[0 * 2 + 0] + Ip[1 * 2 + 1] * Ip[0 * 2 + 1]) * y[4];
  B[2 * m + 2] =  2 * (Ip[1 * 2 + 0] * Ip[1 * 2 + 0] + Ip[1 * 2 + 1] * Ip[1 * 2 + 1]) * y[4];

  B[4 * m + 1] = -2 * ((Ip[0 * 2 + 0] * Ip[0 * 2 + 0] + Ip[0 * 2 + 1] * Ip[0 * 2 + 1]) * (coef * y[0] - y[1]) + (Ip[0 * 2 + 0] * Ip[1 * 2 + 0] + Ip[0 * 2 + 1] * Ip[1 * 2 + 1]) * (coef * y[0] - y[2]));
  B[4 * m + 2] = -2 * ((Ip[1 * 2 + 0] * Ip[0 * 2 + 0] + Ip[1 * 2 + 1] * Ip[0 * 2 + 1]) * (coef * y[0] - y[1]) + (Ip[1 * 2 + 0] * Ip[1 * 2 + 0] + Ip[1 * 2 + 1] * Ip[1 * 2 + 1]) * (coef * y[0] - y[2]));

  B[0 * m + 4] =  2 * coef2 * y[0] - 2 * coef * (((Ip[0 * 2 + 0] * Ip[0 * 2 + 0] + Ip[0 * 2 + 1] * Ip[0 * 2 + 1]) * (coef * y[0] - y[1]) + (Ip[0 * 2 + 0] * Ip[1 * 2 + 0] + Ip[0 * 2 + 1] * Ip[1 * 2 + 1]) * (coef * y[0] - y[2]))
                  + ((Ip[1 * 2 + 0] * Ip[0 * 2 + 0] + Ip[1 * 2 + 1] * Ip[0 * 2 + 1]) * (coef * y[0] - y[1]) + (Ip[1 * 2 + 0] * Ip[1 * 2 + 0] + Ip[1 * 2 + 1] * Ip[1 * 2 + 1]) * (coef * y[0] - y[2])));
  B[1 * m + 4] =  2 * ((Ip[0 * 2 + 0] * Ip[0 * 2 + 0] + Ip[0 * 2 + 1] * Ip[0 * 2 + 1]) * (coef * y[0] - y[1]) + (Ip[0 * 2 + 0] * Ip[1 * 2 + 0] + Ip[0 * 2 + 1] * Ip[1 * 2 + 1]) * (coef * y[0] - y[2]));
  B[2 * m + 4] =  2 * ((Ip[1 * 2 + 0] * Ip[0 * 2 + 0] + Ip[1 * 2 + 1] * Ip[0 * 2 + 1]) * (coef * y[0] - y[1]) + (Ip[1 * 2 + 0] * Ip[1 * 2 + 0] + Ip[1 * 2 + 1] * Ip[1 * 2 + 1]) * (coef * y[0] - y[2]));

  for (i = 0; i < m; i++)
  {
    x[i] = b[i];
    for (j = 0; j < m; j++)
      x[i] += (A[j * m + i] + B[j * m + i]) * y[j];
  }

  /*   for( i = 0 ; i < m ; ++i ) */
  /*     printf("b[%i] =  %14.7e\n",i,b[i]); printf("\n"); */

  /*    for( i = 0 ; i < m ; ++i ) */
  /*     printf("y[%i] =  %14.7e\n",i,y[i]); printf("\n"); */

  for (i = 0; i < m; i++)
    x[i] = (A[0 * m + i] + B[0 * m + i]) * y[0] + (A[1 * m + i] + B[1 * m + i]) * y[1] + (A[2 * m + i] + B[2 * m + i]) * y[2] + (A[3 * m + i] + B[3 * m + i]) * y[3] + (A[4 * m + i] + B[4 * m + i]) * y[4] + b[i];

  for (i = 0 ; i < m ; ++i)
  {
    G[i] = sqrt(x[i] * x[i] +  y[i] * y[i]) - (x[i] + y[i]);
  }


  /*  for( i = 0 ; i < m ; ++i ){ */
  /*     for( j = 0 ; j < m ; ++j ){ */
  /*       printf("A[%i,%i] =  %14.7e\t",i,j,A[j*m+i]); */
  /*     } */
  /*     printf("\n"); */
  /*   } */
  /*   printf("\n"); */
  /*   for( i = 0 ; i < m ; ++i ){ */
  /*     for( j = 0 ; j < m ; ++j ){ */
  /*       printf("B[%i,%i] =  %14.7e\t",i,j,B[j*m+i]); */
  /*     } */
  /*     printf("\n"); */
  /*   } */
  /*   printf("\n"); */

  /*   for( i = 0 ; i < m ; ++i ) */
  /*     printf("x[%i] =  %14.7e\n",i,x[i]); printf("\n"); */
  /*   for( i = 0 ; i < m ; ++i ) */
  /*     printf("G[%i] =  %14.7e\n",i,sqrt(x[i]*x[i] +  y[i]*y[i]) - (x[i] + y[i])); printf("\n"); */




  free(A);
  free(B);
}



/* Compute Jacobian of function G */
void Compute_JacG_FB(int m, double *JacG , double *y , double *C , double *x , double *b, double *Ip , double *IP, double *I3, double coef3 , double coef21, double coef)
{

  double *A, *B, *r, *delta;
  double coef2;
  int i , j , mm;
  mm = m * m;

  A     = (double*)malloc(mm * sizeof(double));
  B     = (double*)malloc(mm * sizeof(double));
  delta = (double*)malloc(mm * sizeof(double));
  r     = (double*)malloc(m * sizeof(double));

  for (i = 0 ; i < m ; ++i)
  {
    r[i] = 0.;
    for (j = 0 ; j < m ; ++j)
      A[j * m + i] = B[j * m + i] = delta[j * m + i] = 0.;
  }

  coef2 = coef * coef;

  /* Ip = matrix_inv2(Ip`)*/
  /* IP = matrix_inv2(Ip`)*mup */
  /* I3 = matrix_inv2(Ip)*e3 */

  A[0 * m + 0] =  C[0 * 3 + 0] + (C[1 * 3 + 0] * IP[0] + C[2 * 3 + 0] * IP[1]);
  A[1 * m + 0] = -(C[1 * 3 + 0] * Ip[0 * 2 + 0] + C[2 * 3 + 0] * Ip[0 * 2 + 1]);
  A[2 * m + 0] = -(C[1 * 3 + 0] * Ip[1 * 2 + 0] + C[2 * 3 + 0] * Ip[1 * 2 + 1]);

  A[0 * m + 1] =  -(C[0 * 3 + 1] * Ip[0 * 2 + 0] + C[0 * 3 + 2] * Ip[0 * 2 + 1] + (C[1 * 3 + 1] * IP[0] + C[2 * 3 + 1] * IP[1]) * Ip[0 * 2 + 0] + (C[1 * 3 + 2] * IP[0] + C[2 * 3 + 2] * IP[1]) * Ip[0 * 2 + 1]);
  A[0 * m + 2] =  -(C[0 * 3 + 1] * Ip[1 * 2 + 0] + C[0 * 3 + 2] * Ip[1 * 2 + 1] + (C[1 * 3 + 1] * IP[0] + C[2 * 3 + 1] * IP[1]) * Ip[1 * 2 + 0] + (C[1 * 3 + 2] * IP[0] + C[2 * 3 + 2] * IP[1]) * Ip[1 * 2 + 1]);

  A[1 * m + 1] = Ip[0 * 2 + 0] * (C[1 * 3 + 1] * Ip[0 * 2 + 0] + C[2 * 3 + 1] * Ip[0 * 2 + 1]) + Ip[0 * 2 + 1] * (C[1 * 3 + 2] * Ip[0 * 2 + 0] + C[2 * 3 + 2] * Ip[0 * 2 + 1]);
  A[2 * m + 1] = Ip[0 * 2 + 0] * (C[1 * 3 + 1] * Ip[1 * 2 + 0] + C[2 * 3 + 1] * Ip[1 * 2 + 1]) + Ip[0 * 2 + 1] * (C[1 * 3 + 2] * Ip[1 * 2 + 0] + C[2 * 3 + 2] * Ip[1 * 2 + 1]);

  A[1 * m + 2] = Ip[1 * 2 + 0] * (C[1 * 3 + 1] * Ip[0 * 2 + 0] + C[2 * 3 + 1] * Ip[0 * 2 + 1]) + Ip[1 * 2 + 1] * (C[1 * 3 + 2] * Ip[0 * 2 + 0] + C[2 * 3 + 2] * Ip[0 * 2 + 1]);
  A[2 * m + 2] = Ip[1 * 2 + 0] * (C[1 * 3 + 1] * Ip[1 * 2 + 0] + C[2 * 3 + 1] * Ip[1 * 2 + 1]) + Ip[1 * 2 + 1] * (C[1 * 3 + 2] * Ip[1 * 2 + 0] + C[2 * 3 + 2] * Ip[1 * 2 + 1]);

  A[3 * m + 1] = -I3[0];
  A[3 * m + 2] = -I3[1];

  A[0 * m + 3] = coef * (1 - I3[0] - I3[1]);
  A[1 * m + 3] = I3[0];
  A[2 * m + 3] = I3[1];

  B[0 * m + 1] = -2 * ((Ip[0 * 2 + 0] * Ip[0 * 2 + 0] + Ip[0 * 2 + 1] * Ip[0 * 2 + 1]) + (Ip[0 * 2 + 0] * Ip[1 * 2 + 0] + Ip[0 * 2 + 1] * Ip[1 * 2 + 1])) * coef * y[4];
  B[0 * m + 2] = -2 * ((Ip[1 * 2 + 0] * Ip[0 * 2 + 0] + Ip[1 * 2 + 1] * Ip[0 * 2 + 1]) + (Ip[1 * 2 + 0] * Ip[1 * 2 + 0] + Ip[1 * 2 + 1] * Ip[1 * 2 + 1])) * coef * y[4];

  B[1 * m + 1] =  2 * (Ip[0 * 2 + 0] * Ip[0 * 2 + 0] + Ip[0 * 2 + 1] * Ip[0 * 2 + 1]) * y[4];
  B[1 * m + 2] =  2 * (Ip[1 * 2 + 0] * Ip[0 * 2 + 0] + Ip[1 * 2 + 1] * Ip[0 * 2 + 1]) * y[4];
  B[2 * m + 1] =  2 * (Ip[0 * 2 + 0] * Ip[1 * 2 + 0] + Ip[0 * 2 + 1] * Ip[1 * 2 + 1]) * y[4];
  B[2 * m + 2] =  2 * (Ip[1 * 2 + 0] * Ip[1 * 2 + 0] + Ip[1 * 2 + 1] * Ip[1 * 2 + 1]) * y[4];

  B[4 * m + 1] = -2 * ((Ip[0 * 2 + 0] * Ip[0 * 2 + 0] + Ip[0 * 2 + 1] * Ip[0 * 2 + 1]) * (coef * y[0] - y[1]) + (Ip[0 * 2 + 0] * Ip[1 * 2 + 0] + Ip[0 * 2 + 1] * Ip[1 * 2 + 1]) * (coef * y[0] - y[2]));
  B[4 * m + 2] = -2 * ((Ip[1 * 2 + 0] * Ip[0 * 2 + 0] + Ip[1 * 2 + 1] * Ip[0 * 2 + 1]) * (coef * y[0] - y[1]) + (Ip[1 * 2 + 0] * Ip[1 * 2 + 0] + Ip[1 * 2 + 1] * Ip[1 * 2 + 1]) * (coef * y[0] - y[2]));

  B[0 * m + 4] =  2 * coef2 * y[0] - 2 * coef * (((Ip[0 * 2 + 0] * Ip[0 * 2 + 0] + Ip[0 * 2 + 1] * Ip[0 * 2 + 1]) * (coef * y[0] - y[1]) + (Ip[0 * 2 + 0] * Ip[1 * 2 + 0] + Ip[0 * 2 + 1] * Ip[1 * 2 + 1]) * (coef * y[0] - y[2]))
                  + ((Ip[1 * 2 + 0] * Ip[0 * 2 + 0] + Ip[1 * 2 + 1] * Ip[0 * 2 + 1]) * (coef * y[0] - y[1]) + (Ip[1 * 2 + 0] * Ip[1 * 2 + 0] + Ip[1 * 2 + 1] * Ip[1 * 2 + 1]) * (coef * y[0] - y[2])));
  B[1 * m + 4] =  2 * ((Ip[0 * 2 + 0] * Ip[0 * 2 + 0] + Ip[0 * 2 + 1] * Ip[0 * 2 + 1]) * (coef * y[0] - y[1]) + (Ip[0 * 2 + 0] * Ip[1 * 2 + 0] + Ip[0 * 2 + 1] * Ip[1 * 2 + 1]) * (coef * y[0] - y[2]));
  B[2 * m + 4] =  2 * ((Ip[1 * 2 + 0] * Ip[0 * 2 + 0] + Ip[1 * 2 + 1] * Ip[0 * 2 + 1]) * (coef * y[0] - y[1]) + (Ip[1 * 2 + 0] * Ip[1 * 2 + 0] + Ip[1 * 2 + 1] * Ip[1 * 2 + 1]) * (coef * y[0] - y[2]));

  /* for (i=0;i<m;i++){ */
  /*     x[i] = b[i]; */
  /*     for (j=0;j<m;j++) */
  /*       x[i] += (A[j*m+i]+ B[j*m+i])*y[j]; */
  /*   } */

  for (i = 0; i < m; i++)
    x[i] = (A[0 * m + i] + B[0 * m + i]) * y[0] + (A[1 * m + i] + B[1 * m + i]) * y[1] + (A[2 * m + i] + B[2 * m + i]) * y[2] + (A[3 * m + i] + B[3 * m + i]) * y[3] + (A[4 * m + i] + B[4 * m + i]) * y[4] + b[i];

  for (i = 0; i < m; i++)
    delta[i * m + i] = 1.;


  for (i = 0; i < m; i++)
  {
    r[i] = sqrt(x[i] * x[i] +  y[i] * y[i]);
    for (j = 0; j < m; j++)
    {
      if (r[i])
        JacG[j * m + i] = ((A[j * m + i] + B[j * m + i]) * x[i] + delta[j * m + i] * y[i]) / r[i] - (A[j * m + i] + B[j * m + i] + delta[j * m + i]);
      else
        JacG[j * m + i] =  - (A[j * m + i] + B[j * m + i] + delta[j * m + i]);
    }
  }

  /*    for( i = 0 ; i < m ; ++i ) */
  /*     printf("b[%i] =  %14.7e\n",i,b[i]); printf("\n"); */

  /*    for( i = 0 ; i < m ; ++i ) */
  /*     printf("y[%i] =  %14.7e\n",i,y[i]); printf("\n"); */

  /*   for( i = 0 ; i < m ; ++i ){ */
  /*     for( j = 0 ; j < m ; ++j ){ */
  /*       printf("A[%i,%i] =  %14.7e\t",i,j,A[j*m+i]); */
  /*     } */
  /*     printf("\n"); */
  /*   } */
  /*   printf("\n"); */
  /*   for( i = 0 ; i < m ; ++i ){ */
  /*     for( j = 0 ; j < m ; ++j ){ */
  /*       printf("B[%i,%i] =  %14.7e\t",i,j,B[j*m+i]); */
  /*     } */
  /*     printf("\n"); */
  /*   } */
  /*   printf("\n"); */

  /*   for( i = 0 ; i < m ; ++i ) */
  /*     printf("x[%i] =  %14.7e\n",i,x[i]); printf("\n"); */
  /*   for( i = 0 ; i < m ; ++i ) */
  /*     printf("G[%i] =  %14.7e\n",i,sqrt(x[i]*x[i] +  y[i]*y[i]) - (x[i] + y[i])); printf("\n"); */



  /*  for( i = 0 ; i < m ; ++i ){ */
  /*     //    printf("r[%i] =  %14.7e\n",i,r[i]); printf("\n"); */
  /*     for( j = 0 ; j < m ; ++j ){ */
  /*       printf("JacG[%i,%i] =  %14.7e\t",i,j,JacG[j*m+i]); */
  /*     } */
  /*     printf("\n"); */
  /*   } */
  /*   printf("\n"); */


  free(A);
  free(B);
  free(r);
  free(delta);
}

/*  Inverse Matrix 2x2  */
void matrix_inv2(double *a, double *b)
{
  double det;
  det = a[0 * 2 + 0] * a[1 * 2 + 1] - a[1 * 2 + 0] * a[0 * 2 + 1];

  b[0 * 2 + 0] =  a[1 * 2 + 1] / det;
  b[0 * 2 + 1] = -a[0 * 2 + 1] / det;
  b[1 * 2 + 0] = -a[1 * 2 + 0] / det;
  b[1 * 2 + 1] =  a[0 * 2 + 0] / det;

}

/* Lineserach */
void Linesearch_FB(int n, double *zz, double *ww, double *www, double *b, double *C, double *Ip, double *IP, double *I3, double coef21, double coef3, double mu, double err1)
{

  double alpha, err2, qs, a1;
  double *G, *zzzz, *wwww;

  int i, incx, incy;

  zzzz = (double*)malloc(n * sizeof(double));
  wwww = (double*)malloc(n * sizeof(double));
  G    = (double*)malloc(n * sizeof(double));

  for (i = 0 ; i < n ; ++i)
    zzzz[i] = wwww[i] = G[i] = 0.;

  incx =  1;
  incy =  1;
  a1 = -1.;
  qs = 1.;

  alpha = 1.;
  while (alpha > 0.05)
  {

    DCOPY(n , zz , incx , zzzz , incy);
    DAXPY(n , alpha , www , incx , zzzz , incy);

    Compute_G_FB(n , G , zzzz , C , wwww , b , Ip , IP , I3 , coef21 , coef3 , mu);
    err2 = DNRM2(n, G , incx);

    if (err2 < err1) break;
    alpha = alpha * 0.5;
  }
  err1 = err2;

  DCOPY(n , zzzz , incx , zz , incy);
  DCOPY(n , wwww , incx , ww , incy);

  free(zzzz);
  free(wwww);
  free(G);

}
