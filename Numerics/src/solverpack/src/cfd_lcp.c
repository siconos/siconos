
/*!\file cfd_lcp.c

   This file allows the formulation in the LCP (Linear  Complementary Problem) form of a contact problem with friction.

*/



/*!\fn void cfd_lcp (int *dim_n,double *mumu,double *vecc,double *qc,double *MM,double *q)

   cfd_lcp subroutine allows the formulation in the LCP (Linear  Complementary Problem) form of a contact problem with friction..

   \param dim_n On enter  a pointer over integers, the number of normal variables after the condensation of a contact friction problem.
   \param mumu On enter a pointer over doubles, the friction coefficient.
   \param vecc On enter a pointer over doubles containing the components of a double matrix (2*dim_n,2*dim_n) with a fortran90 allocation.
   \param qc On enter a pointer over doubles containing the components of a double vector (2*dim_n).
   \param MM On return a pointer over doubles containing the components of a double matrix (3*dim_n,3*dim_n) with a fortran90 allocation.
   \param q On return a pointer over doubles, a double vector (3*dim_n).

   \author Nineb Sheherazade.

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


void cfd_lcp(int *dim_n, double *mumu, double *vecc, double *qc, double *MM, double *q)
{
  FILE *f10, *f14;
  int n = *dim_n, i, j;
  double mu = *mumu;
  double Kc[2 * n][2 * n], M[3 * n][3 * n];


  for (i = 0; i < 2 * n; i++)
    for (j = 0; j < 2 * n; j++)
      Kc[i][j] = vecc[i * 2 * n + j];

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      M[i][j] = Kc[i][j];
      M[i][n + j] = -Kc[i][n + j];
      M[i][2 * n + j] = Kc[i][n + j];
      M[n + i][j] = -Kc[n + i][j] + mu * Kc[i][j];
      M[n + i][n + j] = Kc[n + i][n + j] - mu * Kc[i][n + j];
      M[n + i][2 * n + j] = mu * Kc[i][n + j] - Kc[n + i][n + j];
      M[2 * n + i][j] = Kc[n + i][j] + mu * Kc[i][j];
      M[2 * n + i][n + j] = -Kc[n + i][n + j] - mu * Kc[i][n + j];
      M[2 * n + i][2 * n + j] = mu * Kc[i][n + j] + Kc[n + i][n + j];
    }

  for (i = 0; i < 3 * n; i++)
    for (j = 0; j < 3 * n; j++)
      MM[i * 3 * n + j] = 0.;

  for (i = 0; i < 3 * n; i++)
    for (j = 0; j < 3 * n; j++)
      MM[i + 3 * n * j] = M[i][j]; //fortran compatibility


  for (i = 0; i < n; i++)
  {
    q[i] = qc[i];
    q[n + i] = mu * qc[i] - qc[n + i];
    q[2 * n + i] = mu * qc[i] + qc[n + i];
  }

}
