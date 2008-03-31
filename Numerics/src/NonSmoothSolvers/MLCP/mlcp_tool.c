#include "mlcp_tool.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
/*
 *if sW2V[i]==0
 *  v[i] not null w2[i] null
 *else
 *  v[i] null and w2[i] not null
 */
void mlcp_buildM(int * zw, double * M, double * Mref, int n, int m)
{
  int col, i;
  int npm = n + m;
  double * Aux;
  double * AuxRef;
  memcpy(M, Mref, n * npm * sizeof(double));
  Aux = M + n * npm;
  AuxRef = Mref + n * npm;

  for (col = 0; col < m; col++)
  {
    if (zw[col] == 0)
    {
      memcpy(Aux, AuxRef, npm * sizeof(double));
    }
    else
    {
      for (i = 0; i < npm; i++) Aux[i] = 0;
      /*memcpy(Aux,sColNul,npm*sizeof(double));*/
      Aux[n + col] = -1;
      /*M[(n+col)*npm+col+n]=-1;*/
    }
    Aux = Aux + n + m;
    AuxRef = AuxRef + n + m;
  }

}

void   mlcp_fillSolution(double * z1, double * z2, double * w1, double * w2, int n, int m, int* zw, double * Q)
{
  int lin;
  for (lin = 0; lin < n; lin++)
  {
    z1[lin] = Q[lin];
    w1[lin] = 0;
  }

  for (lin = 0; lin < m; lin++)
  {
    if (zw[lin] == 0)
    {
      w2[lin] = 0;
      z2[lin] = Q[n + lin];
    }
    else
    {
      z2[lin] = 0;
      w2[lin] = Q[n + lin];
    }
  }
}
