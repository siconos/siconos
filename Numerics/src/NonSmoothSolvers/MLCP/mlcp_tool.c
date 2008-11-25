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
void mlcp_buildM(int * zw, double * M, double * Mref, int n, int m, int NbLines)
{
  int col, i;
  int npm = n + m;
  double * Aux;
  double * AuxRef;
  /*First, copy the n first collums.*/
  memcpy(M, Mref, n * NbLines * sizeof(double));
  Aux = M + n * NbLines;
  AuxRef = Mref + n * NbLines;

  for (col = 0; col < m; col++)
  {
    if (zw[col] == 0)
    {
      memcpy(Aux, AuxRef, NbLines * sizeof(double));
    }
    else
    {
      for (i = 0; i < NbLines; i++) Aux[i] = 0;
      /*memcpy(Aux,sColNul,npm*sizeof(double));*/
      Aux[(NbLines - m) + col] = -1;
      /*M[(n+col)*npm+col+n]=-1;*/
    }
    Aux = Aux + NbLines;
    AuxRef = AuxRef + NbLines;
  }

}

void   mlcp_fillSolution(double * z1, double * z2, double * w1, double * w2, int n, int m, int NbLines, int* zw, double * Q)
{
  int lin;
  for (lin = 0; lin < n; lin++)
  {
    z1[lin] = Q[lin];
  }
  for (lin = 0; lin < NbLines; lin++)
    w1[lin] = 0;

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

void   mlcp_DisplaySolution(double * z1, double * z2, double * w1, double * w2, int n, int m, int Nblines)
{
  int lin;
  printf("z1:\n");
  for (lin = 0; lin < n; lin++)
    printf("z1[,%d]=%e\n", lin, z1[lin]);

  for (lin = 0; lin < Nblines - m; lin++)
    printf("w1[%d]=%e\n", lin, w1[lin]);

  printf("w2,z2:\n");
  for (lin = 0; lin < m; lin++)
    printf("z2[%d],w2[%d],=%e\t%e\n", lin, lin, z2[lin], w2[lin]);


}
