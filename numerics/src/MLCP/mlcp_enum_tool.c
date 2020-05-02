/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

#include "mlcp_enum_tool.h"
#include <stdio.h>
#include <string.h>
#include "numerics_verbose.h"

static unsigned long long int sCurrentEnum = 0;
static unsigned long long int sCmpEnum = 0;
static unsigned long long int sNbCase = 0;
static double sProgress = 0;

unsigned long long int computeNbCase(int M)
{
  unsigned long long int nbCase = 1;
  for(int cmp = 0; cmp < M; cmp++)
    nbCase = nbCase << 1;
  return nbCase;
}

void initEnum(int M)
{
  /*  sCurrentEnum = 0;*/
  numerics_printf_verbose(1,"----- initEnum -- problem size :%i", M);
  numerics_printf_verbose(1,"----- initEnum -- currentEnum :%i", (int)sCurrentEnum);

  sCmpEnum = 0;

  sNbCase  = computeNbCase(M);

  sProgress = 0;
  numerics_printf_verbose(1,"----- initEnum -- number of cases :%i", (int)sNbCase);
}

static void affectW2V(int * zw, int size)
{
  unsigned long  int aux = sCurrentEnum;
  for(int i = 0; i < size; i++)
  {
    zw[i] = aux & 1;
    aux = aux >> 1;
  }

  if(verbose > 1)
  {
    for(int i = 0; i < size; i++)
      printf("zw[%d]=%d \t", i, zw[i]);
    printf("\n");
  }

}

int nextEnum(int * W2V, int size)
{
  if(sCmpEnum == sNbCase)
    return 0;
  if(sCurrentEnum >= sNbCase)
    sCurrentEnum = 0;

  numerics_printf_verbose(1,"----- nextEnum -- try enum :%d", (int)sCurrentEnum);

  affectW2V(W2V, size);
  sCurrentEnum++;
  sCmpEnum++;
  if(verbose && sCmpEnum > (unsigned long int)sProgress * sNbCase)
  {
    sProgress += 0.001;
    numerics_printf_verbose(1,"progress %f %d / %d", sProgress, (int) sCurrentEnum,  sNbCase);
  }

  return 1;
}
/*
 *if zw[i]==0
 *  v[i] >=0 and w_i[i]=0
 *else
 *  v[i] =0 and w_i[i] >=0
 */

void mlcp_enum_build_M(int * zw, double * M, double * Mref, int n, int m, int NbLines)
{

  /*First, copy the n first collums.*/
  memcpy(M, Mref, n * NbLines * sizeof(double));

  double * current_col_M = M + n * NbLines;
  double * current_col_MRef = Mref + n * NbLines;

  for(int col = 0; col < m; col++)
  {
    if(zw[col] == 0)
    {
      /* copy current column of MRef in M */
      memcpy(current_col_M, current_col_MRef, NbLines * sizeof(double));
    }
    else
    {
      for(int i = 0; i < NbLines; i++) current_col_M[i] = 0;
      /*memcpy(current_col_M, sColNul, npm*sizeof(double));*/
      current_col_M[(NbLines - m) + col] = -1;
      /*M[(n+col)*npm+col+n]=-1;*/
    }
    current_col_M = current_col_M + NbLines;
    current_col_MRef = current_col_MRef + NbLines;
  }

}
void mlcp_enum_build_M_Block(int * zw, double * M, double * Mref, int n, int m, int NbLines, int *indexInBlock)
{
  int col, i;
  double * current_col_M;
  memcpy(M, Mref, NbLines * NbLines * sizeof(double));
  //current_col_M=M+n*NbLines;
  //current_col_MRef = Mref+n*NbLines;

  for(col = 0; col < m; col++)
  {
    if(zw[col])
    {
      current_col_M = M + indexInBlock[col] * NbLines;
      for(i = 0; i < NbLines; i++) current_col_M[i] = 0;
      current_col_M[indexInBlock[col]] = -1;
    }
  }
}

void   mlcp_enum_fill_solution(double * z1, double * z2, double * w1, double * w2, int n, int m, int NbLines, int* zw, double * Q)
{
  int lin;
  for(lin = 0; lin < n; lin++)
  {
    z1[lin] = Q[lin];
  }
  for(lin = 0; lin < NbLines; lin++)
    w1[lin] = 0;

  for(lin = 0; lin < m; lin++)
  {
    if(zw[lin] == 0)
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

void   mlcp_enum_fill_solution_Block(double * z, double * w, int n, int m, int NbLines, int* zw, double * Q, int *indexInBlock)
{
  int lin;
  for(lin = 0; lin < NbLines; lin++)
  {
    z[lin] = Q[lin];
  }
  for(lin = 0; lin < NbLines; lin++)
    w[lin] = 0;

  for(lin = 0; lin < m; lin++)
  {
    if(zw[lin])
    {
      z[indexInBlock[lin]] = 0;
      w[indexInBlock[lin]] = Q[indexInBlock[lin]];
    }
  }
}

void   mlcp_enum_display_solution(double * z1, double * z2, double * w1, double * w2, int n, int m, int Nblines)
{
  int lin;
  printf("z1:\n");
  for(lin = 0; lin < n; lin++)
    printf("z1[,%d]=%.15e\n", lin, z1[lin]);

  for(lin = 0; lin < Nblines - m; lin++)
    printf("w1[%d]=%.15e\n", lin, w1[lin]);

  printf("z2,w2:\n");
  for(lin = 0; lin < m; lin++)
    printf("z2[%d],w2[%d],=%.15e\t%.15e\n", lin, lin, z2[lin], w2[lin]);


}

void   mlcp_enum_display_solution_Block(double * z, double * w, int n, int m, int Nblines, int *indexInBlock)
{
  int lin;
  int curCompIndex = 0;

  for(lin = 0; lin < Nblines; lin++)
    if(indexInBlock && indexInBlock[curCompIndex] == lin)
    {
      printf("z[%d],w[%d],=%.15e\t%.15e. (complementarity cond)\n", lin, lin, z[lin], w[lin]);
      curCompIndex++;
    }
    else
    {
      printf("z[%d],w[%d],=%.15e\t%.15e. (equality cond)\n", lin, lin, z[lin], w[lin]);
    }
}


void mlcp_enum_build_indexInBlock(MixedLinearComplementarityProblem* problem, int *indexInBlock)
{
  int numBlock = 0;
  int n = problem->n; /* Equalities */
  int m = problem->m; /* Inequalities */
  int curCompl = 0;
  int curComplLine = 0;
  int lin;
  while(problem->blocksRows[numBlock] < n + m)
  {
    if(!problem->blocksIsComp[numBlock])
    {
      for(lin = problem->blocksRows[numBlock]; lin < problem->blocksRows[numBlock + 1]; lin++)
      {
        curComplLine++;
      }
    }
    else
    {
      for(lin = problem->blocksRows[numBlock]; lin < problem->blocksRows[numBlock + 1]; lin++)
      {
        indexInBlock[curCompl] = curComplLine;
        curCompl++;
        curComplLine++;
      }
    }
    numBlock++;
  }
}
