/* Siconos-Numerics, Copyright INRIA 2005-2011.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
#ifndef MLCP_PROBLEM_C
#define MLCP_PROBLEM_C


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "MixedLinearComplementarityProblem.h"


void mixedLinearComplementarity_display(MixedLinearComplementarityProblem* p)
{
  int n = p->n;
  int m = p->m;
  printf("MLCP DISPLAY:\n-------------\n");
  printf("n :%d m: %d\n", p->n, p->m);


  printf(p->isStorageType1 ? "using (M)\n" : "not using (M)\n");
  printf(p->isStorageType2 ? "using (ABCD)\n" : "not using (ABCD)\n");
  if (p->blocksRows)
  {
    printf("blocks are:\n");
    int NumBlock = 0;
    while (p->blocksRows[NumBlock] < n + m)
    {
      if (p->blocksIsComp[NumBlock])
      {
        printf("->block of complementarity condition (type %d), from line %d, to line %d.\n", p->blocksIsComp[NumBlock], p->blocksRows[NumBlock], p->blocksRows[NumBlock + 1] - 1);
      }
      else
      {
        printf("->block of equality type (type %d), from line %d, to line %d.\n", p->blocksIsComp[NumBlock], p->blocksRows[NumBlock], p->blocksRows[NumBlock + 1] - 1);
      }
      NumBlock++;
    }
  }

  if (p->M)
  {
    printf("M matrix:\n");
    display(p->M);
  }
  else
    printf("No M matrix:\n");

  if (p->q)
  {
    printf("q matrix:\n");
    displayMat(p->q, n + m, 1, 0);
  }
  else
    printf("No q matrix:\n");

  if (p->A)
  {
    printf("A matrix:\n");
    displayMat(p->A, n, n, 0);
  }
  else
  {
    printf("No A matrix:\n");
    if (!p->M->storageType)
    {
      printf("A matrix from M:\n");
      displayMat(p->M->matrix0, n, n, n + m);
    }
  }
  if (p->B)
  {
    printf("B matrix:\n");
    displayMat(p->B, m, m, 0);
  }
  else
  {
    printf("No B matrix:\n");
    if (!p->M->storageType)
    {
      printf("B matrix from M:\n");
      displayMat(p->M->matrix0 + n * (n + m) + n, m, m, n + m);
    }
  }

  if (p->C)
  {
    printf("C matrix:\n");
    displayMat(p->C, n, m, 0);
  }
  else
  {
    printf("No C matrix:\n");
    if (!p->M->storageType)
    {
      printf("C matrix from M:\n");
      displayMat(p->M->matrix0 + n * (n + m), n, m, n + m);
    }
  }

  if (p->D)
  {
    printf("D matrix:\n");
    displayMat(p->D, m, n, 0);
  }
  else
  {
    printf("No D matrix:\n");
    if (!p->M->storageType)
    {
      printf("D matrix from M:\n");
      displayMat(p->M->matrix0 + n, m, n, n + m);
    }
  }
  if (p->a)
  {
    printf("a matrix:\n");
    displayMat(p->a, n, 1, 0);
  }
  else
    printf("No a matrix:\n");
  if (p->b)
  {
    printf("b matrix:\n");
    displayMat(p->b, m, 1, 0);
  }
  else
    printf("No b matrix:\n");

}
int mixedLinearComplementarity_printInFile(MixedLinearComplementarityProblem* problem, FILE* file)
{
  int info = 0;
  if (! problem)
  {
    fprintf(stderr, "Numerics, MixedLinearComplementarityProblem printInFile failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  int i, j;
  fprintf(file, "%d\n", problem->isStorageType1);
  fprintf(file, "%d\n", problem->isStorageType2);
  int n = problem->n;
  fprintf(file, "%d\n", n);
  int m = problem->m;
  fprintf(file, "%d\n", m);

  if (problem->isStorageType1)
  {
    int nbBlocks = 0;
    fprintf(file, "%d ", problem->blocksRows[nbBlocks]);
    while (problem->blocksRows[nbBlocks] < (n + m))
    {
      nbBlocks++;
      fprintf(file, "%d ", problem->blocksRows[nbBlocks]);
    }
    fprintf(file, "\n");
    for (i = 0; i < nbBlocks; i++)
    {
      fprintf(file, "%d ", problem->blocksIsComp[i]);
    }
    fprintf(file, "\n");
    printInFile(problem->M, file);

    for (i = 0; i < problem->M->size1; i++)
    {
      fprintf(file, "%32.24e ", problem->q[i]);
    }
    fprintf(file, "\n");
    /* return 1; */
    /* if (problem->isStorageType2)  */
    /* { */
    /*   printf("Numerics, MixedLinearComplementarityProblem printInFile only Storage1 has been printed.\n"); */
    /* } */
  }
  if (problem->isStorageType2)
  {

    for (i = 0; i < problem->n; i++)
    {
      for (j = 0; j < problem->n; j++)
      {
        fprintf(file, "%32.24e ", problem->A[i + j * n]);
      }
      fprintf(file, "\n");
    }

    for (i = 0; i < problem->m; i++)
    {
      for (j = 0; j < problem->m; j++)
      {
        fprintf(file, "%32.24e ", problem->B[i + j * m]);
      }
      fprintf(file, "\n");
    }
    for (i = 0; i < problem->n; i++)
    {
      for (j = 0; j < problem->m; j++)
      {
        fprintf(file, "%32.24e ", problem->C[i + j * n]);
      }
      fprintf(file, "\n");
    }
    for (i = 0; i < problem->m; i++)
    {
      for (j = 0; j < problem->n; j++)
      {
        fprintf(file, "%32.24e ", problem->D[i + j * m]);
      }
      fprintf(file, "\n");
    }

    for (i = 0; i < problem->n; i++)
    {
      fprintf(file, "%32.24e ", problem->a[i]);
    }
    fprintf(file, "\n");
    for (i = 0; i < problem->m; i++)
    {
      fprintf(file, "%32.24e ", problem->b[i]);
    }

  }
  return info;

}

int mixedLinearComplementarity_newFromFile(MixedLinearComplementarityProblem* problem, FILE* file)
{
  int info = 0;
  if (! problem)
  {
    fprintf(stderr, "Numerics, MixedLinearComplementarityProblem printInFile failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  int i, j;
  int st1, st2;
  fscanf(file, "%d\n", &st1);
  problem->isStorageType1 = st1;
  fscanf(file, "%d\n", &st2);
  problem->isStorageType2 = st2;

  int n ;
  fscanf(file, "%d\n", &n);
  problem->n = n;
  int m;
  fscanf(file, "%d\n", &m);
  problem->m = m;



  if (problem->isStorageType1)
  {
    int * blocksRows = (int *)malloc((n + m + 1) * sizeof(int));
    int nbBlocks = 0;
    fscanf(file, "%d ", &(blocksRows[nbBlocks]));
    while (blocksRows[nbBlocks] < (n + m))
    {
      nbBlocks++;
      fscanf(file, "%d ", &(blocksRows[nbBlocks]));
    }
    problem->blocksRows = (int *)malloc((nbBlocks + 1) * sizeof(int));
    //fscanf(file,"\n");
    for (i = 0; i <= nbBlocks; i++)
    {
      problem->blocksRows[i] = blocksRows[i];
    }
    free(blocksRows);
    problem->blocksIsComp = (int *)malloc((nbBlocks) * sizeof(int));

    //fprintf(file,"\n");
    for (i = 0; i < nbBlocks; i++)
    {
      fscanf(file, "%d ", &(problem->blocksIsComp[i]));
    }

    //fprintf(file,"\n");
    problem->M = (NumericsMatrix *)malloc(sizeof(NumericsMatrix));
    newFromFile(problem->M, file);
    problem->q = (double *) malloc(problem->M->size1 * sizeof(double));

    for (i = 0; i < problem->M->size1; i++)
    {
      fscanf(file, "%lf ", &(problem->q[i]));
    }
    //fprintf(file,"\n");
    /* return 1; */
    /* if (problem->isStorageType2)  */
    /* { */
    /*   printf("Numerics, MixedLinearComplementarityProblem printInFile only Storage1 has been printed.\n"); */
    /* } */
  }
  if (problem->isStorageType2)
  {

    problem->A = (double*)malloc(n * n * sizeof(double));
    problem->B = (double*)malloc(m * m * sizeof(double));
    problem->C = (double*)malloc(n * m * sizeof(double));
    problem->D = (double*)malloc(m * n * sizeof(double));
    problem->a = (double*)malloc(n * sizeof(double));
    problem->b = (double*)malloc(m * sizeof(double));

    for (i = 0; i < problem->n; i++)
    {
      for (j = 0; j < problem->n; j++)
      {
        fscanf(file, "%lf ", &(problem->A[i + j * n]));
      }
      /* fscanf(file,"\n"); */
    }

    for (i = 0; i < problem->m; i++)
    {
      for (j = 0; j < problem->m; j++)
      {
        fscanf(file, "%lf ", &(problem->B[i + j * m]));
      }
      /* fprintf(file,"\n"); */
    }
    for (i = 0; i < problem->n; i++)
    {
      for (j = 0; j < problem->m; j++)
      {
        fscanf(file, "%lf ", &(problem->C[i + j * n]));
      }
      /* fprintf(file,"\n"); */
    }
    for (i = 0; i < problem->m; i++)
    {
      for (j = 0; j < problem->n; j++)
      {
        fscanf(file, "%lf ", &(problem->D[i + j * m]));
      }
      /* fprintf(file,"\n"); */
    }

    for (i = 0; i < problem->n; i++)
    {
      fscanf(file, "%lf ", &(problem->a[i]));
    }
    /* fprintf(file,"\n"); */
    for (i = 0; i < problem->m; i++)
    {
      fscanf(file, "%lf ", &(problem->b[i]));
    }

  }
  return info;

}


int mixedLinearComplementarity_newFromFileOld(MixedLinearComplementarityProblem* problem, FILE* MLCPfile)
{
  int n = 0, m = 0, NbLines = 0;
  int i, j,  n2, m2, isol;
  char val[128];

  double *vecA, *vecB, *vecC, *vecD, *vecM, *vecQ;
  double *a, *b;
  int nread;

  nread = fscanf(MLCPfile , "%d" , &n);
  nread = fscanf(MLCPfile , "%d" , &m);
  nread = fscanf(MLCPfile , "%d" , &NbLines);

  n2 = n * n;
  m2 = m * m;
  isol = 1;

  vecM = (double*)malloc((n + m) * (NbLines) * sizeof(double));

  vecQ = (double*)malloc((NbLines) * sizeof(double));
  vecA = (double*)malloc(n * (NbLines - m) * sizeof(double));
  vecB = (double*)malloc(m2 * sizeof(double));
  vecC = (double*)malloc((NbLines - m) * m * sizeof(double));
  vecD = (double*)malloc(m * n * sizeof(double));
  a    = (double*)malloc((NbLines - m) * sizeof(double));
  b    = (double*)malloc(m * sizeof(double));

  problem->blocksRows = (int*)malloc(3 * sizeof(int));
  problem->blocksIsComp = (int*)malloc(2 * sizeof(int));
  problem->blocksRows[0] = 0;
  problem->blocksRows[1] = n;
  problem->blocksRows[2] = n + m;
  problem->blocksIsComp[0] = 0;
  problem->blocksIsComp[1] = 1;





  problem->M = (NumericsMatrix *)malloc(sizeof(NumericsMatrix));


  NumericsMatrix * M = problem->M;

  M->storageType = 0;
  M->matrix0 = vecM;

  problem->isStorageType1 = 1; // Both problems seems to be stored
  problem->isStorageType2 = 1; // Both problems seems to be stored

  problem->q = vecQ;
  problem->A = vecA;
  problem->B = vecB;
  problem->C = vecC;
  problem->D = vecD;
  problem->a = a;
  problem->b = b;
  problem->blocksRows[1] = n;
  problem->blocksRows[2] = n + m;
  problem->n = n;
  problem->m = m;

  M->size0 = NbLines;
  M->size1 = n + m;



  for (i = 0 ; i < NbLines - m ; ++i)
  {
    for (j = 0 ; j < n ; ++j)
    {
      nread = fscanf(MLCPfile, "%s", val);
      vecA[(NbLines - m)*j + i ] = atof(val);
      vecM[(NbLines)*j + i ] = atof(val);
    }
  }
  for (i = 0 ; i < m ; ++i)
  {
    for (j = 0 ; j < m ; ++j)
    {
      nread = fscanf(MLCPfile, "%s", val);
      vecB[ m * j + i ] = atof(val);
      /*  vecM[ n*(m+n)+(n+m)*j+n+i ] = atof(val);*/
      vecM[ n * (NbLines) + (NbLines)*j + (NbLines - m) + i ] = atof(val);

    }
  }
  for (i = 0 ; i < NbLines - m ; ++i)
  {
    for (j = 0 ; j < m ; ++j)
    {
      nread = fscanf(MLCPfile, "%s", val);
      vecC[(NbLines - m)*j + i ] = atof(val);
      vecM[(NbLines) * (n + j) + i ] = atof(val);
    }
  }
  for (i = 0 ; i < m ; ++i)
  {
    for (j = 0 ; j < n ; ++j)
    {
      nread = fscanf(MLCPfile, "%s", val);
      vecD[ m * j + i ] = atof(val);
      vecM[(NbLines)*j + i + (NbLines - m) ] = atof(val);
    }
  }

  for (i = 0 ; i < NbLines - m ; ++i)
  {
    nread = fscanf(MLCPfile , "%s" , val);
    a[i] = atof(val);
    vecQ[i] = atof(val);
  }
  for (i = 0 ; i < m ; ++i)
  {
    nread = fscanf(MLCPfile , "%s" , val);
    b[i] = atof(val);
    vecQ[i + NbLines - m] = atof(val);
  }

  return 0;
}

int mixedLinearComplementarity_newFromFilename(MixedLinearComplementarityProblem* problem, char* filename)
{
  int info = 0;
  FILE * file = fopen(filename, "r");

  info = mixedLinearComplementarity_newFromFile(problem, file);

  fclose(file);
  return info;
}




void freeMixedLinearComplementarityProblem(MixedLinearComplementarityProblem* problem)
{
  if (problem->isStorageType1)
  {
    freeNumericsMatrix(problem->M);
    free(problem->M);
    free(problem->q);
    free(problem->blocksRows);
    free(problem->blocksIsComp);
  }
  if (problem->isStorageType2)
  {
    free(problem->A);
    free(problem->B);
    free(problem->C);
    free(problem->D);
    free(problem->a);
    free(problem->b);
  }
  free(problem);
  problem = NULL;
}
#endif
