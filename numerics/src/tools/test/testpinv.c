#include "pinv.h"
#include "cond.h"
#include <stdlib.h>
#include <stdio.h>
#include "NumericsMatrix.h"
#include <string.h>
int main(void)
{
  int n = 4;
  int m = 5;
  int info = -1;

  double * W = (double*)malloc(n * m * sizeof(double));
  double * Wpinv = (double*)malloc(n * m * sizeof(double));
  double * Wpinvtest = (double*)malloc(m * n * sizeof(double));

  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      W[i + j * n] = 0.0;
    }
  }
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < n; j++)
    {
      Wpinv[i + j * m] = 0.0;
      Wpinvtest[i + j * m] = 0.0;
    }
  }
  W[0 + 0 * n] = 1.0;
  W[0 + 4 * n] = 2.0;
  W[1 + 2 * n] = 3.0;
  W[3 + 1 * n] = 4.0;

  printf("Original Matrix W\n");
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      printf("%8.6e\t", W[i + j * n]) ;
    }
    printf("\n");
  }

  NumericsMatrix *Wnum = NM_new();
  Wnum->storageType = 0;
  Wnum-> size0 = n;
  Wnum-> size1 = m;
  Wnum->matrix1 = NULL;
  Wnum->matrix2 = NULL;
  Wnum->internalData = NULL;
  Wnum->matrix0 = W;

  FILE * file1 = fopen("dataW.dat", "w");
  NM_write_in_file_scilab(Wnum, file1);
  fclose(file1);

  NumericsMatrix *WnumpInv = NM_new();
  WnumpInv->storageType = 0;
  WnumpInv-> size0 = n;
  WnumpInv-> size1 = m;
  WnumpInv->matrix1 = NULL;
  WnumpInv->matrix2 = NULL;
  WnumpInv->internalData = NULL;
  WnumpInv->matrix0 = Wpinv;

  double tol = 1e-24;
  memcpy(Wpinv, W, n * m * sizeof(double));
  pinv(Wpinv, n, m, tol);
  printf("Winvtest\n");
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < n; j++)
    {
      printf("%8.6e\t", Wpinvtest[i + j * m]) ;
    }
    printf("\n");
  }

  Wpinvtest[0 + 0 * m] = 0.2;
  Wpinvtest[4 + 0 * m] = 0.4;
  Wpinvtest[2 + 1 * m] = 1.0 / 3.0;
  Wpinvtest[1 + 3 * m] = 1.0 / 4.0;

  printf("Winvtest\n");
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < n; j++)
    {
      printf("%8.6e\t", Wpinvtest[i + j * m]) ;
    }
    printf("\n");
  }
  printf("Pseudo inverseWinv\n");
  double err = 0.0;
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < n; j++)
    {
      printf("%8.6e\t", Wpinv[i + j * m]) ;
      err += (Wpinv[i + j * m] - Wpinvtest[i + j * m]) * (Wpinv[i + j * m] - Wpinvtest[i + j * m]);
    }
    printf("\n");
  }

  if (err < 1e-16) info = 0 ;


  printf("--------------------------\n");
  printf("test with transpose matrix\n");
  printf("--------------------------\n");

  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < n; j++)
    {
      W[i + j * m] = W[j + i * n];
    }
  }
  W[4 + 0 * m] = 2.0;
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      Wpinv[i + j * n] = 0.0;
      Wpinvtest[i + j * n] = 0.0;
    }
  }
  Wpinvtest[0 + 0 * n] = 0.2;
  Wpinvtest[0 + 4 * n] = 0.4;
  Wpinvtest[1 + 2 * n] = 1.0 / 3.0;
  Wpinvtest[3 + 1 * n] = 1.0 / 4.0;

  memcpy(Wpinv, W, n * m * sizeof(double));

  pinv(Wpinv, m, n, tol);
  printf("Winvtest\n");

  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      printf("%8.6e\t", Wpinvtest[i + j * n]) ;
    }
    printf("\n");
  }
  printf("Pseudo inverseWinv\n");
  err = 0.0;
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      printf("%8.6e\t", Wpinv[i + j * n]) ;
      err += (Wpinv[i + j * n] - Wpinvtest[i + j * n]) * (Wpinv[i + j * n] - Wpinvtest[i + j * n]);
    }
    printf("\n");
  }

  if (err < 1e-16) info = 0 ;


  FILE * file2 = fopen("dataWPseudoInverse.dat", "w");
  NM_write_in_file_scilab(WnumpInv, file2);
  fclose(file2);





  free(Wpinvtest);
  NM_free(Wnum);
  NM_free(WnumpInv);
  free(Wnum);
  free(WnumpInv);

  printf("-----------------------------------\n");
  printf("test with nearly 10*identity matrix\n");
  printf("-----------------------------------\n");
  n = 4;
  m = 4;
  W = (double*)malloc(n * m * sizeof(double));
  Wpinv = (double*)malloc(n * m * sizeof(double));
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      W[i + j * n] = 0.0;
    }
    W[i + i * n] = 10.0;
  }
  W[1 + 1 * n] = 1e-18;
  tol = 1e-16;
  memcpy(Wpinv, W, n * m * sizeof(double));
  pinv(Wpinv, n, m, tol);


  free(W);
  free(Wpinv);


  return info;

}
