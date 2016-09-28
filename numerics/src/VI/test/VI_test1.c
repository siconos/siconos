#include <stdio.h>
#include <stdlib.h>
#include "NonSmoothDrivers.h"
#include "stdlib.h"
#include "VariationalInequality.h"
#include "SolverOptions.h"
#include "VI_cst.h"
#include "VariationalInequality_Solvers.h"

#pragma GCC diagnostic ignored "-Wmissing-prototypes"

void Ftest(void * viIn, int n, double *x, double *F)
{
  int i;
  VariationalInequality * vi = (VariationalInequality* ) viIn;
  for (i =0; i< vi->size ; i++)
  {
    F[i] = x[i]-i+4;
  }
}
void PXtest(void *viIn, double *x, double *PX)
{
  VariationalInequality * vi = (VariationalInequality* ) viIn;
  int i;
  for (i =0; i< vi->size ; i++)
  {
    PX[i] = x[i];
    if (PX[i] < 1.0) PX[i]=1.0;
  }
}




int main(void)
{

  VariationalInequality vi;
  variationalInequality_clear(&vi);

  vi.size=10;
  //vi.Callback = (CallbackVI *)malloc(sizeof(CallbackVI));

  vi.env = &vi;

  vi.F = &Ftest;
  vi.ProjectionOnX = &PXtest ;

  /* Call the callback */
  double x[10], F[10], PX[10];
  int i, n=10;
  for (i =0; i< n ; i++)
  {
    x[i] = i-5;
  }
  vi.F(&vi,n,x,F);
  vi.ProjectionOnX(&vi,x,PX);
  for (i =0; i< n ; i++)
  {
    printf("x[%i]=%f\t",i,x[i]);    printf("F[%i]=%f\t",i,F[i]);    printf("PX[%i]=%f\n",i,PX[i]);
  }
  SolverOptions * options = (SolverOptions *) malloc(sizeof(SolverOptions));

  int info = variationalInequality_setDefaultSolverOptions(options, SICONOS_VI_EG);
  options->dparam[0]=1e-10;
  

  info = variationalInequality_driver(&vi,
                                      x,
                                      F,
                                      options);

  for (i =0; i< n ; i++)
  {
    printf("x[%i]=%f\t",i,x[i]);    printf("w[%i]=F[%i]=%f\n",i,i,F[i]);
  }

  solver_options_delete(options);
  free(options);

  return info;
}
