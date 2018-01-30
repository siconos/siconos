#include <stdio.h>
#include <stdlib.h>
#include "NonSmoothDrivers.h"
#include "stdlib.h"
#include "ConvexQP.h"
#include "SolverOptions.h"
#include "ConvexQP_cst.h"
#include "ConvexQP_Solvers.h"
#include "NumericsMatrix.h"
#include "NumericsSparseMatrix.h"
#include "CSparseMatrix.h"
#include "numerics_verbose.h"
#pragma GCC diagnostic ignored "-Wmissing-prototypes"


void PXtest(void *cqpIn, double *x, double *PX)
{
  ConvexQP * cqp = (ConvexQP* ) cqpIn;
  int i;
  for (i =0; i< cqp->m ; i++)
  {
    PX[i] = x[i];
    if (PX[i] < 3.0) PX[i]=3.0;
  }
}

int main(void)
{

  ConvexQP cqp;

  convexQP_clear(&cqp);
  cqp.size=10;
  cqp.m=10;
  //cqp.Callback = (CallbackCQP *)malloc(sizeof(CallbackCQP));

  cqp.env = &cqp;

   NumericsMatrix * M  = NM_create(NM_SPARSE,cqp.size, cqp.size);
  NM_triplet_alloc(M,0);
  M->matrix2->origin= NSM_TRIPLET;

  for (int k =0; k< cqp.size; k++)
  {
    NM_zentry(M, k, k, 1);
  }
  /* NM_display(M); */


  double * q = (double *) malloc(cqp.size*sizeof(double));
  for (int k =0; k< cqp.size; k++)
  {
    q[k]=-k-1;
  }


  printf("test step 1\n");
  cqp.M = M;
  cqp.ProjectionOnC = &PXtest ;
  cqp.q = q;
  convexQP_display(&cqp);


  /* Call the callback */
  double x[10], w[10], PX[10];
  int i, n=10;
  for (i =0; i< n ; i++)
  {
    x[i] = i-5;
  }


  cqp.ProjectionOnC(&cqp,x,PX);
  for (i =0; i< n ; i++)
  {
    printf("x[%i]=%f\t",i,x[i]);     printf("PX[%i]=%f\n",i,PX[i]);
  }
  SolverOptions * options = (SolverOptions *) malloc(sizeof(SolverOptions));

  verbose=1;
  int info = convexQP_ProjectedGradient_setDefaultSolverOptions(options);

  options->dparam[0]=1e-12;
  options->dparam[3]=1.0;

  convexQP_ProjectedGradient(&cqp, x, w, &info, options);


  for (i =0; i< n ; i++)
  {
    printf("x[%i]=%f\t",i,x[i]);    printf("w[%i]=w[%i]=%f\n",i,i,w[i]);
  }

  solver_options_delete(options);
  free(options);

  return info;
}
