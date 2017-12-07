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
#include "SparseMatrix.h"
#include "numerics_verbose.h"
#pragma GCC diagnostic ignored "-Wmissing-prototypes"


void PXtest(void *cqpIn, double *x, double *PX)
{
  ConvexQP * cqp = (ConvexQP* ) cqpIn;
  int i;
  for (i =0; i< cqp->size ; i++)
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
  //cqp.Callback = (CallbackCQP *)malloc(sizeof(CallbackCQP));

  cqp.env = &cqp;

  NumericsMatrix * M  = NM_create(NM_SPARSE,cqp.size, cqp.size);
  NM_triplet_alloc(M,0);
  M->matrix2->origin= NS_TRIPLET;

  for (int k =0; k< cqp.size; k++)
  {
    NM_zentry(M, k, k, 1);
  }
  NM_display(M);


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
  double x[10], u[10], xsi[10], PX[10];
  int i, n=cqp.size;
  for (i =0; i< n ; i++)
  {
    x[i] = i-5;
    u[i] = 0.0;
    xsi[i] =0.0;
  }


  cqp.ProjectionOnC(&cqp,x,PX);

  for (i =0; i< n ; i++)
  {
    printf("x[%i]=%f\t",i,x[i]);     printf("PX[%i]=%f\n",i,PX[i]);
  }
  for (i =0; i< n ; i++)
  {
    printf("q[%i]=%f\t",i,q[i]);
  }
  SolverOptions * options = (SolverOptions *) malloc(sizeof(SolverOptions));

  verbose=1;
  int info = convexQP_ADMM_setDefaultSolverOptions(options);

  options->dparam[SICONOS_DPARAM_TOL]=1e-14;
  //options->iparam[0]=30;
  options->dparam[SICONOS_CONVEXQP_ADMM_RHO]=1.0;
  printf("test step 1\n");
  convexQP_ADMM(&cqp, x, u, xsi, &info, options);
  //convexQP_ProjectedGradient(&cqp, x, w, &info, options);


  for (i =0; i< n ; i++)
  {
    printf("x[%i]=%f\t",i,x[i]);    printf("u[%i]=%f\t",i,u[i]); printf("xsi[%i]=%f\n",i,xsi[i]);
  }

  solver_options_delete(options);
  free(options);

  return info;
}
