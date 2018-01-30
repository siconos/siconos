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

//#define DEBUG_NOCOLOR
#define DEBUG_MESSAGES
#define DEBUG_STDOUT
#include "debug.h"


void PXtest(void *cqpIn, double *x, double *PX)
{
  ConvexQP * cqp = (ConvexQP* ) cqpIn;
  int i;
  for (i =0; i< cqp->m ; i++)
  {
    PX[i] = x[i];
    if (PX[i] < 4.0) PX[i]=4.0;
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
  M->matrix2->origin= NSM_TRIPLET;

  for (int k =0; k< cqp.size; k++)
  {
    NM_zentry(M, k, k, 1);
  }
  DEBUG_EXPR(NM_display(M));


  double * q = (double *) malloc(cqp.size*sizeof(double));
  for (int k =0; k< cqp.size; k++)
  {
    q[k]=-k-1;
  }


  cqp.m=5;
  NumericsMatrix * A  = NM_create(NM_SPARSE,cqp.m, cqp.size);
  NM_triplet_alloc(A,0);
  A->matrix2->origin= NSM_TRIPLET;

  for (int k =0; k< cqp.m; k++)
  {
    NM_zentry(A, k, k, 1);
  }
  DEBUG_EXPR(NM_display(A));


  double * b = (double *) malloc(cqp.size*sizeof(double));
  for (int k =0; k< cqp.m; k++)
  {
    b[k]=1.0;
  }

  printf("test step 1\n");
  cqp.M = M;
  cqp.ProjectionOnC = &PXtest ;
  cqp.q = q;
  cqp.A = A;
  cqp.b = b;
  convexQP_display(&cqp);


  /* Call the callback */
  double z[10], u[10], xi[10], w[10], PX[10];
  int i, n=cqp.size;
  for (i =0; i< n ; i++)
  {
    z[i] = i-5;
    u[i] = 0.0;
    xi[i] =0.0;
  }


  cqp.ProjectionOnC(&cqp,z,PX);

  for (i =0; i< n ; i++)
  {
    printf("z[%i]=%f\t",i,z[i]);     printf("PX[%i]=%f\n",i,PX[i]);
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
  options->iparam[SICONOS_CONVEXQP_ADMM_IPARAM_ACCELERATION]=SICONOS_CONVEXQP_ADMM_NO_ACCELERATION;
  printf("test step 1\n");
  convexQP_ADMM(&cqp, z, w, xi, u, &info, options);
  //convexQP_ProjectedGradient(&cqp, x, w, &info, options);


  for (i =0; i< n ; i++)
  {
    printf("z[%i]=%f\t",i,z[i]);printf("w[%i]=%f\n",i,w[i]);
  }
  for (i =0; i< cqp.m ; i++)
  {
    printf("u[%i]=%f\t",i,u[i]); printf("xi[%i]=%f\n",i,xi[i]);
  }

  solver_options_delete(options);
  free(options);

  return info;
}
