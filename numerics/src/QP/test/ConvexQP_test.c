#include "ConvexQP.h"
#include "stdlib.h"
#include "NumericsMatrix.h"
#include "NumericsSparseMatrix.h"
#include "CSparseMatrix.h"
#include "ConvexQP_cst.h"
#include "ConvexQP_Solvers.h"

#include "SolverOptions.h"
#include "numerics_verbose.h"


//#define DEBUG_NOCOLOR
//#define DEBUG_MESSAGES
//#define DEBUG_STDOUT
#include "debug.h"



static void PXtest_0(void *cqpIn, double *x, double *PX)
{
  ConvexQP * cqp = (ConvexQP* ) cqpIn;
  printf("Size of cqp :%i\n", cqp->size);
  int i;
  for (i =0; i< cqp->size ; i++)
  {
    PX[i] = x[i];
    if (PX[i] <0) PX[i]=0.0;
  }
}




static int test_0(void)
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
    q[k]=k;
  }


  printf("test step 1\n");
  cqp.M = M;
  cqp.ProjectionOnC = &PXtest_0 ;
  cqp.q = q;
  cqp.A = NULL;
  cqp.b = NULL;

  /* convexQP_display(&cqp); */

  /* Call the callback */
  double x[10],  PX[10];
  int i, n=cqp.size;
  for (i =0; i< n ; i++)
  {
    x[i] = i-5;
  }
  cqp.ProjectionOnC(&cqp,x,PX);
  for (i =0; i< n ; i++)
  {
    printf("x[%i]=%f\t",i,x[i]);   printf("PX[%i]=%f\n",i,PX[i]);
  }

  return 0;
}


static void PXtest_1(void *cqpIn, double *x, double *PX)
{
  ConvexQP * cqp = (ConvexQP* ) cqpIn;
  int i;
  for (i =0; i< cqp->m ; i++)
  {
    PX[i] = x[i];
    if (PX[i] < 3.0) PX[i]=3.0;
  }
}

static int test_1(void)
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
  cqp.ProjectionOnC = &PXtest_1 ;
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

  /* verbose=1; */
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

static void PXtest_2(void *cqpIn, double *x, double *PX)
{
  ConvexQP * cqp = (ConvexQP* ) cqpIn;
  int i;
  for (i =0; i< cqp->m ; i++)
  {
    PX[i] = x[i];
    if (PX[i] < 3.0) PX[i]=3.0;
  }
}

static int test_2(void)
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
  /* NM_display(M); */


  double * q = (double *) malloc(cqp.size*sizeof(double));
  for (int k =0; k< cqp.size; k++)
  {
    q[k]=-k-1;
  }


  printf("test step 1\n");
  cqp.M = M;
  cqp.ProjectionOnC = &PXtest_2 ;
  cqp.q = q;
  /* convexQP_display(&cqp); */


  /* Call the callback */
  double z[10], w[10], u[10], xi[10], PX[10];
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

  /* verbose=1; */
  int info = convexQP_ADMM_setDefaultSolverOptions(options);

  options->dparam[SICONOS_DPARAM_TOL]=1e-14;
  //options->iparam[0]=30;
  options->dparam[SICONOS_CONVEXQP_ADMM_RHO]=1.0;
  printf("test step 1\n");
  convexQP_ADMM(&cqp, z, w, xi, u, &info, options);
  //convexQP_ProjectedGradient(&cqp, x, w, &info, options);


  for (i =0; i< n ; i++)
  {
    printf("z[%i]=%f\t",i,z[i]); printf("w[%i]=%f\t",i,w[i]);    printf("u[%i]=%f\t",i,u[i]); printf("xi[%i]=%f\n",i,xi[i]);
  }

  solver_options_delete(options);
  free(options);

  return info;
}

static void PXtest_3(void *cqpIn, double *x, double *PX)
{
  ConvexQP * cqp = (ConvexQP* ) cqpIn;
  int i;
  for (i =0; i< cqp->m ; i++)
  {
    PX[i] = x[i];
    if (PX[i] < 4.0) PX[i]=4.0;
  }
}

static int test_3(void)
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
  /* DEBUG_EXPR(NM_display(M)); */


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
  /* DEBUG_EXPR(NM_display(A)); */


  double * b = (double *) malloc(cqp.size*sizeof(double));
  for (int k =0; k< cqp.m; k++)
  {
    b[k]=1.0;
  }

  printf("test step 1\n");
  cqp.M = M;
  cqp.ProjectionOnC = &PXtest_3 ;
  cqp.q = q;
  cqp.A = A;
  cqp.b = b;
  /* convexQP_display(&cqp); */


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

  /* verbose=1; */
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


static void PXtest_4(void *cqpIn, double *x, double *PX)
{
  ConvexQP * cqp = (ConvexQP* ) cqpIn;
  int i;
  for (i =0; i< cqp->m ; i++)
  {
    PX[i] = x[i];
    if (PX[i] < 4.0) PX[i]=4.0;
  }
}

static int test_4(void)
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
  /* DEBUG_EXPR(NM_display(M)); */


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
  /* DEBUG_EXPR(NM_display(A)); */


  double * b = (double *) malloc(cqp.size*sizeof(double));
  for (int k =0; k< cqp.m; k++)
  {
    b[k]=1.0;
  }

  printf("test step 1\n");
  cqp.M = M;
  cqp.ProjectionOnC = &PXtest_4 ;
  cqp.q = q;
  cqp.A = A;
  cqp.b = b;
  /* convexQP_display(&cqp); */


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

  /* verbose=1; */
  int info = convexQP_ADMM_setDefaultSolverOptions(options);

  options->dparam[SICONOS_DPARAM_TOL]=1e-14;
  //options->iparam[0]=30;
  options->dparam[SICONOS_CONVEXQP_ADMM_RHO]=1.0;
  options->iparam[SICONOS_CONVEXQP_ADMM_IPARAM_ACCELERATION] = SICONOS_CONVEXQP_ADMM_ACCELERATION;
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

static void PXtest_5(void *cqpIn, double *x, double *PX)
{
  ConvexQP * cqp = (ConvexQP* ) cqpIn;
  int i;
  for (i =0; i< cqp->m ; i++)
  {
    PX[i] = x[i];
    if (PX[i] < 4.0) PX[i]=4.0;
  }
}

static int test_5(void)
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
  /* DEBUG_EXPR(NM_display(M)); */


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
  /* DEBUG_EXPR(NM_display(A)); */


  double * b = (double *) malloc(cqp.size*sizeof(double));
  for (int k =0; k< cqp.m; k++)
  {
    b[k]=1.0;
  }

  printf("test step 1\n");
  cqp.M = M;
  cqp.ProjectionOnC = &PXtest_5 ;
  cqp.q = q;
  cqp.A = A;
  cqp.b = b;
  /* convexQP_display(&cqp); */


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

  /* verbose=1; */
  int info = convexQP_ADMM_setDefaultSolverOptions(options);

  options->dparam[SICONOS_DPARAM_TOL]=1e-14;
  //options->iparam[0]=30;
  options->dparam[SICONOS_CONVEXQP_ADMM_RHO]=1.0;
  options->iparam[SICONOS_CONVEXQP_ADMM_IPARAM_ACCELERATION] = SICONOS_CONVEXQP_ADMM_ACCELERATION_AND_RESTART;
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

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  int i=0;
  printf("start test #%i\n",i);
  int info = test_0();
  if (!info)
  {
    printf("end test #%i successful\n",i);
  }
  else
  {
    printf("end test #%i  not  successful\n",i);
  }

  i++;
  printf("start test #%i\n",i);
  int info_test = test_1();
  info += info_test;
  if (!info_test)
  {
    printf("end test #%i successful\n",i);
  }
  else
  {
    printf("end test #%i  not  successful\n",i);
  }

  i++;
  printf("start test #%i\n",i);
  info_test = test_2();


#ifndef WITH_MUMPS
  info += info_test;
#endif
  if (!info_test)
  {
    printf("end test #%i successful\n",i);
  }
  else
  {
#ifndef WITH_MUMPS
    printf("end test #%i  not  successful\n",i);
#else
    printf("end test #%i  not  successful (as predicted with mumps)\n",i);
#endif
  }



  i++;
  printf("start test #%i ConvexQP_ADDM\n",i);
  info_test = test_3();
  info += info_test;
  if (!info_test)
  {
    printf("end test #%i successful\n",i);
  }
  else
  {
    printf("end test #%i  not  successful\n",i);
  }


  i++;
  printf("start test #%i ConvexQP_ADDM_ACCELERATION\n",i);
  info_test = test_4();
  info += info_test;
  if (!info)
  {
    printf("end test #%i successful\n",i);
  }
  else
  {
    printf("end test #%i  not  successful\n",i);
  }

  i++;
  printf("start test #%i ConvexQP_ADDM_ACCELERATION_AND_RESTART\n",i);
  info_test = test_5();
  info += info_test;
  if (!info)
  {
    printf("end test #%i successful\n",i);
  }
  else
  {
    printf("end test #%i  not  successful\n",i);
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return info;
}
