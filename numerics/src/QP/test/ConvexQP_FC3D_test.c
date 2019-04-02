#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "NonSmoothDrivers.h"
#include "projectionOnCone.h"
#include "ConvexQP.h"
#include "FrictionContactProblem.h"
#include "SolverOptions.h"
#include "ConvexQP_cst.h"
#include "ConvexQP_Solvers.h"
#include "SiconosBlas.h"
#include "NumericsMatrix.h"
#include "numerics_verbose.h"
#pragma GCC diagnostic ignored "-Wmissing-prototypes"

typedef struct {
  ConvexQP * cqp;
  FrictionContactProblem * fc3d;
} Problems;


static void PXtest_0(void *cqpIn, double *x, double *PX)
{
  ConvexQP * cqp = (ConvexQP *) cqpIn;
  Problems* pb = (Problems *)cqp->env;
  FrictionContactProblem * fc3d = pb->fc3d;
  //frictionContact_display(fc3d);

  int contact =0;
  int nc = fc3d->numberOfContacts;
  int nLocal =  fc3d->dimension;
  int n = nc * nLocal;

  cblas_dcopy(n , x , 1 , PX, 1);

  for (contact = 0 ; contact < nc ; ++contact)
  {
    int pos = contact * nLocal;
    projectionOnCone(&PX[pos], fc3d->mu[contact]);
  }
}



static int test_0(void)
{
  ConvexQP cqp;
  convexQP_clear(&cqp);
  //cqp.env = &cqp;
  cqp.ProjectionOnC = &PXtest_0;

  SolverOptions * options = (SolverOptions *) malloc(sizeof(SolverOptions));

  verbose=1;
  int info = convexQP_ProjectedGradient_setDefaultSolverOptions(options);
  options->dparam[0]=1e-13;

  char filename[50] = "./data/FC3D_Example1_SBM.dat";
  FILE * finput  =  fopen(filename, "r");
  FrictionContactProblem* problem = (FrictionContactProblem *)malloc(sizeof(FrictionContactProblem));
  info = frictionContact_newFromFile(problem, finput);
  // frictionContact_display(problem);


  Problems *pb= (Problems *)malloc(sizeof(Problems));
  cqp.env = pb;


  pb->cqp = &cqp;
  pb->fc3d = problem;
  frictionContact_display(pb->fc3d);

  int n = problem->numberOfContacts * problem->dimension;
  cqp.size=n;

  cqp.M=problem->M;
  cqp.q=problem->q;
  cqp.m=problem->numberOfContacts * problem->dimension;


  cqp.A=NULL;


  double *z = (double*)calloc(n, sizeof(double));
  double *w = (double*)calloc(n, sizeof(double));


  PXtest_0(&cqp, z,w);

  convexQP_ProjectedGradient(&cqp, z, w, &info, options);
  int i =0;
  for (i =0; i< n ; i++)
  {
    printf("x[%i]=%f\t",i,z[i]);    printf("w[%i]=F[%i]=%f\n",i,i,w[i]);
  }

  solver_options_delete(options);
  free(options);
  free(problem);
  free(z);
  free(w);

  return info;
}




static void PXtest_1(void *cqpIn, double *x, double *PX)
{
  ConvexQP * cqp = (ConvexQP *) cqpIn;
  Problems* pb = (Problems *)cqp->env;
  FrictionContactProblem * fc3d = pb->fc3d;
  //frictionContact_display(fc3d);

  int contact =0;
  int nc = fc3d->numberOfContacts;
  int nLocal =  fc3d->dimension;
  int n = nc * nLocal;

  cblas_dcopy(n , x , 1 , PX, 1);

  for (contact = 0 ; contact < nc ; ++contact)
  {
    int pos = contact * nLocal;
    projectionOnCone(&PX[pos], fc3d->mu[contact]);
  }
}



static int test_1(void)
{
  ConvexQP cqp;
  convexQP_clear(&cqp);
  //cqp.env = &cqp;
  cqp.ProjectionOnC = &PXtest_1;

  SolverOptions * options = (SolverOptions *) malloc(sizeof(SolverOptions));

  verbose=1;
  int info = convexQP_ADMM_setDefaultSolverOptions(options);
  options->iparam[SICONOS_CONVEXQP_ADMM_IPARAM_ACCELERATION]=SICONOS_CONVEXQP_ADMM_NO_ACCELERATION;
  options->dparam[0]=1e-13;
  options->dparam[3]=0.8;

  char filename[50] = "./data/FC3D_Example1_SBM.dat";
  FILE * finput  =  fopen(filename, "r");
  FrictionContactProblem* problem = (FrictionContactProblem *)malloc(sizeof(FrictionContactProblem));
  info = frictionContact_newFromFile(problem, finput);
  // frictionContact_display(problem);


  Problems *pb= (Problems *)malloc(sizeof(Problems));
  cqp.env = pb;


  pb->cqp = &cqp;
  pb->fc3d = problem;
  frictionContact_display(pb->fc3d);

  int n = problem->numberOfContacts * problem->dimension;
  cqp.size=n;

  cqp.M=problem->M;
  cqp.q=problem->q;
  cqp.m=problem->numberOfContacts * problem->dimension;


  cqp.A=NULL;


  double *z = (double*)calloc(n, sizeof(double));
  double *w = (double*)calloc(n, sizeof(double));
  double *u = (double*)calloc(n, sizeof(double));
  double *xi = (double*)calloc(n, sizeof(double));

  PXtest_1(&cqp, z,w);

  convexQP_ADMM(&cqp, z, w, xi, u, &info, options);

  int i =0;
  printf("--------- \n\n");
  for (i =0; i< n ; i++)
  {
    printf("x[%i]=%f\t",i,z[i]);    printf("w[%i]=A^Txi[%i]=%f\n",i,i,w[i]);
  }
  printf("--------- \n\n");
  for (i= 0; i< n ; i++)
  {
    printf("xi[%i]=%f\t",i,xi[i]);    printf("u[%i]=%f\n",i,u[i]);
  }
  if (!info)
  {
    printf("test successful, residual = %g\n", options->dparam[1]);
  }
  else
  {
    printf("test unsuccessful, residual = %g\n", options->dparam[1]);
  }
  solver_options_delete(options);
  free(options);
  free(problem);
  free(z);
  free(w);

  return info;
}



static void PXtest_2(void *cqpIn, double *x, double *PX)
{
  ConvexQP * cqp = (ConvexQP *) cqpIn;
  Problems* pb = (Problems *)cqp->env;
  FrictionContactProblem * fc3d = pb->fc3d;
  //frictionContact_display(fc3d);

  int contact =0;
  int nc = fc3d->numberOfContacts;
  int nLocal =  fc3d->dimension;
  int n = nc * nLocal;

  cblas_dcopy(n , x , 1 , PX, 1);

  for (contact = 0 ; contact < nc ; ++contact)
  {
    int pos = contact * nLocal;
    projectionOnCone(&PX[pos], fc3d->mu[contact]);
  }
}



static int test_2(void)
{
  ConvexQP cqp;
  convexQP_clear(&cqp);
  //cqp.env = &cqp;
  cqp.ProjectionOnC = &PXtest_2;

  SolverOptions * options = (SolverOptions *) malloc(sizeof(SolverOptions));

  verbose=1;
  int info = convexQP_ADMM_setDefaultSolverOptions(options);
  options->dparam[0]=1e-13;
  options->dparam[3]=0.8;
  options->iparam[SICONOS_CONVEXQP_ADMM_IPARAM_ACCELERATION] = SICONOS_CONVEXQP_ADMM_ACCELERATION;

  char filename[50] = "./data/FC3D_Example1_SBM.dat";
  FILE * finput  =  fopen(filename, "r");
  FrictionContactProblem* problem = (FrictionContactProblem *)malloc(sizeof(FrictionContactProblem));
  info = frictionContact_newFromFile(problem, finput);
  // frictionContact_display(problem);


  Problems *pb= (Problems *)malloc(sizeof(Problems));
  cqp.env = pb;


  pb->cqp = &cqp;
  pb->fc3d = problem;
  frictionContact_display(pb->fc3d);

  int n = problem->numberOfContacts * problem->dimension;
  cqp.size=n;

  cqp.M=problem->M;
  cqp.q=problem->q;
  cqp.m=problem->numberOfContacts * problem->dimension;


  cqp.A=NULL;


  double *z = (double*)calloc(n, sizeof(double));
  double *w = (double*)calloc(n, sizeof(double));
  double *u = (double*)calloc(n, sizeof(double));
  double *xi = (double*)calloc(n, sizeof(double));

  PXtest_2(&cqp, z,w);

  convexQP_ADMM(&cqp, z, w, xi, u, &info, options);

  int i =0;
  printf("--------- \n\n");
  for (i =0; i< n ; i++)
  {
    printf("x[%i]=%f\t",i,z[i]);    printf("w[%i]=A^Txi[%i]=%f\n",i,i,w[i]);
  }
  printf("--------- \n\n");
  for (i= 0; i< n ; i++)
  {
    printf("xi[%i]=%f\t",i,xi[i]);    printf("u[%i]=%f\n",i,u[i]);
  }
  if (!info)
  {
    printf("test successful, residual = %g\n", options->dparam[1]);
  }
  else
  {
    printf("test unsuccessful, residual = %g\n", options->dparam[1]);
  }
  solver_options_delete(options);
  free(options);
  free(problem);
  free(z);
  free(w);

  return info;
}



static void PXtest_3(void *cqpIn, double *x, double *PX)
{
  ConvexQP * cqp = (ConvexQP *) cqpIn;
  Problems* pb = (Problems *)cqp->env;
  FrictionContactProblem * fc3d = pb->fc3d;
  //frictionContact_display(fc3d);

  int contact =0;
  int nc = fc3d->numberOfContacts;
  int nLocal =  fc3d->dimension;
  int n = nc * nLocal;

  cblas_dcopy(n , x , 1 , PX, 1);

  for (contact = 0 ; contact < nc ; ++contact)
  {
    int pos = contact * nLocal;
    projectionOnCone(&PX[pos], fc3d->mu[contact]);
  }
}



static int test_3(void)
{
  ConvexQP cqp;
  convexQP_clear(&cqp);
  //cqp.env = &cqp;
  cqp.ProjectionOnC = &PXtest_3;

  SolverOptions * options = (SolverOptions *) malloc(sizeof(SolverOptions));

  verbose=1;
  int info = convexQP_ADMM_setDefaultSolverOptions(options);
  options->dparam[0]=1e-13;
  options->dparam[3]=0.8;
  options->iparam[SICONOS_CONVEXQP_ADMM_IPARAM_ACCELERATION] = SICONOS_CONVEXQP_ADMM_ACCELERATION_AND_RESTART;

  char filename[50] = "./data/FC3D_Example1_SBM.dat";
  FILE * finput  =  fopen(filename, "r");
  FrictionContactProblem* problem = (FrictionContactProblem *)malloc(sizeof(FrictionContactProblem));
  info = frictionContact_newFromFile(problem, finput);
  // frictionContact_display(problem);


  Problems *pb= (Problems *)malloc(sizeof(Problems));
  cqp.env = pb;


  pb->cqp = &cqp;
  pb->fc3d = problem;
  frictionContact_display(pb->fc3d);

  int n = problem->numberOfContacts * problem->dimension;
  cqp.size=n;

  cqp.M=problem->M;
  cqp.q=problem->q;
  cqp.m=problem->numberOfContacts * problem->dimension;


  cqp.A=NULL;


  double *z = (double*)calloc(n, sizeof(double));
  double *w = (double*)calloc(n, sizeof(double));
  double *u = (double*)calloc(n, sizeof(double));
  double *xi = (double*)calloc(n, sizeof(double));

  PXtest_3(&cqp, z,w);

  convexQP_ADMM(&cqp, z, w, xi, u, &info, options);

  int i =0;
  printf("--------- \n\n");
  for (i =0; i< n ; i++)
  {
    printf("x[%i]=%f\t",i,z[i]);    printf("w[%i]=A^Txi[%i]=%f\n",i,i,w[i]);
  }
  printf("--------- \n\n");
  for (i= 0; i< n ; i++)
  {
    printf("xi[%i]=%f\t",i,xi[i]);    printf("u[%i]=%f\n",i,u[i]);
  }
  if (!info)
  {
    printf("test successful, residual = %g\n", options->dparam[1]);
  }
  else
  {
    printf("test unsuccessful, residual = %g\n", options->dparam[1]);
  }
  solver_options_delete(options);
  free(options);
  free(problem);
  free(z);
  free(w);

  return info;
}


int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  int i=0;
  printf("start test #%i ConvexQP_PG_FC3D \n",i);
  int info = test_0();
  if (!info)
  {
    printf("end test #%i sucessfull\n",i);
  }
  else
  {
    printf("end test #%i  not  sucessfull\n",i);
  }

  i++;
  printf("start test #%i ConvexQP_ADMM_FC3D\n",i);
  info += test_1();
  if (!info)
  {
    printf("end test #%i sucessfull\n",i);
  }
  else
  {
    printf("end test #%i  not  sucessfull\n",i);
  }
  i++;
  printf("start test #%i ConvexQP_ADMM_ACCELERATION_FC3D\n",i);
  info += test_2();
  if (!info)
  {
    printf("end test #%i sucessfull\n",i);
  }
  else
  {
    printf("end test #%i  not  sucessfull\n",i);
  }
  i++;
  printf("start test #%i ConvexQP_ADMM_ACCELERATION_AND_RESTART_FC3D\n",i);
  info += test_3();
  if (!info)
  {
    printf("end test #%i sucessfull\n",i);
  }
  else
  {
    printf("end test #%i  not  sucessfull\n",i);
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return info;
}
