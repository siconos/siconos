#include <math.h>                           // for sqrt
#include <stdio.h>                          // for printf, fopen, NULL, FILE
#include <stdlib.h>                         // for free, malloc, calloc
#include "FrictionContactProblem.h"         // for FrictionContactProblem
#include "NonSmoothDrivers.h"               // for variationalInequality_driver
#include "NumericsFwd.h"                    // for FrictionContactProblem
#include "NumericsMatrix.h"                 // for NM_gemv
#include "SolverOptions.h"                  // for solver_options_delete
#include "VI_cst.h"                         // for SICONOS_VI_EG, SICONOS_VI...
#include "VariationalInequality.h"          // for VariationalInequality
#include "VariationalInequality_Solvers.h"  // for variationalInequality_set...
#include "projectionOnCone.h"               // for projectionOnCone
#include "SiconosBlas.h"   // for cblas_dcopy

#pragma GCC diagnostic ignored "-Wmissing-prototypes"

typedef struct {
  VariationalInequality * vi;
  FrictionContactProblem * fc3d;
} Problems;


static void Ftest_0(void * self, int n_unused, double *x, double *F)
{
  VariationalInequality * vi = (VariationalInequality *) self;
  Problems* pb = (Problems *)vi->env;             
  FrictionContactProblem * fc3d = pb->fc3d;
  //frictionContact_display(fc3d);
 
  int nc = fc3d->numberOfContacts;
  int nLocal =  fc3d->dimension;
  int n = nc * nLocal; 

  cblas_dcopy(n , fc3d->q , 1 , F, 1);
  NM_gemv(1.0, fc3d->M, x, 1.0, F);
  int contact =0;
 
  for (contact = 0 ; contact < nc ; ++contact)
  {
    int pos = contact * nLocal;
    double  normUT = sqrt(F[pos + 1] * F[pos + 1] + F[pos + 2] * F[pos + 2]);
    F[pos] +=  (fc3d->mu[contact] * normUT);
  }
}


static void PXtest_0(void *viIn, double *x, double *PX)
{
  VariationalInequality * vi = (VariationalInequality *) viIn;
  Problems* pb = (Problems *)vi->env;             
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
  VariationalInequality vi;
  variationalInequality_clear(&vi);
  //vi.self = &vi;
  vi.F = &Ftest_0;
  vi.ProjectionOnX = &PXtest_0;
  vi.normVI = 0.0;
  vi.istheNormVIset = 0;
  vi.set = NULL;
  vi.nabla_F = NULL;
  SolverOptions * options = (SolverOptions *) malloc(sizeof(SolverOptions));
  int info = variationalInequality_setDefaultSolverOptions(options, SICONOS_VI_EG);
  options->dparam[0]=1e-8;

  FILE * finput  =  fopen("./data/FC3D_Example1_SBM.dat", "r");
  FrictionContactProblem* problem = (FrictionContactProblem *)malloc(sizeof(FrictionContactProblem));
  info = frictionContact_newFromFile(problem, finput);
  // frictionContact_display(problem);

  Problems *pb= (Problems *)malloc(sizeof(Problems));
  vi.env = pb;


  pb->vi = &vi;
  pb->fc3d = problem;
  frictionContact_display(pb->fc3d);

  int n = problem->numberOfContacts * problem->dimension;
  vi.size=n;

  double *x = (double*)calloc(n, sizeof(double));
  double *w = (double*)calloc(n, sizeof(double));

  PXtest_0(&vi, x,w);

 
  info = variationalInequality_driver(&vi, 
                                      x, 
                                      w, 
                                      options);
  int i =0;
  for (i =0; i< n ; i++)
  {
    printf("x[%i]=%f\t",i,x[i]);    printf("w[%i]=F[%i]=%f\n",i,i,w[i]);
  }
 
  solver_options_delete(options);
  free(options);
  free(problem);
  free(x);
  free(w);

  return info;
}


static void Ftest_1(void * self, int n_unused, double *x, double *F)
{
  VariationalInequality * vi = (VariationalInequality *) self;
  Problems* pb = (Problems *)vi->env;             
  FrictionContactProblem * fc3d = pb->fc3d;
  //frictionContact_display(fc3d);
 
  int nc = fc3d->numberOfContacts;
  int nLocal =  fc3d->dimension;
  int n = nc * nLocal; 

  cblas_dcopy(n , fc3d->q , 1 , F, 1);
  NM_gemv(1.0, fc3d->M, x, 1.0, F);
  int contact =0;
 
  for (contact = 0 ; contact < nc ; ++contact)
  {
    int pos = contact * nLocal;
    double  normUT = sqrt(F[pos + 1] * F[pos + 1] + F[pos + 2] * F[pos + 2]);
    F[pos] +=  (fc3d->mu[contact] * normUT);
  }
}


static void PXtest_1(void *viIn, double *x, double *PX)
{
  VariationalInequality * vi = (VariationalInequality *) viIn;
  Problems* pb = (Problems *)vi->env;             
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
  VariationalInequality vi;
  variationalInequality_clear(&vi);
  //vi.env = &vi;
  vi.F = &Ftest_1;
  vi.ProjectionOnX = &PXtest_1;
  vi.normVI = 0.0;
  vi.istheNormVIset = 0;
  vi.set = NULL;
  vi.nabla_F = NULL;
  SolverOptions * options = (SolverOptions *) malloc(sizeof(SolverOptions));
  int info = variationalInequality_setDefaultSolverOptions(options, SICONOS_VI_FPP);
  options->dparam[0]=1e-8;
  
  FILE * finput  =  fopen("./data/FC3D_Example1_SBM.dat", "r");
  FrictionContactProblem* problem = (FrictionContactProblem *)malloc(sizeof(FrictionContactProblem));
  info = frictionContact_newFromFile(problem, finput);
  // frictionContact_display(problem);

  Problems *pb= (Problems *)malloc(sizeof(Problems));
  vi.env = pb;


  pb->vi = &vi;
  pb->fc3d = problem;
  frictionContact_display(pb->fc3d);

  int n = problem->numberOfContacts * problem->dimension;
  vi.size=n;

  double *x = (double*)calloc(n, sizeof(double));
  double *w = (double*)calloc(n, sizeof(double));

  PXtest_1(&vi, x,w);

 
  info = variationalInequality_driver(&vi, 
                                      x, 
                                      w, 
                                      options);
  int i =0;
  for (i =0; i< n ; i++)
  {
    printf("x[%i]=%f\t",i,x[i]);    printf("w[%i]=F[%i]=%f\n",i,i,w[i]);
  }
 
  solver_options_delete(options);
  free(options);
  free(problem);
  free(x);
  free(w);

  return info;
}


static void Ftest_2(void * self, int n_unused, double *x, double *F)
{
  VariationalInequality * vi = (VariationalInequality *) self;
  Problems* pb = (Problems *)vi->env;             
  FrictionContactProblem * fc3d = pb->fc3d;
  //frictionContact_display(fc3d);
 
  int nc = fc3d->numberOfContacts;
  int nLocal =  fc3d->dimension;
  int n = nc * nLocal; 

  cblas_dcopy(n , fc3d->q , 1 , F, 1);
  NM_gemv(1.0, fc3d->M, x, 1.0, F);
  int contact =0;
 
  for (contact = 0 ; contact < nc ; ++contact)
  {
    int pos = contact * nLocal;
    double  normUT = sqrt(F[pos + 1] * F[pos + 1] + F[pos + 2] * F[pos + 2]);
    F[pos] +=  (fc3d->mu[contact] * normUT);
  }
}


static void PXtest_2(void *viIn, double *x, double *PX)
{
  VariationalInequality * vi = (VariationalInequality *) viIn;
  Problems* pb = (Problems *)vi->env;             
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
  VariationalInequality vi;
  variationalInequality_clear(&vi);
  //vi.env = &vi;
  vi.F = &Ftest_2;
  vi.ProjectionOnX = &PXtest_2;
  vi.normVI = 0.0;
  vi.istheNormVIset = 0;
  vi.set = NULL;
  vi.nabla_F = NULL;
  SolverOptions * options = (SolverOptions *) malloc(sizeof(SolverOptions));
  int info = variationalInequality_setDefaultSolverOptions(options, SICONOS_VI_HP);
  options->dparam[0]=1e-8;
  //char filename[50] = "./data/Confeti-ex13-Fc3D-SBM.dat";
  char filename[50] = "./data/FC3D_Example1_SBM.dat";
  
  FILE * finput  =  fopen(filename, "r");
  FrictionContactProblem* problem = (FrictionContactProblem *)malloc(sizeof(FrictionContactProblem));
  info = frictionContact_newFromFile(problem, finput);
  // frictionContact_display(problem);

  Problems *pb= (Problems *)malloc(sizeof(Problems));
  vi.env = pb;


  pb->vi = &vi;
  pb->fc3d = problem;
  frictionContact_display(pb->fc3d);

  int n = problem->numberOfContacts * problem->dimension;
  vi.size=n;

  double *x = (double*)calloc(n, sizeof(double));
  double *w = (double*)calloc(n, sizeof(double));

  PXtest_2(&vi, x,w);

 
  info = variationalInequality_driver(&vi, 
                                      x, 
                                      w, 
                                      options);
  int i =0;
  for (i =0; i< n ; i++)
  {
    printf("x[%i]=%f\t",i,x[i]);    printf("w[%i]=F[%i]=%f\n",i,i,w[i]);
  }
 
  solver_options_delete(options);
  free(options);
  free(problem);
  free(x);
  free(w);

  return info;
}



int main(void)
{

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
  info += test_1();
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
  int info_fail = test_2();
  
  if (!info_fail)
  {
    printf("end test #%i successful\n",i);
  }
  else
  {
    printf("end test #%i  not  successful (as predicted)\n",i);
  }
  return info;
}
