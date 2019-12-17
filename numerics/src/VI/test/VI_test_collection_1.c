/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2019 INRIA.
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
#include <stdio.h>                  // for printf
#include "NonSmoothDrivers.h"       // for variationalInequality_driver
#include "NumericsFwd.h"            // for VariationalInequality, SolverOptions
#include "SolverOptions.h"          // for solver_options_delete, solver_opti...
#include "VI_cst.h"                 // for SICONOS_VI_HP, SICONOS_VI_EG, SIC...
#include "VariationalInequality.h"  // for VariationalInequality, variationa...

static void Ftest_0(void * viIn, int n, double *x, double *F)
{
  int i;
  VariationalInequality * vi = (VariationalInequality*) viIn;
  printf("Size of vi :%i\n", vi->size);

  for(i =0; i< vi->size ; i++)
  {
    F[i] = x[i];
  }
}
static void PXtest_0(void *viIn, double *x, double *PX)
{
  VariationalInequality * vi = (VariationalInequality*) viIn;
  printf("Size of vi :%i\n", vi->size);
  int i;
  for(i =0; i< vi->size ; i++)
  {
    PX[i] = x[i];
    if(PX[i] <0) PX[i]=0.0;
  }
}




static int test_0(void)
{

  VariationalInequality vi;
  variationalInequality_clear(&vi);

  vi.size=10;
  //vi.Callback = (CallbackVI *)malloc(sizeof(CallbackVI));

  vi.env = &vi;

  vi.F = &Ftest_0;
  vi.ProjectionOnX = &PXtest_0 ;

  /* Call the callback */
  double x[10], F[10], PX[10];
  int i, n=10;
  for(i =0; i< n ; i++)
  {
    x[i] = i-5;
  }
  vi.F(&vi,n,x,F);
  vi.ProjectionOnX(&vi,x,PX);
  for(i =0; i< n ; i++)
  {
    printf("x[%i]=%f\t",i,x[i]);
    printf("F[%i]=%f\t",i,F[i]);
    printf("PX[%i]=%f\n",i,PX[i]);
  }

  return 0;

}

static void Ftest_1(void * viIn, int n, double *x, double *F)
{
  int i;
  VariationalInequality * vi = (VariationalInequality*) viIn;
  for(i =0; i< vi->size ; i++)
  {
    F[i] = x[i]-i+4;
  }
}
static void PXtest_1(void *viIn, double *x, double *PX)
{
  VariationalInequality * vi = (VariationalInequality*) viIn;
  int i;
  for(i =0; i< vi->size ; i++)
  {
    PX[i] = x[i];
    if(PX[i] < 1.0) PX[i]=1.0;
  }
}




static int test_1(void)
{

  VariationalInequality vi;
  variationalInequality_clear(&vi);

  vi.size=10;
  //vi.Callback = (CallbackVI *)malloc(sizeof(CallbackVI));

  vi.env = &vi;

  vi.F = &Ftest_1;
  vi.ProjectionOnX = &PXtest_1 ;

  /* Call the callback */
  double x[10], F[10], PX[10];
  int i, n=10;
  for(i =0; i< n ; i++)
  {
    x[i] = i-5;
  }
  vi.F(&vi,n,x,F);
  vi.ProjectionOnX(&vi,x,PX);
  for(i =0; i< n ; i++)
  {
    printf("x[%i]=%f\t",i,x[i]);
    printf("F[%i]=%f\t",i,F[i]);
    printf("PX[%i]=%f\n",i,PX[i]);
  }
  SolverOptions * options = solver_options_create(SICONOS_VI_EG);
  options->dparam[SICONOS_DPARAM_TOL]=1e-10;


  int info = variationalInequality_driver(&vi,
                                          x,
                                          F,
                                          options);

  for(i =0; i< n ; i++)
  {
    printf("x[%i]=%f\t",i,x[i]);
    printf("w[%i]=F[%i]=%f\n",i,i,F[i]);
  }

  solver_options_delete(options);

  return info;
}

static void Ftest_2(void * viIn, int n, double *x, double *F)
{
  int i;
  VariationalInequality * vi = (VariationalInequality*) viIn;
  for(i =0; i< vi->size ; i++)
  {
    F[i] = x[i]-i+4;
  }
}
static  void PXtest_2(void *viIn, double *x, double *PX)
{
  VariationalInequality * vi = (VariationalInequality*) viIn;
  int i;
  for(i =0; i< vi->size ; i++)
  {
    PX[i] = x[i];
    if(PX[i] < 1.0) PX[i]=1.0;
  }
}




static int test_2(void)
{

  VariationalInequality vi;
  variationalInequality_clear(&vi);

  vi.size=10;
  //vi.Callback = (CallbackVI *)malloc(sizeof(CallbackVI));

  vi.env = &vi;

  vi.F = &Ftest_2;
  vi.ProjectionOnX = &PXtest_2 ;

  /* Call the callback */
  double x[10], F[10], PX[10];
  int i, n=10;
  for(i =0; i< n ; i++)
  {
    x[i] = i-5;
  }
  vi.F(&vi,n,x,F);
  vi.ProjectionOnX(&vi,x,PX);
  for(i =0; i< n ; i++)
  {
    printf("x[%i]=%f\t",i,x[i]);
    printf("F[%i]=%f\t",i,F[i]);
    printf("PX[%i]=%f\n",i,PX[i]);
  }
  SolverOptions * options = solver_options_create(SICONOS_VI_FPP);
  options->dparam[SICONOS_DPARAM_TOL]=1e-10;

  int info = variationalInequality_driver(&vi,
                                          x,
                                          F,
                                          options);

  for(i =0; i< n ; i++)
  {
    printf("x[%i]=%f\t",i,x[i]);
    printf("w[%i]=F[%i]=%f\n",i,i,F[i]);
  }

  solver_options_delete(options);

  return info;
}


static void Ftest_3(void * viIn, int n, double *x, double *F)
{
  int i;
  VariationalInequality * vi = (VariationalInequality*) viIn;
  for(i =0; i< vi->size ; i++)
  {

    F[i] = x[i]-i+4;
  }
}
static void PXtest_3(void *viIn, double *x, double *PX)
{
  VariationalInequality * vi = (VariationalInequality*) viIn;
  int i;
  for(i =0; i< vi->size ; i++)
  {
    PX[i] = x[i];
    if(PX[i] < 1.0) PX[i]=1.0;
  }
}




static int test_3(void)
{

  VariationalInequality vi;
  variationalInequality_clear(&vi);

  vi.size=10;
  //vi.Callback = (CallbackVI *)malloc(sizeof(CallbackVI));

  vi.env = &vi;

  vi.F = &Ftest_3;
  vi.ProjectionOnX = &PXtest_3 ;

  /* Call the callback */
  double x[10], F[10], PX[10];
  int i, n=10;
  for(i =0; i< n ; i++)
  {
    x[i] = i-5;
  }
  vi.F(&vi,n,x,F);
  vi.ProjectionOnX(&vi,x,PX);
  for(i =0; i< n ; i++)
  {
    printf("x[%i]=%f\t",i,x[i]);
    printf("F[%i]=%f\t",i,F[i]);
    printf("PX[%i]=%f\n",i,PX[i]);
  }
  SolverOptions * options = solver_options_create(SICONOS_VI_HP);
  options->dparam[SICONOS_DPARAM_TOL]=1e-02;
  options->iparam[SICONOS_IPARAM_MAX_ITER]=100000;


  int info = variationalInequality_driver(&vi,
                                          x,
                                          F,
                                          options);

  for(i =0; i< n ; i++)
  {
    printf("x[%i]=%f\t",i,x[i]);
    printf("w[%i]=F[%i]=%f\n",i,i,F[i]);
  }

  solver_options_delete(options);

  return info;
}

static void Ftest_4(void * viIn, int n, double *x, double *F)
{
  int i;
  VariationalInequality * vi = (VariationalInequality*) viIn;
  for(i =0; i< vi->size ; i++)
  {
    F[i] = x[i]-4;
  }
}
static void PXtest_4(void *viIn, double *x, double *PX)
{
  VariationalInequality * vi = (VariationalInequality*) viIn;
  int i;
  for(i =0; i< vi->size ; i++)
  {
    PX[i] = x[i];
    if(PX[i] < 0.0) PX[i]=0.0;
  }
}




static int test_4(void)
{

  VariationalInequality vi;
  variationalInequality_clear(&vi);

  vi.size=1;
  //vi.Callback = (CallbackVI *)malloc(sizeof(CallbackVI));

  vi.env = &vi;

  vi.F = &Ftest_4;
  vi.ProjectionOnX = &PXtest_4 ;

  /* Call the callback */
  double x[1], F[1], PX[1];
  int i, n=1;
  for(i =0; i< n ; i++)
  {
    x[i] = i-5;
  }
  vi.F(&vi,n,x,F);
  vi.ProjectionOnX(&vi,x,PX);
  for(i =0; i< n ; i++)
  {
    printf("x[%i]=%f\t",i,x[i]);
    printf("F[%i]=%f\t",i,F[i]);
    printf("PX[%i]=%f\n",i,PX[i]);
  }
  SolverOptions * options = solver_options_create(SICONOS_VI_HP);

  options->dparam[SICONOS_DPARAM_TOL]=1e-10;
  options->iparam[SICONOS_IPARAM_MAX_ITER]=50000000;


  int info = variationalInequality_driver(&vi,
                                          x,
                                          F,
                                          options);

  for(i =0; i< n ; i++)
  {
    printf("x[%i]=%f\t",i,x[i]);
    printf("w[%i]=F[%i]=%f\n",i,i,F[i]);
  }

  solver_options_delete(options);
  options = NULL;

  return info;
}

int main(void)
{

  int i=0;
  printf("start test #%i\n",i);
  int info = test_0();
  if(!info)
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
  if(!info)
  {
    printf("end test #%i successful\n",i);
  }
  else
  {
    printf("end test #%i  not  successful\n",i);
  }

  i++;
  printf("start test #%i\n",i);
  info += test_2();
  if(!info)
  {
    printf("end test #%i successful\n",i);
  }
  else
  {
    printf("end test #%i  not  successful\n",i);
  }


  i++;
  printf("start test #%i \n",i);
  info += test_3();
  if(!info)
  {
    printf("end test #%i successful\n",i);
  }
  else
  {
    printf("end test #%i  not  successful\n",i);
  }


  i++;
  printf("start test #%i \n",i);
  info += test_4();
  if(!info)
  {
    printf("end test #%i successful\n",i);
  }
  else
  {
    printf("end test #%i  not  successful\n",i);
  }


  return info;



}
