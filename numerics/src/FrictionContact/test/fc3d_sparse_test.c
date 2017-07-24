/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <errno.h>


#include "NonSmoothDrivers.h"
#include "frictionContact_test_function.h"
#include "SolverOptions.h"
#include "FrictionContactProblem.h"
#include "Friction_cst.h"
#include "fc3d_Solvers.h"
#include "numerics_verbose.h"

#define GET_ITER(X) X.iparam[1] ? X.iparam[1] : X.iparam[7]

static double set_tol(int solver_id, char* problem)
{
  if (solver_id == SICONOS_FRICTION_3D_NSN_FB)
  {
    return 1e-12;
  }
  if ((solver_id == SICONOS_FRICTION_3D_SOCLCP) && !strcmp(problem, "data/NESpheres_30_1.dat"))
  {
    return 1e-10;
  }
  return 100 * DBL_EPSILON;
}

static bool skip(int solver_id, char* problem)
{
  if (problem)
  {
    if ((solver_id == SICONOS_FRICTION_3D_NSGS) && !strcmp(problem, "data/Rover4396.dat"))
    {
      return true;
    }
  }
  return false;
}

static int solve_sparse(int solver_id, FrictionContactProblem* FC, double* r, double* u, char* problem)
{
  if (skip(solver_id, problem)) return 0;
  SolverOptions SO;
  fc3d_setDefaultSolverOptions(&SO, solver_id);
  SO.dparam[0] = set_tol(solver_id, problem);
  int info;
  if (FC)
    info = fc3d_driver(FC, r, u, &SO);
  else
  {
    FILE *finput = fopen(problem, "r");
    if (!finput)
    {
      int _errno = errno;
      fprintf(stderr, "%s :: unable to open file %s\n", __func__, problem);
      return _errno;
    }
    FrictionContactProblem* problem = (FrictionContactProblem*)malloc(sizeof(FrictionContactProblem));
    frictionContact_newFromFile(problem, finput);
    fclose(finput);
    int NC = problem->numberOfContacts;
    int dim = problem->dimension;
    int n = NC*dim;
    double *reaction = (double*)calloc(n, sizeof(double));
    double *velocity = (double*)calloc(n, sizeof(double));

    NumericsMatrix* W = NM_create(NM_SPARSE, n, n);
    NM_copy_to_sparse(problem->M, W);
    NM_free(problem->M);
    free(problem->M);
    problem->M = W;

    info = fc3d_driver(problem, reaction, velocity, &SO);

    free(reaction);
    free(velocity);

    freeFrictionContactProblem(problem);

  }

  if (info)
  {
    fprintf(stderr, "Solver %s (sparse) FAILED with error %d on problem %s. Residual is %e\n", solver_options_id_to_name(solver_id), info, problem, SO.dparam[1]);
  }
  else
  {
    printf("Solver %s (sparse) succeded with %d iterations on problem %s\n", solver_options_id_to_name(solver_id), GET_ITER(SO), problem);
  }

  solver_options_delete(&SO);
  return info;
}

static int solve_dense(int solver_id, FrictionContactProblem* FC, double* r, double* u, char* problem)
{
  if (skip(solver_id, problem)) return 0;
  SolverOptions SO;
  fc3d_setDefaultSolverOptions(&SO, solver_id);
  SO.dparam[0] = set_tol(solver_id, problem);
  int info;
  if (FC)
    info = fc3d_driver(FC, r, u, &SO);
  else
  {
    FILE *finput = fopen(problem, "r");
    if (!finput)
    {
      int _errno = errno;
      fprintf(stderr, "%s :: unable to open file %s\n", __func__, problem);
      return _errno;
    }
    FrictionContactProblem* problem = (FrictionContactProblem*)malloc(sizeof(FrictionContactProblem));
    frictionContact_newFromFile(problem, finput);
    fclose(finput);
    int NC = problem->numberOfContacts;
    int dim = problem->dimension;
    int n = NC*dim;
    double *reaction = (double*)calloc(n, sizeof(double));
    double *velocity = (double*)calloc(n, sizeof(double));

    info = fc3d_driver(problem, reaction, velocity, &SO);

    free(reaction);
    free(velocity);

    freeFrictionContactProblem(problem);
  }

  if (info)
  {
    fprintf(stderr, "Solver %s (dense) FAILED with error %d on problem %s. Residual is %e\n", solver_options_id_to_name(solver_id), info, problem, SO.dparam[1]);
  }
  else
  {
    printf("Solver %s (dense) succeded with %d iterations on problem %s\n", solver_options_id_to_name(solver_id), GET_ITER(SO), problem);
  }
  solver_options_delete(&SO);
  return info;
}

int main(void)
{
  int total_info = 0;

//  numerics_set_verbose(1);

  double q[] = { -1, 1, 3, -1, 1, 3, -1, 1, 3};
  double mu[] = {0.1, 0.1, 0.1};

  double Wdata[81] = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};

  NumericsMatrix* tmpM = NM_create_from_data(NM_DENSE, 9, 9, Wdata);
  NumericsMatrix* W = NM_create(NM_SPARSE, 9, 9);
  NM_copy_to_sparse(tmpM, W);

  int solvers_to_test[] = {SICONOS_FRICTION_3D_NSGS,
                           SICONOS_FRICTION_3D_NSN_AC,
                           SICONOS_FRICTION_3D_NSN_AC_TEST,
                           SICONOS_FRICTION_3D_NSN_FB,
                           SICONOS_FRICTION_3D_NSN_NM,
                           SICONOS_FRICTION_3D_SOCLCP,
                           SICONOS_FRICTION_3D_ACLMFP,
                           SICONOS_FRICTION_3D_PROX,
                           SICONOS_FRICTION_3D_HP,
                           SICONOS_FRICTION_3D_FPP,
                           SICONOS_FRICTION_3D_EG,
                           SICONOS_FRICTION_3D_VI_FPP,
                           SICONOS_FRICTION_3D_VI_EG};

  char* filetests[] = {"data/Rover1039.dat",
                       "data/Rover1040.dat",
                       "data/Rover1041.dat",
                       "data/Rover11035.dat",
                       "data/Rover11211.dat",
                       "data/Rover3865.dat",
                       "data/Rover4144.dat",
                       "data/Rover4396.dat",
                       "data/Rover4493.dat",
                       "data/Rover4516.dat",
                       "data/Rover4609.dat",
                       "data/Rover4613.dat",
                       "data/Rover4622.dat",
                       "data/Rover9770.dat",
                       "data/NESpheres_10_1.dat",
                       "data/GFC3D_OneContact.dat",
                       "data/NESpheres_30_1.dat"};

  for (size_t s = 0; s < sizeof(solvers_to_test)/sizeof(int); ++s)
  {
    int solver_id = solvers_to_test[s];

    FrictionContactProblem* FC = frictionContactProblem_new(3, 3, W, q, mu);
    double r[9] = {0.};
    double u[9] = {0.};

    int infos1 = solve_sparse(solver_id, FC, r, u, "dummy");
    total_info = total_info ? total_info : infos1;
    FC->M = NULL;
    FC->q = NULL;
    FC->mu = NULL;
    freeFrictionContactProblem(FC);

    if (solver_id != SICONOS_FRICTION_3D_HP)
    {
      for (size_t i = 0; i < sizeof(filetests)/sizeof(char*); ++i)
      {
        int infos2 = solve_sparse(solver_id, NULL, NULL, NULL, filetests[i]);
        total_info = total_info ? total_info : infos2;
      }
    }

    FrictionContactProblem* FCdense = frictionContactProblem_new(3, 3, tmpM, q, mu);
    double rdense[9] = {0.};
    double udense[9] = {0.};

    int infod = solve_dense(solver_id, FCdense, rdense, udense, "dummy");
    total_info = total_info ? total_info : infod;

    FCdense->M = NULL;
    FCdense->q = NULL;
    FCdense->mu = NULL;
    freeFrictionContactProblem(FCdense);

    if (solver_id != SICONOS_FRICTION_3D_HP)
    {
      for (size_t i = 0; i < sizeof(filetests)/sizeof(char*); ++i)
      {
        int infod2 = solve_dense(solver_id, NULL, NULL, NULL, filetests[i]);
        total_info = total_info ? total_info : infod2;
      }
    }

  }

  NM_free(W);
  tmpM->matrix0 = NULL;
  NM_free(tmpM);
  free(W);
  free(tmpM);

  return total_info;
}
