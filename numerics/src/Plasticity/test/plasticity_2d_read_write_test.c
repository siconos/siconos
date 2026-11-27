// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include <errno.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "NonSmoothDrivers.h"
#include "NumericsMatrix.h"
#include "PlasticityProblem.h"
#include "SolverOptions.h"
#include "numerics_verbose.h"

/* write a problem, read it back and compare all the fields (M as dense) */
static int round_trip(PlasticityProblem* problem, const char* filename) {
  FILE* f = fopen(filename, "w");
  if (!f) return 1;
  plasticity_printInFile(problem, f);
  fclose(f);
  PlasticityProblem* r = plasticity_new_from_filename(filename);
  int n = problem->dimension * problem->numberOfCones;
  int nc = problem->numberOfCones;
  int err = 0;
  err += r->dimension != problem->dimension;
  err += r->numberOfCones != nc;
  err += r->model_type != problem->model_type;
  err += r->stress_coordinates != problem->stress_coordinates;
  for (int i = 0; i < n; i++) err += r->q[i] != problem->q[i];
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      err += NM_get_value(r->M, i, j) != NM_get_value(problem->M, i, j);
  if (problem->model_type == PLASTICITY_MODEL_VON_MISES) {
    err += r->model.von_mises->hypothesis != problem->model.von_mises->hypothesis;
    for (int i = 0; i < nc; i++)
      err += r->model.von_mises->sigma_y[i] != problem->model.von_mises->sigma_y[i];
  } else {
    Plasticity_DruckerPrager_model *a = problem->model.drucker_prager,
                                   *b = r->model.drucker_prager;
    for (int i = 0; i < nc; i++)
      err += (a->eta[i] != b->eta[i]) + (a->theta[i] != b->theta[i]);
    err += (a->cohesion == NULL) != (b->cohesion == NULL);
    if (a->cohesion && b->cohesion)
      for (int i = 0; i < nc; i++) err += a->cohesion[i] != b->cohesion[i];
  }
  printf("round trip %s: %d differences\n", filename, err);
  plasticityProblem_free(r);
  return err ? 1 : 0;
}

/* Von Mises plane stress in Voigt coordinates, Drucker-Prager with cohesion in Mandel
 * coordinates: the optional keyword lines of the file format */
static int test_keywords(void) {
  int info = 0;
  double Md[36] = {0};
  for (int i = 0; i < 6; i++) Md[i + 6 * i] = 2.0 + i;
  Md[1] = Md[6] = 0.5;
  {
    PlasticityProblem* p = plasticityProblem_new();
    p->dimension = 3;
    p->numberOfCones = 2;
    double* Mcopy = (double*)malloc(sizeof(Md));
    memcpy(Mcopy, Md, sizeof(Md));
    p->M = NM_create_with_data(NM_DENSE, 6, 6, Mcopy);
    p->q = (double*)malloc(6 * sizeof(double));
    for (int i = 0; i < 6; i++) p->q[i] = 0.1 * (i - 2.5);
    double* sy = (double*)malloc(2 * sizeof(double));
    sy[0] = 210.;
    sy[1] = 105.;
    p->model_type = PLASTICITY_MODEL_VON_MISES;
    p->model.von_mises = plasticity_VonMises_model_new(sy);
    p->model.von_mises->hypothesis = PLASTICITY_VONMISES_PLANE_STRESS;
    p->stress_coordinates = PLASTICITY_STRESS_COORDINATES_CARTESIAN_VOIGT;
    info += round_trip(p, "plasticity_2d_rw_vonmises_plane_stress.dat");
    plasticityProblem_free(p);
  }
  {
    PlasticityProblem* p = plasticityProblem_new();
    p->dimension = 3;
    p->numberOfCones = 2;
    double* Mcopy = (double*)malloc(sizeof(Md));
    memcpy(Mcopy, Md, sizeof(Md));
    p->M = NM_create_with_data(NM_DENSE, 6, 6, Mcopy);
    p->q = (double*)malloc(6 * sizeof(double));
    for (int i = 0; i < 6; i++) p->q[i] = 0.1 * (i - 2.5);
    double* eta = (double*)malloc(2 * sizeof(double));
    double* theta = (double*)malloc(2 * sizeof(double));
    eta[0] = eta[1] = 0.3;
    theta[0] = theta[1] = 0.1;
    p->model_type = PLASTICITY_MODEL_DRUCKER_PRAGER;
    p->model.drucker_prager = plasticity_DruckerPrager_model_new(eta, theta);
    p->model.drucker_prager->cohesion = (double*)malloc(2 * sizeof(double));
    p->model.drucker_prager->cohesion[0] = 1.5;
    p->model.drucker_prager->cohesion[1] = 2.5;
    p->stress_coordinates = PLASTICITY_STRESS_COORDINATES_CARTESIAN_MANDEL;
    info += round_trip(p, "plasticity_2d_rw_druckerprager_cohesion.dat");
    plasticityProblem_free(p);
  }
  return info;
}

int main(void) {
  int total_info = 0;

  //  numerics_set_verbose(1);

  double q[] = {-1, 1, 3, -1, 1, 3, -1, 1, 3};
  double theta[] = {0.1, 0.1, 0.1};
  double eta[] = {0.1, 0.1, 0.1};

  double Wdata[81] = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                      0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
                      0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};

  NumericsMatrix* tmpM = NM_create_with_data(NM_DENSE, 9, 9, Wdata);
  NumericsMatrix* W = NM_copy_to_sparse(tmpM, 1e-16);

  PlasticityProblem* PLASTICITY_2D = plasticity2DProblem_new_with_data(3, 3, W, q, eta, theta);
  double stress[9] = {0.};
  double plastic_strain_rate[9] = {0.};
  plasticity2D_display(PLASTICITY_2D);
  plasticity2D_printInFilename(PLASTICITY_2D, "plasticity_2d_example1.dat");

  PLASTICITY_2D->M = NULL;
  PLASTICITY_2D->q = NULL;
  PLASTICITY_2D->model.drucker_prager->eta = NULL;
  PLASTICITY_2D->model.drucker_prager->theta = NULL;
  plasticity2DProblem_free(PLASTICITY_2D);

  PlasticityProblem* PLASTICITY_2D_r =
      plasticity2D_new_from_filename("plasticity_2d_example1.dat");
  plasticity2D_display(PLASTICITY_2D_r);
  plasticity2DProblem_free(PLASTICITY_2D_r);

  NM_clear(W);
  tmpM->dense = NULL;
  NM_clear(tmpM);
  free(W);
  free(tmpM);

  total_info += test_keywords();
  return total_info;
}
