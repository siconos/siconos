/*
   Use this command to compile the example:
   cl xp_example2.c api/gdxcc.c api/optcc.c api/gamsxcc.c -Iapi
   */

/*
   This program performs the following steps:
   1. Generate a gdx file with demand data
   2. Calls GAMS to solve a simple transportation model
   (The GAMS model writes the solution to a gdx file)
   3. The solution is read from the gdx file
   */

/* GAMS stuff */

#define _XOPEN_SOURCE 700

#include "NumericsMatrix.h"
#include "NumericsSparseMatrix.h"
#include "GlobalFrictionContactProblem.h"
#include "gfc3d_Solvers.h"

#ifdef HAVE_GAMS_C_API

#include "GAMSlink.h"

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>


#include "sanitizer.h"

#define DEBUG_STDOUT
#define DEBUG_MESSAGES
#include "debug.h"

#define NB_APPROX 10

static ptrdiff_t SN_rm_normal_part_on_H(ptrdiff_t i, ptrdiff_t j, double val, void* env)
{
  if (j%3 == 0)
  {
    return 0;
  }
  else
  {
    return 1;
  }
}

static void FC3D_gams_generate_first_constraints(NumericsMatrix* Akmat, double* mus)
{
  unsigned nb_contacts = (unsigned)Akmat->size1/3;
  assert(nb_contacts*3 == (unsigned)Akmat->size1);
  unsigned nb_approx = (unsigned)Akmat->size0/nb_contacts;
  assert(nb_approx*nb_contacts == (unsigned)Akmat->size0);
  unsigned offset_row = 0;
  CSparseMatrix* triplet_mat = Akmat->matrix2->triplet;

  double angle = 2*M_PI/(NB_APPROX + 1);
  DEBUG_PRINTF("angle: %g\n", angle);

  for (unsigned j = 0; j < nb_contacts; ++j)
  {
    double mu = mus[j];
    for (unsigned i = 0; i < nb_approx; ++i)
    {
      cs_entry(triplet_mat, i + offset_row, 3*j, mu);
      cs_entry(triplet_mat, i + offset_row, 3*j + 1, cos(i*angle));
      cs_entry(triplet_mat, i + offset_row, 3*j + 2, sin(i*angle));
    }
    offset_row += nb_approx;
  }
}

static int gfc3d_AVI_gams_base(GlobalFrictionContactProblem* problem, double *reaction, double *velocity, SolverOptions* options, const char* solverName)
{

  assert(problem);
  assert(problem->numberOfContacts > 0);
  assert(problem->M);
  assert(problem->q);

  /* Handles to the Option objects */
  optHandle_t Optr = NULL;
  optHandle_t solverOptPtr = NULL;

  int status;
  char sysdir[GMS_SSSIZE], model[GMS_SSSIZE], msg[GMS_SSSIZE], template_filename[GMS_SSSIZE], hdf5_filename[GMS_SSSIZE], log_filename[GMS_SSSIZE], base_filename[GMS_SSSIZE];
  const char defModel[] = SPACE_CONC(GAMS_MODELS_SHARE_DIR, "/fc_vi.gms");
  const char defGAMSdir[] = GAMS_DIR;

  int size = problem->dimension*problem->numberOfContacts;

  NumericsMatrix Htmat;
  NM_fill(&Htmat, NM_SPARSE, problem->H->size0, problem->H->size1, NULL);
  NumericsMatrix Emat;
  NM_null(&Emat);
  Emat.storageType = NM_SPARSE;
  NM_sparse(&Emat);
  Emat.size0 = size;
  Emat.size1 = size;

  Emat.matrix2->triplet = cs_spalloc(size, size, problem->numberOfContacts, 1, 1);

  for (unsigned i = 0; i < size; i += 3)
  {
    cs_entry(Emat.matrix2->triplet, i, i, 1.);
  }

  NumericsMatrix Akmat;
  NM_null(&Akmat);
  Akmat.storageType = NM_SPARSE;
  NM_sparse(&Akmat);
  Akmat.size0 = NB_APPROX*problem->numberOfContacts;
  Akmat.size1 = size;
  Akmat.matrix2->triplet = cs_spalloc(NB_APPROX*problem->numberOfContacts, size, NB_APPROX*problem->numberOfContacts*3, 1, 1);
  CSparseMatrix* Ak_triplet = Akmat.matrix2->triplet;

  FC3D_gams_generate_first_constraints(&Akmat, problem->mu);


  SN_Gams_set_dirs((SN_GAMSparams*)options->solverParameters, defModel, defGAMSdir, model, sysdir, "/fc_vi.gms");

  const char* filename = GAMSP_get_filename(options->solverParameters);

  /* Logger starting  */
  if (filename)
  {
    strncpy(hdf5_filename, filename, sizeof(hdf5_filename));
    strncpy(base_filename, filename, sizeof(base_filename));

    const char* suffix = GAMSP_get_filename_suffix(options->solverParameters);

    if (suffix)
    {
      strncat(hdf5_filename, "_", sizeof(hdf5_filename) - strlen(hdf5_filename) - 1);
      strncat(hdf5_filename, suffix, sizeof(hdf5_filename) - strlen(hdf5_filename) - 1);
      strncat(base_filename, "_", sizeof(base_filename) - strlen(base_filename) - 1);
      strncat(base_filename, suffix, sizeof(base_filename) - strlen(base_filename) - 1);
    }
  }
  else
  {
    strncpy(hdf5_filename, "logger", sizeof(hdf5_filename));
    strncpy(base_filename, "output", sizeof(base_filename));
  }
  if (! optCreateD (&Optr, sysdir, msg, sizeof(msg))) {
    printf("Could not create opt object: %s\n", msg);
    return 1;
  }

  if (! optCreateD (&solverOptPtr, sysdir, msg, sizeof(msg))) {
    printf("Could not create opt object: %s\n", msg);
    return 1;
  }

  getGamsSolverOpt(solverOptPtr, sysdir, solverName);
  optSetDblStr(solverOptPtr, "convergence_tolerance", options->dparam[0]);
//  strncpy(msg, "./", sizeof(deffile));
  strncpy(msg, solverName, sizeof(msg));
  strncat(msg, ".opt", sizeof(msg) - strlen(msg) - 1);

  FILE* f = fopen("jams.opt", "w");
  if (f)
  {
    char contents[] = "subsolveropt 1";
    fprintf(f, "%s\n", contents);
    fclose(f);
  }
  else
  {
    printf("Failed to create jams.opt!\n");
  }
  getGamsOpt(Optr, sysdir);
  if (strcmp(solverName, "path"))
  {
    optSetStrStr(Optr, "emp", solverName);
  }

  optHandle_t Opts[] = {Optr, solverOptPtr};
  SN_Gams_set_options((SN_GAMSparams*)options->solverParameters, Opts);

  optWriteParameterFile(solverOptPtr, msg);

  optSetIntStr(Optr, "Keep", 0);

  NM_copy_to_sparse(problem->H, &Htmat);
  cs_fkeep(NM_csc(&Htmat), &SN_rm_normal_part_on_H, NULL);

  cblas_dcopy(size, problem->b, 1, reaction, 1);
  for (unsigned i = 0; i < size; i += 3)
  {
    reaction[i] = 0.;
  }

  SN_GAMS_gdx* gdx_data = (SN_GAMS_gdx*)malloc(sizeof(SN_GAMS_gdx));
  gdx_data->mat_for_gdx = NULL;
  gdx_data->vec_for_gdx = NULL;
  gdx_data->vec_from_gdx = NULL;

  SN_GAMS_add_NM_to_gdx(gdx_data, problem->M, "M");
  SN_GAMS_add_NM_to_gdx(gdx_data, problem->H, "H");
  SN_GAMS_add_NM_to_gdx(gdx_data, &Htmat, "Ht");
  SN_GAMS_add_NM_to_gdx(gdx_data, &Emat, "E");
  SN_GAMS_add_NM_to_gdx(gdx_data, &Akmat, "Ak");
  
  SN_GAMS_add_NV_to_gdx(gdx_data, problem->q, "q", problem->M->size0);
  SN_GAMS_add_NV_to_gdx(gdx_data, problem->b, "b", size);
  SN_GAMS_add_NV_to_gdx(gdx_data, reaction, "bt", size);

   SN_GAMS_add_NV_from_gdx(gdx_data, reaction, "reaction", size);
   SN_GAMS_add_NV_from_gdx(gdx_data, velocity, "velocity", problem->M->size0);

   unsigned iter = 1;
   filename_datafiles(iter, options->solverId, base_filename, sizeof(template_filename), template_filename, log_filename);
   optSetStrStr(Optr, "LogFile", log_filename);
   status = SN_gams_solve(iter, Optr, sysdir, model, template_filename, options, gdx_data);


  SN_free_SN_GAMS_gdx(gdx_data);
  free(gdx_data);
  NM_free(&Htmat);
  NM_free(&Emat);
  NM_free(&Akmat);
  optFree(&Optr);
  optFree(&solverOptPtr);
  return status;
}

void gfc3d_AVI_gams_path(GlobalFrictionContactProblem* problem, double* reaction, double* velocity, int *info, SolverOptions* options)
{
  *info = gfc3d_AVI_gams_base(problem, reaction, velocity, options, "path");
}

void gfc3d_AVI_gams_pathvi(GlobalFrictionContactProblem* problem, double* reaction, double* velocity, int *info, SolverOptions* options)
{
  *info = gfc3d_AVI_gams_base(problem, reaction, velocity, options, "pathvi");
}

#else

void gfc3d_AVI_gams_path(GlobalFrictionContactProblem* problem, double *reaction, double *velocity, int *info, SolverOptions* options)
{
  printf("fc3d_gams :: gams was not enabled at compile time!\n");
  exit(EXIT_FAILURE);
}

void gfc3d_AVI_gams_pathvi(GlobalFrictionContactProblem* problem, double *reaction, double *velocity, int *info, SolverOptions* options)
{
  printf("fc3d_gams :: gams was not enabled at compile time!\n");
  exit(EXIT_FAILURE);
}
#endif
