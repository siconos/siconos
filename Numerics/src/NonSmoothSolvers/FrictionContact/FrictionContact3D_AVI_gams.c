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
#include "FrictionContactProblem.h"
#include "FrictionContact3D_Solvers.h"
#include "FrictionContact3D_compute_error.h"

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

static int SN_rm_normal_part(int i, int j, double val, void* env)
{
  if (i%3 == 0)
  {
    return 0;
  }
  else
  {
    return 1;
  }
}

static int frictionContact3D_AVI_gams_base(FrictionContactProblem* problem, double *reaction, double *velocity, SolverOptions* options, const char* solverName)
{

  assert(problem);
  assert(problem->numberOfContacts > 0);
  assert(problem->M);
  assert(problem->q);

  /* Handles to the GAMSX, GDX, and Option objects */
  gamsxHandle_t Gptr = NULL;
  idxHandle_t Xptr = NULL;
  optHandle_t Optr = NULL;
  optHandle_t solverOptPtr = NULL;

  int status;
  char sysdir[GMS_SSSIZE], model[GMS_SSSIZE], msg[GMS_SSSIZE];
  const char defModel[] = SPACE_CONC(GAMS_MODELS_SHARE_DIR, "/fc_vi-condensed.gms");
  const char defGAMSdir[] = GAMS_DIR;

  int size = problem->dimension*problem->numberOfContacts;

  NumericsMatrix Wmat;
  NumericsMatrix Emat;
  NumericsMatrix Akmat;

  NM_null(&Wmat);
  NM_null(&Emat);
  NM_null(&Akmat);

  DEBUG_PRINT("FC3D_AVI_GAMS :: seeting basic directories\n");
  SN_Gams_set_dirs(options->solverParameters, defModel, defGAMSdir, model, sysdir, "/fc_vi-condensed.gms");

  /* Create objects */
  DEBUG_PRINT("FC3D_AVI_GAMS :: creating gamsx object\n");
  if (! gamsxCreateD (&Gptr, sysdir, msg, sizeof(msg))) {
    printf("Could not create gamsx object: %s\n", msg);
    return 1;
  }

  DEBUG_PRINT("FC3D_AVI_GAMS :: creating gdx object\n");
  if (! idxCreateD (&Xptr, sysdir, msg, sizeof(msg))) {
    printf("Could not create gdx object: %s\n", msg);
    return 1;
  }

  DEBUG_PRINT("FC3D_AVI_GAMS :: creating opt object\n");
  if (! optCreateD (&Optr, sysdir, msg, sizeof(msg))) {
    printf("Could not create opt object: %s\n", msg);
    return 1;
  }

  DEBUG_PRINT("FC3D_AVI_GAMS :: creating solveropt object\n");
  if (! optCreateD (&solverOptPtr, sysdir, msg, sizeof(msg))) {
    printf("Could not create solveropt object: %s\n", msg);
    return 1;
  }

  getGamsSolverOpt(solverOptPtr, sysdir, solverName);
  optSetDblStr(solverOptPtr, "convergence_tolerance", options->dparam[0]);
//  strncpy(msg, "./", sizeof(deffile));
  strncpy(msg, solverName, sizeof(msg));
  strncat(msg, ".opt", sizeof(msg));
  optWriteParameterFile(solverOptPtr, msg);

  FILE* f = fopen("jams.opt", "w");
  if (f)
  {
    char contents[] = "subsolveropt 1";
    fprintf(f, contents);
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

  idxOpenWrite(Xptr, "fc3d_avi-condensed.gdx", "Siconos/Numerics NM_to_GDX", &status);
  if (status)
    idxerrorR(status, "idxOpenWrite");
  DEBUG_PRINT("FC3D_AVI_GAMS :: fc3d_avi-condensed.gdx opened\n");

  if ((status=NM_to_GDX(Xptr, "W", "W matrix", problem->M))) {
    printf("Model data not written\n");
    goto TERMINATE;
  }
  DEBUG_PRINT("FC3D_AVI_GAMS :: W matrix written\n");

  Wmat.storageType = NM_SPARSE;
  NM_copy_to_sparse(problem->M, &Wmat);
  cs_fkeep(NM_csc(&Wmat), &SN_rm_normal_part, NULL);
  DEBUG_PRINT("FC3D_AVI_GAMS :: Wt matrix constructed\n");

  cblas_dcopy(size, problem->q, 1, reaction, 1);
  for (unsigned i = 0; i < size; i += 3)
  {
    reaction[i] = 0.;
  }

  Emat.storageType = NM_SPARSE;
  NM_sparse(&Emat);
  Emat.size0 = size;
  Emat.size1 = size;

  Emat.matrix2->triplet = cs_spalloc(size, size, problem->numberOfContacts, 1, 1);

  for (int i = 0; i < size; i += 3)
  {
    cs_entry(Emat.matrix2->triplet, i, i, 1.);
  }

  Akmat.storageType = NM_SPARSE;
  NM_sparse(&Akmat);
  Akmat.size0 = NB_APPROX*problem->numberOfContacts;
  Akmat.size1 = size;
  Akmat.matrix2->triplet = cs_spalloc(NB_APPROX*problem->numberOfContacts, size, NB_APPROX*problem->numberOfContacts*3, 1, 1);
  double vec[3];
  double angle = 2*M_PI/(NB_APPROX + 1);
  DEBUG_PRINTF("angle: %g\n", angle);
  unsigned offset_row = 0;
  for (unsigned j = 0; j < problem->numberOfContacts; ++j)
  {
    vec[0] = problem->mu[j];
    for (unsigned i = 0; i < NB_APPROX; ++i)
    {
      vec[1] = cos(i*angle);
      vec[2] = sin(i*angle);
      cs_entry(Akmat.matrix2->triplet, i + offset_row, 3*j, vec[0]);
      cs_entry(Akmat.matrix2->triplet, i + offset_row, 3*j + 1, vec[1]);
      cs_entry(Akmat.matrix2->triplet, i + offset_row, 3*j + 2, vec[2]);
    }
    offset_row += NB_APPROX;
  }

  if ((status=NM_to_GDX(Xptr, "Wt", "Wt matrix", &Wmat))) {
    printf("Model data not written\n");
    goto TERMINATE;
  }
  DEBUG_PRINT("FC3D_AVI_GAMS :: Wt matrix written\n");

  if ((status=NM_to_GDX(Xptr, "E", "E matrix", &Emat))) {
    printf("Model data not written\n");
    goto TERMINATE;
  }
  DEBUG_PRINT("FC3D_AVI_GAMS :: E matrix written\n");

  if ((status=NM_to_GDX(Xptr, "Ak", "Ak matrix", &Akmat))) {
    printf("Model data not written\n");
    goto TERMINATE;
  }
  DEBUG_PRINT("FC3D_AVI_GAMS :: Ak matrix written\n");

  if ((status=NV_to_GDX(Xptr, "q", "q vector", problem->q, size))) {
    printf("Model data not written\n");
    goto TERMINATE;
  }
  DEBUG_PRINT("FC3D_AVI_GAMS :: q vector written\n");

  if ((status=NV_to_GDX(Xptr, "qt", "qt vector", reaction, size))) {
    printf("Model data not written\n");
    goto TERMINATE;
  }
  DEBUG_PRINT("FC3D_AVI_GAMS :: qt vector written\n");

  if (idxClose(Xptr))
    idxerrorR(idxGetLastError(Xptr), "idxClose");

  if ((status=CallGams(Gptr, Optr, sysdir, model))) {
    printf("Call to GAMS failed\n");
    goto TERMINATE;
  }


  /************************************************
   * Read back solution
   ************************************************/
  idxOpenRead(Xptr, "fc3d_avi-condensed_sol.gdx", &status);
  if (status)
    idxerrorR(status, "idxOpenRead");

  /* GAMS does not set a value to 0 ... --xhub */
  memset(reaction, 0, size*sizeof(double));
  if ((status=GDX_to_NV(Xptr, "reaction", reaction, size))) {
    printf("Model data not read\n");
    goto TERMINATE;
  }

  memset(velocity, 0, size*sizeof(double));
  if ((status=GDX_to_NV(Xptr, "velocity", velocity, size))) {
    printf("Model data not read\n");
    goto TERMINATE;
  }

  DEBUG_PRINT_VEC(reaction, size);
  DEBUG_PRINT_VEC(velocity, size);

  if (idxClose(Xptr))
    idxerrorR(idxGetLastError(Xptr), "idxClose");

TERMINATE:
  optFree(&Optr);
  idxFree(&Xptr);
  optFree(&solverOptPtr);
  gamsxFree(&Gptr);
  freeNumericsMatrix(&Wmat);
  freeNumericsMatrix(&Emat);
  freeNumericsMatrix(&Akmat);

  status = FrictionContact3D_compute_error(problem, reaction, velocity, options->dparam[0], options, &(options->dparam[1]));
  DEBUG_PRINTF("FrictionContact3D_AVI_gams :: error = %g\n", options->dparam[1]);
  return status;
}

void frictionContact3D_AVI_gams_path(FrictionContactProblem* problem, double* reaction, double* velocity, int *info, SolverOptions* options)
{
  *info = frictionContact3D_AVI_gams_base(problem, reaction, velocity, options, "path");
}

void frictionContact3D_AVI_gams_pathvi(FrictionContactProblem* problem, double* reaction, double* velocity, int *info, SolverOptions* options)
{
  *info = frictionContact3D_AVI_gams_base(problem, reaction, velocity, options, "pathvi");
}

#else

void frictionContact3D_AVI_gams_path(FrictionContactProblem* problem, double *reaction, double *velocity, int *info, SolverOptions* options)
{
  printf("frictionContact3D_gams :: gams was not enabled at compile time!\n");
  exit(EXIT_FAILURE);
}

void frictionContact3D_AVI_gams_pathvi(FrictionContactProblem* problem, double *reaction, double *velocity, int *info, SolverOptions* options)
{
  printf("frictionContact3D_gams :: gams was not enabled at compile time!\n");
  exit(EXIT_FAILURE);
}
#endif
