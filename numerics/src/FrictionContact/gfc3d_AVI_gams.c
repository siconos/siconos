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

#include "NumericsMatrix.h"
#include "GlobalFrictionContactProblem.h"
#include "gfc3d_Solvers.h"

#ifdef HAVE_GAMS_C_API

#include "GAMSlink.h"

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>

#include "sanitizer.h"

#define DEBUG_STDOUT
#define DEBUG_MESSAGES
#include "debug.h"

static ptrdiff_t SN_rm_normal_part(ptrdiff_t i, ptrdiff_t j, double val, void* env)
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

static int gfc3d_AVI_gams_base(GlobalFrictionContactProblem* problem, double *reaction, double *velocity, SolverOptions* options, const char* solverName)
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
  const char defModel[] = SPACE_CONC(GAMS_MODELS_SHARE_DIR, "/fc_vi.gms");
  const char defGAMSdir[] = GAMS_DIR;

  int size = problem->dimension*problem->numberOfContacts;

  NumericsMatrix Htmat;
  fillNumericsMatrix(&Htmat, NM_SPARSE, problem->H->size0, problem->H->size1, NULL);

  SN_Gams_set_dirs(options->solverParameters, defModel, defGAMSdir, model, sysdir, "/fc_vi.gms");

  /* Create objects */
  if (! gamsxCreateD (&Gptr, sysdir, msg, sizeof(msg))) {
    printf("Could not create gamsx object: %s\n", msg);
    return 1;
  }

  if (! idxCreateD (&Xptr, sysdir, msg, sizeof(msg))) {
    printf("Could not create gdx object: %s\n", msg);
    return 1;
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
  optWriteParameterFile(solverOptPtr, msg);

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

  idxOpenWrite(Xptr, "fc3d_avi.gdx", "Siconos/Numerics NM_to_GDX", &status);
  if (status)
    idxerrorR(status, "idxOpenWrite");
  DEBUG_PRINT("GFC3D_AVI_GAMS :: fc3d_avi.gdx opened");

  if ((status=NM_to_GDX(Xptr, "M", "M matrix", problem->M))) {
    printf("Model data not written\n");
    goto TERMINATE;
  }
  DEBUG_PRINT("FC3D_AVI_GAMS :: M matrix written");

  if ((status=NM_to_GDX(Xptr, "H", "H matrix", problem->H))) {
    printf("Model data not written\n");
    goto TERMINATE;
  }
  DEBUG_PRINT("FC3D_AVI_GAMS :: H matrix written");

  NM_copy_to_sparse(problem->H, &Htmat);
  cs_fkeep(NM_csc(&Htmat), &SN_rm_normal_part, NULL);

  cblas_dcopy(size, problem->b, 1, reaction, 1);
  for (unsigned i = 0; i < size; i += 3)
  {
    reaction[i] = 0.;
  }

  if ((status=NM_to_GDX(Xptr, "Ht", "Ht matrix", &Htmat))) {
    printf("Model data not written\n");
    goto TERMINATE;
  }

  if ((status=NV_to_GDX(Xptr, "q", "q vector", problem->q, size))) {
    printf("Model data not written\n");
    goto TERMINATE;
  }

  if ((status=NV_to_GDX(Xptr, "b", "b vector", problem->b, size))) {
    printf("Model data not written\n");
    goto TERMINATE;
  }

  if ((status=NV_to_GDX(Xptr, "bt", "bt vector", reaction, size))) {
    printf("Model data not written\n");
    goto TERMINATE;
  }

  if (idxClose(Xptr))
    idxerrorR(idxGetLastError(Xptr), "idxClose");

  if ((status=CallGams(Gptr, Optr, sysdir, model))) {
    printf("Call to GAMS failed\n");
    goto TERMINATE;
  }


  /************************************************
   * Read back solution
   ************************************************/
  idxOpenRead(Xptr, "fc3d_avi_sol.gdx", &status);
  if (status)
    idxerrorR(status, "idxOpenRead");

  if ((status=GDX_to_NV(Xptr, "reaction", reaction, size))) {
    printf("Model data not read\n");
    goto TERMINATE;
  }

  if ((status=GDX_to_NV(Xptr, "velocities", reaction, size))) {
    printf("Model data not read\n");
    goto TERMINATE;
  }

  if (idxClose(Xptr))
    idxerrorR(idxGetLastError(Xptr), "idxClose");

TERMINATE:
  optFree(&Optr);
  optFree(&solverOptPtr);
  idxFree(&Xptr);
  gamsxFree(&Gptr);
  freeNumericsMatrix(&Htmat);

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
