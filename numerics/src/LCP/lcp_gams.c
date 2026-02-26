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

#include "LCP_Solvers.h"       // for lcp_gams
#include "lcp_cst.h"
#include "NumericsFwd.h"       // for LinearComplementarityProblem, SolverOptions
#include "numerics_verbose.h"

#ifdef HAVE_GAMS_C_API

#include <assert.h>
#include <string.h>

#include "GAMSlink.h"
#include "LinearComplementarityProblem.h"
#include "NumericsMatrix.h"
#include "SolverOptions.h"
#endif

void lcp_gams(LinearComplementarityProblem *problem, double *z, double *w, int *info,
              SolverOptions *options) {
#ifndef HAVE_GAMS_C_API
  numerics_error("lcp_gams", "GAMS API is not enabled. Try to re-compile Siconos with GAMS.");
#else

  if (!options->solverParameters)
    options->solverParameters = createGAMSparams(GAMS_MODELS_SHARE_DIR, GAMS_DIR);

  assert(problem);
  assert(problem->size > 0);
  assert(problem->M);
  assert(problem->q);

  /* Handles to the GAMSX, GDX, and Option objects */
  gamsxHandle_t Gptr = NULL;
  idxHandle_t Xptr = NULL;
  optHandle_t Optr = NULL;
  optHandle_t solverOptPtr = NULL;

  int status = 0;
  char sysdir[GMS_SSSIZE], model[GMS_SSSIZE], msg[GMS_SSSIZE];
  const char defModel[] = SPACE_CONC(GAMS_MODELS_SHARE_DIR, "/lcp.gms");
  const char defGAMSdir[] = GAMS_DIR;

  SN_Gams_set_dirs((SN_GAMSparams *)options->solverParameters, defModel, defGAMSdir, model,
                   sysdir, "/lcp.gms");

  /* Create objects */
  if (!gamsxCreateD(&Gptr, sysdir, msg, sizeof(msg))) {
    printf("Could not create gamsx object: %s\n", msg);
    return;
  }

  if (!idxCreateD(&Xptr, sysdir, msg, sizeof(msg))) {
    printf("Could not create gdx object: %s\n", msg);
    return;
  }

  if (!optCreateD(&solverOptPtr, sysdir, msg, sizeof(msg))) {
    printf("Could not create solver opt object: %s\n", msg);
    return;
  }

  if (!optCreateD(&Optr, sysdir, msg, sizeof(msg))) {
    printf("Could not create opt object: %s\n", msg);
    return;
  }

  char gdxFileName[GMS_SSSIZE];
  char solFileName[GMS_SSSIZE];
  const char *base_name = GAMSP_get_filename(options->solverParameters);
  if (base_name) {
    strncpy(gdxFileName, base_name, sizeof(gdxFileName));
  } else {
    strncpy(gdxFileName, "lcp", sizeof(gdxFileName));
  }

  strncpy(solFileName, gdxFileName, sizeof(solFileName));
  strncat(solFileName, "_sol", sizeof(solFileName) - strlen(solFileName) - 1);

  strncat(gdxFileName, ".gdx", sizeof(gdxFileName) - strlen(gdxFileName) - 1);
  strncat(solFileName, ".gdx", sizeof(solFileName) - strlen(solFileName) - 1);

  getGamsOpt(Optr, sysdir);

  optSetStrStr(Optr, "User1", gdxFileName);
  optSetStrStr(Optr, "User2", solFileName);

  getGamsSolverOpt(solverOptPtr, sysdir, "path");
  optSetDblStr(solverOptPtr, "convergence_tolerance", options->dparam[SICONOS_DPARAM_TOL]);
  optWriteParameterFile(solverOptPtr, "path.opt");

  idxOpenWrite(Xptr, gdxFileName, "Siconos/Numerics NM_to_GDX", &status);
  if (status) idxerror(status, "idxOpenWrite");

  if ((status = NM_to_GDX(Xptr, "M", "M matrix", problem->M))) {
    printf("Model data not written\n");
    goto TERMINATE;
  }

  if ((status = NV_to_GDX(Xptr, "q", "q vector", problem->q, problem->size))) {
    printf("Model data not written\n");
    goto TERMINATE;
  }

  if (idxClose(Xptr)) idxerror(idxGetLastError(Xptr), "idxClose");

  if ((status = CallGams(Gptr, Optr, sysdir, model))) {
    printf("Call to GAMS failed\n");
    goto TERMINATE;
  }

  /************************************************
   * Read back solution
   ************************************************/
  idxOpenRead(Xptr, solFileName, &status);
  if (status) idxerror(status, "idxOpenRead");

  if ((status = GDX_to_NV(Xptr, "sol", z, problem->size))) {
    printf("Model data not read\n");
    goto TERMINATE;
  }

  if (idxClose(Xptr)) idxerror(idxGetLastError(Xptr), "idxClose");

TERMINATE:
  optFree(&Optr);
  optFree(&solverOptPtr);
  idxFree(&Xptr);
  gamsxFree(&Gptr);

  *info = status;

#endif  // HAVE_GAMS_C_API
}

static void lcp_gams_set_default(SolverOptions* options) {
  /* No specific defaults needed */
  (void)options;
}

/* ===========================================================================
 * Solver Registration
 * ===========================================================================
 * This registers SICONOS_LCP_GAMS in the global solver registry.
 */

#include "utils/solver_registry.h"
#include "utils/numerics_errors.h"

static int lcp_gams_init_wrap(void* problem, SolverOptions* options) {
  (void)problem;
  lcp_gams_set_default(options);
  return NUMERICS_OK;
}

static int lcp_gams_solve_wrap(void* problem, double* reaction,
                               double* velocity, SolverOptions* options) {
  int info = NUMERICS_OK;
  lcp_gams((LinearComplementarityProblem*)problem, reaction, velocity, &info, options);
  return info;
}

static void lcp_gams_free_wrap(void* problem, SolverOptions* options) {
  (void)problem;
  (void)options;
}

REGISTER_SOLVER(SICONOS_LCP_GAMS, "LCP_GAMS",
                       "GAMS solver for LCP",
                       lcp_gams_init_wrap,
                       lcp_gams_solve_wrap,
                       lcp_gams_free_wrap,
                       NULL,  /* error function */
                       lcp_gams_set_default,  /* set_default */
                       1000,  /* default_max_iter */
                       1e-6,  /* default_tol */
                       0);     /* is_local_solver */
