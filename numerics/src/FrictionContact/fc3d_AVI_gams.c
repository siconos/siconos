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

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>

#include "SparseMatrix_internal.h"
#include "NumericsMatrix.h"
#include "NumericsSparseMatrix.h"
#include "SolverOptions.h"
#include "FrictionContactProblem.h"
#include "fc3d_Solvers.h"
#include "fc3d_compute_error.h"
#include "projectionOnCone.h"
#include "SiconosCompat.h"

#ifdef HAVE_GAMS_C_API

#include "GAMSlink.h"

#include <math.h>

#include "sanitizer.h"
#include "op3x3.h"

#include "hdf5_logger.h"

#define DEBUG_NOCOLOR
//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "debug.h"

#define NB_APPROX 10

#define TOL_RN 1e-12

#define TOL2 1e-20

#define ETERMINATE 4242

enum { TAKEOFF_CASE, STICKING_CASE, SLIDING_CASE };

#define TOTAL_TIME_USED 2

#define TOTAL_ITER 2
#define LAST_MODEL_STATUS 3
#define LAST_SOLVE_STATUS 4

//#define SMALL_APPROX

#include <fcntl.h>
#include <unistd.h>
#include <errno.h>

static int cp(const char *to, const char *from)
{
    int fd_to, fd_from;
    char buf[4096];
    ssize_t nread;
    int saved_errno;

    fd_from = open(from, O_RDONLY);
    if (fd_from < 0)
        return -1;

    fd_to = open(to, O_WRONLY | O_CREAT | O_EXCL, 0666);
    if (fd_to < 0)
        goto out_error;

    while (nread = read(fd_from, buf, sizeof buf), nread > 0)
    {
        char *out_ptr = buf;
        ssize_t nwritten;

        do {
            nwritten = write(fd_to, out_ptr, nread);

            if (nwritten >= 0)
            {
                nread -= nwritten;
                out_ptr += nwritten;
            }
            else if (errno != EINTR)
            {
                goto out_error;
            }
        } while (nread > 0);
    }

    if (nread == 0)
    {
        if (close(fd_to) < 0)
        {
            fd_to = -1;
            goto out_error;
        }
        close(fd_from);

        /* Success! */
        return 0;
    }

  out_error:
    saved_errno = errno;

    close(fd_from);
    if (fd_to >= 0)
        close(fd_to);

    errno = saved_errno;
    return -1;
}

static inline double rad2deg(double rad) { return rad*180/M_PI; }

static void setDashedOptions(const char* optName, const char* optValue, const char* paramFileName)
{
  FILE* f = fopen(paramFileName, "a");
  if (f)
  {
    fprintf(f, "%s %s\n", optName, paramFileName);
    fclose(f);
  }
  else
  {
    printf("Failed to create option %s with value %s in %s\n", optName, optValue, paramFileName);
  }
}

static CS_INT SN_rm_normal_part(CS_INT i, CS_INT j, double val, void* env)
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

static int FC3D_gams_inner_loop_condensed(unsigned iter, idxHandle_t Xptr, gamsxHandle_t Gptr, optHandle_t Optr, gmoHandle_t gmoPtr, char* sysdir, char* model, const char* base_name, double* restrict reaction, double* restrict velocity, double* restrict tmpq, double* restrict lambda_r, double* restrict lambda_y, NumericsMatrix* W, double* restrict q, NumericsMatrix* Wtmat, NumericsMatrix* Emat, NumericsMatrix* Akmat, SolverOptions* options)
{

  char msg[GMS_SSSIZE];
  int status;
  unsigned size = (unsigned)W->size0;
  double infos[4] = {0.};
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

  DEBUG_PRINT("FC3D_AVI_GAMS :: creating gmo object\n");
  if (! gmoCreateD (&gmoPtr, sysdir, msg, sizeof(msg))) {
    printf("Could not create gmo object: %s\n", msg);
    return 1;
  }

  /* create input and output gdx names*/
  char gdxFileName[GMS_SSSIZE];
  char solFileName[GMS_SSSIZE];
//  char paramFileName[GMS_SSSIZE];

  /* copy the name without extension to creation the   */
//  strncpy(paramFileName, gdxFileName, sizeof(paramFileName));

  strncpy(gdxFileName, base_name, sizeof(gdxFileName));
  strncpy(solFileName, base_name, sizeof(solFileName));
  strncat(solFileName, "_sol", sizeof(solFileName) - strlen(solFileName) - 1);

  strncat(gdxFileName, ".gdx", sizeof(gdxFileName) - strlen(gdxFileName) - 1);
  strncat(solFileName, ".gdx", sizeof(solFileName) - strlen(solFileName) - 1);
//  strncat(paramFileName, ".txt", sizeof(paramFileName));

  /* XXX ParmFile is not a string option */
//  optSetStrStr(Optr, "ParmFile", paramFileName);
//  setDashedOptions("filename", gdxFileName, paramFileName);
   optSetStrStr(Optr, "User1", gdxFileName);
   optSetStrStr(Optr, "User2", solFileName);

  idxOpenWrite(Xptr, gdxFileName, "Siconos/Numerics NM_to_GDX", &status);
  if (status)
    idxerrorR(status, "idxOpenWrite");
  DEBUG_PRINT("FC3D_AVI_GAMS :: fc3d_avi-condensed.gdx opened\n");

  if ((status=NM_to_GDX(Xptr, "W", "W matrix", W))) {
    printf("Model data not written\n");
    infos[1] = (double)-ETERMINATE;
    goto fail;
  }
  DEBUG_PRINT("FC3D_AVI_GAMS :: W matrix written\n");


  cblas_dcopy_msan(size, q, 1, tmpq, 1);
  for (unsigned i = 0; i < size; i += 3)
  {
    tmpq[i] = 0.;
  }


  if ((status=NM_to_GDX(Xptr, "Wt", "Wt matrix", Wtmat))) {
    printf("Model data not written\n");
    infos[1] = (double)-ETERMINATE;
    goto fail;
  }
  DEBUG_PRINT("FC3D_AVI_GAMS :: Wt matrix written\n");

  if ((status=NM_to_GDX(Xptr, "E", "E matrix", Emat))) {
    printf("Model data not written\n");
    infos[1] = (double)-ETERMINATE;
    goto fail;
  }
  DEBUG_PRINT("FC3D_AVI_GAMS :: E matrix written\n");

  if ((status=NM_to_GDX(Xptr, "Ak", "Ak matrix", Akmat))) {
    printf("Model data not written\n");
    infos[1] = (double)-ETERMINATE;
    goto fail;
  }
  DEBUG_PRINT("FC3D_AVI_GAMS :: Ak matrix written\n");

  if ((status=NV_to_GDX(Xptr, "q", "q vector", q, size))) {
    printf("Model data not written\n");
    infos[1] = (double)-ETERMINATE;
    goto fail;
  }
  DEBUG_PRINT("FC3D_AVI_GAMS :: q vector written\n");

  if ((status=NV_to_GDX(Xptr, "qt", "qt vector", tmpq, size))) {
    printf("Model data not written\n");
    infos[1] = (double)-ETERMINATE;
    goto fail;
  }
  DEBUG_PRINT("FC3D_AVI_GAMS :: qt vector written\n");

  if ((status=NV_to_GDX(Xptr, "guess_r", "guess for r", reaction, size))) {
    printf("Model data not written\n");
    infos[1] = (double)-ETERMINATE;
    goto fail;
  }
  DEBUG_PRINT("FC3D_AVI_GAMS :: guess_r vector written\n");

  if ((status=NV_to_GDX(Xptr, "guess_y", "guess for y", velocity, size))) {
    printf("Model data not written\n");
    infos[1] = (double)-ETERMINATE;
    goto fail;
  }
  DEBUG_PRINT("FC3D_AVI_GAMS :: guess_y vector written\n");

  if ((status=NV_to_GDX(Xptr, "guess_lambda_r", "guess for lambda_r", lambda_r, Akmat->size0))) {
    printf("Model data not written\n");
    infos[1] = (double)-ETERMINATE;
    goto fail;
  }
  DEBUG_PRINT("FC3D_AVI_GAMS :: lambda_r vector written\n");

  if ((status=NV_to_GDX(Xptr, "guess_lambda_y", "guess for lambda_y", lambda_y, Akmat->size0))) {
    printf("Model data not written\n");
    infos[1] = (double)-ETERMINATE;
    goto fail;
  }
  DEBUG_PRINT("FC3D_AVI_GAMS :: lambda_y vector written\n");

  if (idxClose(Xptr))
    idxerrorR(idxGetLastError(Xptr), "idxClose");

   cp(gdxFileName, "fc3d_avi-condensed.gdx");

  if ((status=CallGams(Gptr, Optr, sysdir, model))) {
    printf("Call to GAMS failed\n");
    infos[1] = (double)-ETERMINATE;
    goto fail;
  }


  /************************************************
   * Read back solution
   ************************************************/
  idxOpenRead(Xptr, solFileName, &status);
  if (status)
    idxerrorR(status, "idxOpenRead");

  /* GAMS does not set a value to 0 ... --xhub */
  memset(reaction, 0, size*sizeof(double));
  if ((status=GDX_to_NV(Xptr, "reaction", reaction, size))) {
    printf("Model data not read\n");
    infos[1] = (double)-ETERMINATE;
    goto fail;
  }

  memset(velocity, 0, size*sizeof(double));
  if ((status=GDX_to_NV(Xptr, "velocity", velocity, size))) {
    printf("Model data not read\n");
    infos[1] = (double)-ETERMINATE;
    goto fail;
  }

  if ((status=GDX_to_NV(Xptr, "infos", infos, sizeof(infos)/sizeof(double)))) {
    printf("Model data not read\n");
    infos[1] = (double)-ETERMINATE;
    goto fail;
  }

  if (idxClose(Xptr))
    idxerrorR(idxGetLastError(Xptr), "idxClose");

  options->iparam[TOTAL_ITER] += (int)infos[2];
  options->iparam[LAST_MODEL_STATUS] = (int)infos[0];
  options->iparam[LAST_SOLVE_STATUS] = (int)infos[1];
  options->dparam[TOTAL_TIME_USED] += infos[3];
  printf("SolveStat = %d, ModelStat = %d\n", (int)infos[1], (int)infos[0]);
  gmoGetModelStatusTxt(gmoPtr, (int)infos[0], msg);
  DEBUG_PRINTF("%s\n", msg);
  gmoGetSolveStatusTxt(gmoPtr, (int)infos[1], msg);
  DEBUG_PRINTF("%s\n", msg);

fail:
  idxFree(&Xptr);
  gamsxFree(&Gptr);
  gmoFree(&gmoPtr);
  return (int)infos[1];
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

static int fc3d_AVI_gams_base(FrictionContactProblem* problem, double *reaction, double *velocity, SolverOptions* options, const char* solverName)
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
  gmoHandle_t gmoPtr = NULL;

  /* Cleanup previous value */
  options->dparam[TOTAL_TIME_USED] = 0.;
  options->iparam[TOTAL_ITER] = 0;

  int status;
  char sysdir[GMS_SSSIZE], model[GMS_SSSIZE], msg[GMS_SSSIZE], template_filename[GMS_SSSIZE], hdf5_filename[GMS_SSSIZE], log_filename[GMS_SSSIZE], base_filename[GMS_SSSIZE];
  const char defModel[] = SPACE_CONC(GAMS_MODELS_SHARE_DIR, "/fc_vi-condensed.gms");
  const char defGAMSdir[] = GAMS_DIR;

  unsigned size = (unsigned) problem->dimension*problem->numberOfContacts;

  NumericsMatrix Wtmat;
  NumericsMatrix Emat;
  NumericsMatrix Akmat;

  NM_null(&Wtmat);
  NM_null(&Emat);
  NM_null(&Akmat);

  DEBUG_PRINT("FC3D_AVI_GAMS :: seeting basic directories\n");
  SN_Gams_set_dirs((SN_GAMSparams*)options->solverParameters, defModel, defGAMSdir, model, sysdir, "/fc_vi-condensed.gms");

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

  strncat(hdf5_filename, ".hdf5", sizeof(hdf5_filename) - strlen(hdf5_filename) - 1);


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

  optHandle_t Opts[] = {Optr, solverOptPtr};
  SN_Gams_set_options((SN_GAMSparams*)options->solverParameters, Opts);

  if (strcmp(solverName, "path"))
  {
    optSetStrStr(Optr, "emp", solverName);
/*    optSetStrStr(solverOptPtr, "avi_start", "ray_first");
    optSetStrStr(solverOptPtr, "ratio_tester", "expand");
    optSetDblStr(solverOptPtr, "expand_eps", 0.);*/
//    optSetIntStr(solverOptPtr, "scheduler_decompose", 1);
//    optSetStrStr(solverOptPtr, "lemke_factorization_method", "minos_blu");
  }
  else // only for path
  {
/*     optSetIntStr(solverOptPtr, "linear_model_perturb", 0);
    optSetDblStr(solverOptPtr, "proximal_perturbation", 0.);
    optSetStrStr(solverOptPtr, "crash_method", "none");
    optSetIntStr(solverOptPtr, "crash_perturb", 0);
    optSetIntStr(solverOptPtr, "restart_limit", 0);
    optSetStrStr(solverOptPtr, "lemke_start", "first"); */
//    optSetIntStr(solverOptPtr, "output_linear_model", 1);
//    optSetIntStr(solverOptPtr, "output_minor_iterations_frequency", 1);
//    optSetIntStr(solverOptPtr, "output_linear_model", 1);

  }
/*  optSetIntStr(solverOptPtr, "minor_iteration_limit", 100000);
  optSetIntStr(solverOptPtr, "major_iteration_limit", 20);
  optSetDblStr(solverOptPtr, "expand_delta", 1e-10);*/
//  optSetDblStr(solverOptPtr, "convergence_tolerance", 1e-12);
  optSetDblStr(solverOptPtr, "convergence_tolerance", options->dparam[0]);

  optWriteParameterFile(solverOptPtr, msg);

  optSetIntStr(Optr, "Keep", 0);

  Wtmat.storageType = NM_SPARSE;
  NM_copy_to_sparse(problem->M, &Wtmat);
  cs_fkeep(NM_csc(&Wtmat), &SN_rm_normal_part, NULL);
  DEBUG_PRINT("FC3D_AVI_GAMS :: Wt matrix constructed\n");

  Emat.storageType = NM_SPARSE;
  NM_sparse(&Emat);
  Emat.size0 = size;
  Emat.size1 = size;

  Emat.matrix2->triplet = cs_spalloc(size, size, problem->numberOfContacts, 1, 1);

  for (unsigned i = 0; i < size; i += 3)
  {
    cs_entry(Emat.matrix2->triplet, i, i, 1.);
  }

  Akmat.storageType = NM_SPARSE;
  NM_sparse(&Akmat);
  Akmat.size0 = NB_APPROX*problem->numberOfContacts;
  Akmat.size1 = size;
  Akmat.matrix2->triplet = cs_spalloc(NB_APPROX*problem->numberOfContacts, size, NB_APPROX*problem->numberOfContacts*3, 1, 1);
  CSparseMatrix* Ak_triplet = Akmat.matrix2->triplet;

  FC3D_gams_generate_first_constraints(&Akmat, problem->mu);

  double* tmpq = (double*)malloc(size * sizeof(double));
  unsigned size_l = (NB_APPROX+3)*problem->numberOfContacts;
  double* lambda_r = (double*)calloc(size_l, sizeof(double));
  double* lambda_y = (double*)calloc(size_l, sizeof(double));
  double* reaction_old = (double*)calloc(size, sizeof(double));
  double* velocity_old = (double*)calloc(size, sizeof(double));

  double* predicted_angles = (double*)calloc(problem->numberOfContacts, sizeof(double));
  double* delta_angles = (double*)calloc(problem->numberOfContacts, sizeof(double));
  double* real_angles = (double*)calloc(problem->numberOfContacts, sizeof(double));
  double* residual_contact = (double*)calloc(problem->numberOfContacts, sizeof(double));

  /* save what is the current solution:
   * - 0 => r = 0
   * - 1 => r ∈ int K
   * - 2 => r ∈ bdry K \ {0} 
   */
  size_t* type_contact = (size_t*)calloc(problem->numberOfContacts, sizeof(size_t));
  size_t* type_contact_avi = (size_t*)calloc(problem->numberOfContacts, sizeof(size_t));

  bool done = false;
  double total_residual = 0.;
  double old_residual = 1e20;

  double current_nb_approx = NB_APPROX;

  unsigned iter = 0;
  unsigned maxiter = options->iparam[0];

  SN_logh5* logger_s = SN_logh5_init(hdf5_filename, maxiter);

  while (!done && (iter < maxiter))
  {
    iter++;
    total_residual = 0.;
    filename_datafiles(iter, options->solverId, base_filename, sizeof(template_filename), template_filename, log_filename);
    optSetStrStr(Optr, "LogFile", log_filename);

    SN_logh5_new_iter(iter, logger_s);

    SN_logh5_scalar_double(old_residual, "old_residual", logger_s->group);
    SN_logh5_vec_double(size, reaction, "reaction_guess", logger_s->group);
    SN_logh5_vec_double(size, velocity, "velocity_guess", logger_s->group);

    SN_logh5_vec_double(Akmat.size0, lambda_r, "lambda_r", logger_s->group);
    SN_logh5_vec_double(Akmat.size0, lambda_y, "lambda_y", logger_s->group);


    int solverStat = FC3D_gams_inner_loop_condensed(iter, Xptr, Gptr, Optr, gmoPtr, sysdir, model, template_filename, reaction, velocity, tmpq, lambda_r, lambda_y, problem->M, problem->q, &Wtmat, &Emat, &Akmat, options);
    //DEBUG_PRINT_VEC(reaction, size);
    //DEBUG_PRINT_VEC(velocity, size);

    // CSC form should be OK
    SN_logh5_NM(&Akmat, "Akmat", logger_s);
    SN_logh5_vec_double(size, reaction, "reaction_approx", logger_s->group);
    SN_logh5_vec_double(size, velocity, "velocity_approx", logger_s->group);

    SN_logh5_vec_double(problem->numberOfContacts, predicted_angles, "predicted_angles", logger_s->group);

    switch (solverStat)
    {
      case -ETERMINATE:
        {
          goto TERMINATE;
        }
      case gmoSolveStat_Normal:
        {
          /* We are ok here */
          break;
        }
      case gmoSolveStat_Iteration:
        {
          if (verbose > 0)
          {
            printf("Solver failed due to too many iteration\n");
          }
          break;
        }
      case gmoSolveStat_Resource:
      case gmoSolveStat_Solver:
      case gmoSolveStat_User:
      default:
        {
          printf("Unknown Solve Stat return by the solver! Exiting ...\n");
          options->dparam[1] = 1e20;
          goto TERMINATE;
        }
    }

    /************************************************
     * Project on the cone
     ************************************************/

    for (unsigned i3 = 0, i = 0; i3 < size; ++i, i3 += 3)
    {
      /************************************************
       * save contact type based on result from AVI
       ************************************************/
      if (reaction[i3] < TOL_RN)
      {
        type_contact_avi[i] = TAKEOFF_CASE;
      }
      else if (velocity[i3] < 1e-9)
      {
        type_contact_avi[i] = STICKING_CASE;
      }
      else
      {
        type_contact_avi[i] = SLIDING_CASE;
      }

      double mu = problem->mu[i];
      /* Step 1. project r on the cone */

//      DEBUG_PRINTF("contact %d, before projection: theta_r = %.*e\n", i, DECIMAL_DIG, rad2deg(atan2(reaction[i3+2], reaction[i3+1])));
      unsigned proj_type = projectionOnCone(&reaction[i3], mu);
//      DEBUG_PRINTF("contact %d, after projection:  theta_r = %.*e\n", i, DECIMAL_DIG, rad2deg(atan2(reaction[i3+2], reaction[i3+1])));

      double* lri = &reaction[i3];
      if (reaction[i3] < TOL_RN)
      {
        type_contact[i] = TAKEOFF_CASE;
      }
      else if (type_contact_avi[i] == SLIDING_CASE)
      {
        type_contact[i] = SLIDING_CASE;
      }
      else if (proj_type == PROJCONE_INSIDE)
      {
        type_contact[i] = STICKING_CASE;
      }
      else
      {
        assert(proj_type == PROJCONE_BOUNDARY);
        type_contact[i] = SLIDING_CASE;
      }

    }

    memset(lambda_r, 0, size_l * sizeof(double));
    memset(lambda_y, 0, size_l * sizeof(double));

    /************************************************
     * (Re)compute the local velocities
     ************************************************/
    cblas_dcopy(size, problem->q, 1, velocity, 1);
    NM_gemv(1., problem->M, reaction, 1., velocity);

    SN_logh5_vec_double(size, reaction, "reaction_proj", logger_s->group);
    SN_logh5_vec_double(size, velocity, "velocity_proj", logger_s->group);

    /*  BEGIN CLEANUP */
    memset(predicted_angles, 0, problem->numberOfContacts * sizeof(double));
    memset(delta_angles, 0, problem->numberOfContacts * sizeof(double));
    memset(real_angles, 0, problem->numberOfContacts * sizeof(double));
    /* TODO we should zero out the storage and the content, but not deallocate the matrix ! */
    NM_clearSparseStorage(&Akmat);
    /* This is now a upper bound ...  */
    Akmat.matrix2->triplet = cs_spalloc(NB_APPROX*problem->numberOfContacts, size, NB_APPROX*problem->numberOfContacts*3, 1, 1);
    Ak_triplet = Akmat.matrix2->triplet;
    /************************************************
     * Compute the error on each contact point + 
     ************************************************/
    unsigned offset_row = 0;
    double* xtmp = (double*)calloc(size, sizeof(double));
    double workTmp[3];
    for (unsigned i3 = 0, i = 0; i3 < size; ++i, i3 += 3)
    {
      double res = 0.;
      /* Step 2. recompute the local velocities and  */
      double* ri = &reaction[i3];
      double* ui = &velocity[i3];
      DEBUG_PRINTF("Contact %d, old r = [%.*e; %.*e; %.*e]\n", i, DECIMAL_DIG, reaction_old[i3+0], DECIMAL_DIG, reaction_old[i3+1], DECIMAL_DIG, reaction_old[i3+2]);
      DEBUG_PRINTF("Contact %d, new r = [%.*e; %.*e; %.*e]\n", i, DECIMAL_DIG, ri[0], DECIMAL_DIG, ri[1], DECIMAL_DIG, ri[2]);
      DEBUG_PRINTF("Contact %d, del r = [%.*e; %.*e; %.*e]\n", i, DECIMAL_DIG, reaction_old[i3+0]-ri[0], DECIMAL_DIG, reaction_old[i3+1]-ri[1], DECIMAL_DIG, reaction_old[i3+2]-ri[2]);
      assert(i < (unsigned)problem->numberOfContacts);
      double mu = problem->mu[i];
      fc3d_unitary_compute_and_add_error(ri, ui, mu, &res, workTmp);
      residual_contact[i] = sqrt(res);
      DEBUG_EXPR_WE(if (res > old_residual) { printf("Contact %d, res = %g > %g = old_residual\n", i, sqrt(res), old_residual); });
      total_residual += res;
      /* TODO we may want to revisit this, since err < TOL2 should be enough to
       * really reduce the number of constraints ...*/
      /* Well we do not want to mess with the sliding case ( both r and u on
       * the boundaries)*/
      //if ((res < TOL2) && ((ri[0] < TOL_RN) || ((ri[1]*ri[1] + ri[2]*ri[2]) < (1.-10*DBL_EPSILON)*mu*mu * ri[0]*ri[0])))
      if (false)
      {
        DEBUG_PRINTF("Contact %d, res = %g\n", i, sqrt(res));
        DEBUG_EXPR_WE(if (ri[0] < TOL_RN) { printf("ri[0] = %g < %g = tol", ri[0], TOL_RN); });
        DEBUG_EXPR_WE(if ((ri[1]*ri[1] + ri[2]*ri[2]) < mu*mu * ri[0]*ri[0]) { printf("||r_t||^2 = %g < %g = mu^2 r_n^2; diff = %g\n", (ri[1]*ri[1] + ri[2]*ri[2]), mu*mu * ri[0]*ri[0], (ri[1]*ri[1] + ri[2]*ri[2])-(mu*mu * ri[0]*ri[0]));});
      /* 3 hyperplanes (because we don't want a lineality space for now */
        cs_entry(Ak_triplet, offset_row, i3, mu);
        cs_entry(Ak_triplet, offset_row, i3 + 1, 1.);
        cs_entry(Ak_triplet, offset_row, i3 + 2, 0.);

        offset_row++;

        cs_entry(Ak_triplet, offset_row, i3, mu);
        cs_entry(Ak_triplet, offset_row, i3 + 1, -.5);
        cs_entry(Ak_triplet, offset_row, i3 + 2, M_SQRT2/2);

        offset_row++;

        cs_entry(Ak_triplet, offset_row, i3, mu);
        cs_entry(Ak_triplet, offset_row, i3 + 1, -.5);
        cs_entry(Ak_triplet, offset_row, i3 + 2, -M_SQRT2/2);

        offset_row++;

        /* Leave ri as-is  */
      }
//      else if (ri[0] > TOL_RN) // if r= 0 :(
      else if (type_contact[i] == SLIDING_CASE)
      {
        //double delta_angle = atan2(-ri[1]*ui[2] + ui[1]*ri[2], ri[1]*ri[2] + ui[1]*ui[2]);
        double minus_r_angle = atan2(ri[2], ri[1]);
        if (minus_r_angle >= 0.)
        {
          minus_r_angle -= M_PI;
        }
        else if (minus_r_angle < 0.)
        {
          minus_r_angle += M_PI;
        }
        double delta_angle = atan2(ui[2], ui[1]) - minus_r_angle;
        if (delta_angle > M_PI)
        {
          delta_angle -= 2*M_PI;
        }
        else if (delta_angle < -M_PI)
        {
          delta_angle += 2*M_PI;
        }
        delta_angles[i] = rad2deg(delta_angle);
        real_angles[i] = rad2deg(-minus_r_angle);
        if ((fabs(delta_angle) > M_PI/2))
        {
          printf("Contact %d, something bad happened, angle value is %g (rad) or %g (deg)\n", i, delta_angle, rad2deg(delta_angle));
          printf("r = [%g; %g; %g]\tu = [%g; %g; %g]\tres = %g\n", ri[0], ri[1], ri[2], ui[0], ui[1], ui[2], sqrt(res));
          if (((ri[1]*ri[1] + ri[2]*ri[2]) < mu*mu * ri[0]*ri[0]))
          {
            printf("r is the in the interior of the cone ... |r_t| = %g < %g = r_n*mu\n", sqrt((ri[1]*ri[1] + ri[2]*ri[2])), sqrt(mu*mu * ri[0]*ri[0]));
          }
          goto bad_angle;
        }

        /* Ok so here we support that we are in the sliding case */
        DEBUG_PRINTF("contact %d, delta angle = %g, theta_r = %.*e, theta_u = %.*e\n", i, rad2deg(delta_angle), DECIMAL_DIG, rad2deg(atan2(ri[2],ri[1])), 
            DECIMAL_DIG, rad2deg(atan2(ui[2], ui[1])));
        if (fabs(delta_angle) < 1e-12) { printf("Contact %d, delta_angle too small %g; set to 1e-12", i, delta_angle); delta_angle = copysign(1e-12, delta_angle);}

        /* now compute minus the angle, since we want to compute the constraints  */
        unsigned p;
#ifdef SMALL_APPROX
        p = NB_APPROX-1;
#else
        p = NB_APPROX-2;
#endif
        double slice_angle = delta_angle/(p-1);
        DEBUG_PRINTF("contact %d, slice_angle = %g\n", i, rad2deg(slice_angle));

        /* Generate the supporting hyperplanes  */
        double angle = minus_r_angle;
        unsigned offset_row_bck = offset_row;
        for (unsigned row_indx = offset_row; row_indx < p + offset_row; ++row_indx)
        {
          DEBUG_PRINTF("contact %d, row entry %d, picking a point at the angle %g\n", i, row_indx, rad2deg(angle));
          cs_entry(Ak_triplet, row_indx, i3, mu);
          cs_entry(Ak_triplet, row_indx, i3 + 1, cos(angle));
          cs_entry(Ak_triplet, row_indx, i3 + 2, sin(angle));
          angle += slice_angle;
        }
        angle -= slice_angle;
        //assert(fabs(angle - minus_r_angle - angle) < 1e-12);
        if(fabs(angle - minus_r_angle - delta_angle) > 1e-12)
        {
          printf("warning, big difference betwwen original angle and result: %g !\n", angle - minus_r_angle - delta_angle);
        }

        /* Add the last constraint  */
        /* XXX we already have computed those ...  --xhub */
        double middle_point[] = {(cos(minus_r_angle) + cos(angle))/2., (sin(minus_r_angle) + sin(angle))/2.};

#ifdef SMALL_APPROX
        cs_entry(Ak_triplet, p + offset_row, i3, -mu*(hypot(middle_point[0], middle_point[1])));
        cs_entry(Ak_triplet, p + offset_row, i3 + 1, -cos((minus_r_angle+angle)/2));
        cs_entry(Ak_triplet, p + offset_row, i3 + 2, -sin((minus_r_angle+angle)/2));
        /* Update the row index */
        offset_row += p + 1;
#else

        cs_entry(Ak_triplet, p + offset_row, i3, 0.);
        cs_entry(Ak_triplet, p + offset_row+1, i3, 0.);
        if (delta_angle > 0) /* We need to rotate - pi/2 the original angle */
        {
          cs_entry(Ak_triplet, p + offset_row, i3 + 1, cos(minus_r_angle-M_PI/2)); /* XXX there are formulas for this ... */
          cs_entry(Ak_triplet, p + offset_row, i3 + 2, sin(minus_r_angle-M_PI/2));
          cs_entry(Ak_triplet, p + offset_row+1, i3 + 1, cos(angle + M_PI/2));
          cs_entry(Ak_triplet, p + offset_row+1, i3 + 2, sin(angle + M_PI/2));
        }
        else /* We need to rotate of pi/2 */
        {
          cs_entry(Ak_triplet, p + offset_row, i3 + 1, cos(minus_r_angle + M_PI/2));
          cs_entry(Ak_triplet, p + offset_row, i3 + 2, sin(minus_r_angle + M_PI/2));
          cs_entry(Ak_triplet, p + offset_row+1, i3 + 1, cos(angle - M_PI/2));
          cs_entry(Ak_triplet, p + offset_row+1, i3 + 2, sin(angle - M_PI/2));
        }

        offset_row += p + 2;
#endif

        /* update ri =   */
        unsigned pos = (p-1)/2;
        double angle_pos = minus_r_angle + pos*slice_angle;
        predicted_angles[i] = rad2deg(-angle_pos);
        DEBUG_PRINTF("old r[i] = %.*e %.*e %.*e\n", DECIMAL_DIG, ri[0], DECIMAL_DIG, ri[1], DECIMAL_DIG, ri[2]);
        ri[1] = -ri[0]*mu*cos(minus_r_angle + pos*slice_angle);
        ri[2] = -ri[0]*mu*sin(minus_r_angle + pos*slice_angle);
        DEBUG_PRINTF("new r[i] = %g %g %g\n", ri[0], ri[1], ri[2]);
        DEBUG_PRINTF("new \bar{r}[i] = %g %g %g\n", ri[0]*mu/ri[0], ri[1]*mu/ri[0], ri[2]*mu/ri[0]);

//        lambda_r[offset_row_bck + pos] = -(ui[0]*ri[0] + ui[1]*ri[1] + ui[2]*ri[2])/sqrt(1 + mu*mu);
        lambda_r[offset_row_bck + pos] = sqrt(ui[1]*ui[1] + ui[2]*ui[2]);

        lambda_y[offset_row_bck + pos] = ((-ui[1])*ri[1] + (-ui[2])*ri[2])/(ri[0]*mu);
        DEBUG_PRINTF("contact %d, lambda_r = %g; lambda_y = %g\n", i, lambda_r[offset_row_bck + pos], lambda_y[offset_row_bck + pos]);
        DEBUG_PRINTF("contact %d, ui = [%g; %g; %g]\n", i, ui[0], ui[1], ui[2]);
        DEBUG_PRINTF("contact %d, ui_res = [%g; %g; %g]\n", i, -lambda_r[offset_row_bck + pos]*mu + (ui[0] + mu*lambda_y[offset_row_bck + pos]), -lambda_r[offset_row_bck + pos]*cos(angle_pos) + ui[1], -lambda_r[offset_row_bck + pos]*sin(angle_pos) + ui[2]);
        double norm_rt = ri[0]*mu;
        DEBUG_PRINTF("contact %d, yi_res = [%g; %g; %g]\n", i, -lambda_y[offset_row_bck + pos]*mu + ((-ui[1])*ri[1] + (-ui[2])*ri[2])/ri[0], -lambda_y[offset_row_bck + pos]*cos(angle_pos) + ui[1], -lambda_y[offset_row_bck + pos]*sin(angle_pos) + ui[2]);
        /* now we write the estimate for y here */
        ui[0] = ((-ui[1])*ri[1] + (-ui[2])*ri[2])/ri[0];
        ui[1] = ui[0]*ri[1]/ri[0];
        ui[2] = ui[0]*ri[2]/ri[0];
        DEBUG_PRINTF("new y[i] = %g %g %g\n", ui[0], ui[1], ui[2]);
        DEBUG_PRINTF("new \bar{y}[i] = %g %g %g\n", ui[0]*mu/ui[0], ui[1]*mu/ui[0], ui[2]*mu/ui[0]);
        DEBUG_PRINTF(" ||y_t|| = %g <= %g = mu*y_n; diff = %g\n", sqrt(ui[1]*ui[1] + ui[2]*ui[2]), ui[0]*mu, ui[0]*mu-sqrt(ui[1]*ui[1] + ui[2]*ui[2]));
        DEBUG_PRINTF("contact %d, row entry %d, last_entry: angle = %g, norm = %g; coeff = %g, %g, %g\n", i, offset_row, rad2deg(atan2(middle_point[1], middle_point[0])), hypot(middle_point[0], middle_point[1]), -mu*hypot(middle_point[0], middle_point[1]), -middle_point[0], -middle_point[1]);
        double xx[] = {1, -middle_point[0]/((1+hypot(middle_point[0], middle_point[1]))/2), -middle_point[1]/((hypot(middle_point[0], middle_point[1]) + 1)/2)};
        xtmp[i3] = 1./mu; xtmp[i3+1] = xx[1]; xtmp[i3+2] = xx[2];
/*         Akmat.size0 = offset_row + p + 1;
        DEBUG_PRINTF("contact %d, checking feasibility, Ak size = (%d,%d)\n", i, Akmat.size0, Akmat.size1);
        double* xtmp2 = (double*)calloc(Akmat.size0, sizeof(double));
        double* pts = (double*)calloc(size, sizeof(double));
        xtmp[i3] = 1./mu;
        xtmp[i3+1] = xx[1];
        xtmp[i3+2] = xx[2];
        DEBUG_PRINT_VEC(xtmp, size);
        cs_gaxpy(NM_csc(&Akmat), xtmp, xtmp2);
        NM_clearCSC(&Akmat);
        DEBUG_PRINT_VEC(xtmp2, Akmat.size0);
        pts[i3] = mu;
        pts[i3+1] = 0;
        pts[i3+2] = 0;
        DEBUG_PRINT_VEC(pts, size);
        memset(xtmp2, 0, size*sizeof(double));
        cs_gaxpy(NM_csc(&Akmat), pts, xtmp2);
        DEBUG_PRINT_VEC(xtmp2, Akmat.size0);
        NM_clearCSC(&Akmat);
        pts[i3] = 1./mu;
        pts[i3+1] = -xx[1];
        pts[i3+2] = -xx[2];
        printf("%g\n", hypot(middle_point[0], middle_point[1]));
        DEBUG_PRINT_VEC(pts, size);
        memset(xtmp2, 0, size*sizeof(double));
        cs_gaxpy(NM_csc(&Akmat), pts, xtmp2);
        DEBUG_PRINT_VEC(xtmp2, Akmat.size0);
        NM_clearCSC(&Akmat);
        NM_display(&Akmat);
        DEBUG_PRINTF("contact %d, checking feasibility\n", i);
        free(xtmp2);
        free(pts);*/


      }
      else // r = 0, or r in int(cone) but other interactions moved u
bad_angle:
      {
        double offset_angle = atan2(ri[2], ri[1]);
        if (offset_angle >= 0.)
        {
          offset_angle -= M_PI;
        }
        else if (offset_angle < 0.)
        {
          offset_angle += M_PI;
        }

        double slice_angle = 2*M_PI/(NB_APPROX + 1);
        DEBUG_PRINTF("angle: %g\n", slice_angle);
        if (ri[0] < TOL_RN)
        {
          /* Make sure that r = 0 amd not some small negative number ...*/
          ri[0] = 0.;
          ri[1] = 0.;
          ri[2] = 0.;

          /*  Pick 3 contraints as basis */
          /*  The dual variable u is very likely to be non-zero
           *  y is then also likely to be non-zero
           *  Let's make our life easier dans choose as
           *  outer approximation*/
          double offset_angle = -atan2(ui[2], ui[1]);
          unsigned pieces = (NB_APPROX/3);
          double k1 = pieces*slice_angle;
          double k2 = 2*k1;

          /* The transpose of the square submatrix of constraints, in col-major */
          const double mat[9] = {mu, cos(offset_angle),      sin(offset_angle),
                           mu, cos(offset_angle + k1), sin(offset_angle + k1),
                           mu, cos(offset_angle + k2), sin(offset_angle + k2)};
          double b[3];

/*          b[0] = problem->q[i3];
          b[1] = problem->q[i3+1];
          b[2] = problem->q[i3+2];*/

          /* now the dual var: u = A^T \lambda_r */
          b[0] = ui[0];
          b[1] = ui[1];
          b[2] = ui[2];
          solve_3x3_gepp(mat, b);

          assert(isfinite(b[0]));
          assert(isfinite(b[1]));
          assert(isfinite(b[2]));
          /* XXX Project on the dual, review if this is correct and does not
           * hurt  */
          if (b[0] < 0.) { b[0] = 0.; }
          if (b[1] < 0.) { b[1] = 0.; }
          if (b[2] < 0.) { b[2] = 0.; }
          lambda_r[offset_row] = b[0];
          lambda_r[offset_row+pieces] = b[1];
          lambda_r[offset_row+2*pieces] = b[2];
          /* now we write the estimate for y here = (mu \|u_t\|, -mu^2*ut)^T */
          /* Before that, we save ut for later in b[1:3] */
          b[1] = ui[1];
          b[2] = ui[2];

          ui[0] = mu*sqrt(ui[1]*ui[1] + ui[2]*ui[2]);
          ui[1] = -mu*mu*ui[1];
          ui[2] = -mu*mu*ui[2];

          /* now the dual var: A^T \lambda_y = (y_n, u_t) */
          b[0] = ui[0];
          /*  XXX mat is not changed, so we should not have any trouble here */
          solve_3x3_gepp(mat, b);
          assert(isfinite(b[0]));
          assert(isfinite(b[1]));
          assert(isfinite(b[2]));
          /* XXX Project on the dual, review if this is correct and does not
           * hurt  */
          if (b[0] < 0.) { b[0] = 0.; }
          if (b[1] < 0.) { b[1] = 0.; }
          if (b[2] < 0.) { b[2] = 0.; }
          lambda_y[offset_row] = b[0];
          lambda_y[offset_row+pieces] = b[1];
          lambda_y[offset_row+2*pieces] = b[2];

          DEBUG_PRINTF("new y[i] = %g %g %g\n", ui[0], ui[1], ui[2]);
        }
        else if ((ri[1]*ri[1] + ri[2]*ri[2]) < (1.+1e-10)*mu*mu * ri[0]*ri[0]) /* We should have r \in int K and u = y = 0 = dual(y) */
        {
          /* TODO? update r based on the changes in the contact forces  */
          ui[0] = 0.;
          ui[1] = 0.;
          ui[2] = 0.;
          /* lambda_r and lambda_y is already set to 0 */
        }
        else
        {
          printf("bad_angle case\n");
          ri[0] = 0.;
          ri[1] = 0.;
          ri[2] = 0.;

          ui[0] = 0.;
          ui[1] = 0.;
          ui[2] = 0.;
        }

        for (unsigned k = 0; k < NB_APPROX; ++k)
        {
          cs_entry(Ak_triplet, k + offset_row, i3, mu);
          cs_entry(Ak_triplet, k + offset_row, i3 + 1, cos(offset_angle));
          cs_entry(Ak_triplet, k + offset_row, i3 + 2, sin(offset_angle));
          offset_angle += slice_angle;
        }
        offset_row += NB_APPROX;

      }

      /* Update the dimension of Ak */
      Akmat.size0 = offset_row;
    }
    double* xtmp2 = (double*)calloc(Akmat.size0, sizeof(double));
    DEBUG_PRINT_VEC(xtmp, size);
    cs_gaxpy(NM_csc(&Akmat), xtmp, xtmp2);
    DEBUG_PRINT_VEC(xtmp2, Akmat.size0);
    for (unsigned ii = 0; ii < Akmat.size0; ++ii)
    {
      if (xtmp2[ii] < -DBL_EPSILON) printf("constraint violation at row %d, value %g\n", ii, xtmp2[ii]);
    }
    free(xtmp);
    free(xtmp2);

    cblas_dcopy(size, reaction, 1, reaction_old, 1);
    cblas_dcopy(size, velocity, 1, velocity_old, 1);

    total_residual = sqrt(total_residual);
//    optSetDblStr(solverOptPtr, "expand_delta", fmax(1e-13, fmin(1e-10, total_residual*1e-7)));
    optWriteParameterFile(solverOptPtr, msg);
//    optSetDblStr(solverOptPtr, "convergence_tolerance", 1e-12);
    printf("fc3d_AVI_gams :: residual = %g\n", total_residual);
    DEBUG_PRINTF("fc3d_AVI_gams :: residual = %g\n", total_residual);
//    done = (total_residual < options->dparam[0]);
    done = (total_residual < 1e-8);
    if (total_residual > 10*old_residual)
    {
      printf("fc3d_AVI_gams :: failure, new residual %g is bigger than old one %g\n", total_residual, old_residual);
//      goto TERMINATE;
    }
    else
    {
      old_residual = total_residual;
      if (total_residual< .9*old_residual)
      {
//        current_nb_approx = NB_APPROX;
      }
      else
      {
//        current_nb_approx += 5;
      }
    }

    if (!strcmp(solverName, "pathvi"))
    {
      optSetStrStr(solverOptPtr, "avi_start", "regular");
    }

    SN_logh5_vec_double(problem->numberOfContacts, delta_angles, "delta_angles", logger_s->group);
    SN_logh5_vec_double(problem->numberOfContacts, real_angles, "real_angles", logger_s->group);
    SN_logh5_vec_double(problem->numberOfContacts, residual_contact, "residual_contact", logger_s->group);
    /*  XXX buggy here, implement a SN_logh5_vec_integer */
    SN_logh5_vec_uint64(problem->numberOfContacts, type_contact, "type_contact", logger_s->group);
    SN_logh5_vec_uint64(problem->numberOfContacts, type_contact_avi, "type_contact_avi", logger_s->group);
    SN_logh5_end_iter(logger_s);
  }

  /********************************************************
   * Compute the residual and update the local velocities u
   ********************************************************/

  /********************************************************
   * Generate new angles
   ********************************************************/

TERMINATE:

  /* save useful data here. W should also be in the right format now  */
  SN_logh5_scalar_uinteger(iter, "iter", logger_s->file);
  SN_logh5_scalar_uinteger(done, "status", logger_s->file);
  SN_logh5_scalar_uinteger(problem->numberOfContacts, "number_contacts", logger_s->file);
  SN_logh5_vec_double(size, problem->q, "q", logger_s->file);
  SN_logh5_vec_double((unsigned)problem->numberOfContacts, problem->mu, "mu", logger_s->file);
  SN_logh5_NM(problem->M, "W", logger_s);

  /*  Truly, truly, this is the end */
  SN_logh5_end(logger_s);


  optFree(&Optr);
  optFree(&solverOptPtr);
  NM_free(&Wtmat);
  NM_free(&Emat);
  NM_free(&Akmat);

  free(reaction_old);
  free(velocity_old);
  free(tmpq);
  free(lambda_r);
  free(lambda_y);
  free(predicted_angles);
  free(delta_angles); //status = fc3d_compute_error(problem, reaction, velocity, options->dparam[0], options, &(options->dparam[1]));
  free(real_angles);
  free(residual_contact);
  free(type_contact);
  free(type_contact_avi);

  if (done)
  {
    status = 0;
    options->dparam[1] = total_residual;
  }
  else
  {
    status = 1;
    options->dparam[1] = old_residual;
  }

  options->iparam[1] = iter;
  return status;
}

void fc3d_AVI_gams_path(FrictionContactProblem* problem, double* reaction, double* velocity, int *info, SolverOptions* options)
{
  *info = fc3d_AVI_gams_base(problem, reaction, velocity, options, "path");
}

void fc3d_AVI_gams_pathvi(FrictionContactProblem* problem, double* reaction, double* velocity, int *info, SolverOptions* options)
{
  *info = fc3d_AVI_gams_base(problem, reaction, velocity, options, "pathvi");
}

#else

void fc3d_AVI_gams_path(FrictionContactProblem* problem, double *reaction, double *velocity, int *info, SolverOptions* options)
{
  printf("fc3d_gams :: gams was not enabled at compile time!\n");
  exit(EXIT_FAILURE);
}

void fc3d_AVI_gams_pathvi(FrictionContactProblem* problem, double *reaction, double *velocity, int *info, SolverOptions* options)
{
  printf("fc3d_gams :: gams was not enabled at compile time!\n");
  exit(EXIT_FAILURE);
}
#endif
