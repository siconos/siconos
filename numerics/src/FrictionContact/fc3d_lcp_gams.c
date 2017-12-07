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
#include "FrictionContactProblem.h"
#include "fc3d_Solvers.h"
#include "fc3d_compute_error.h"
#include "projectionOnCone.h"

#if 0
//#ifdef HAVE_GAMS_C_API

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

#define MIN_DELTA_ANGLE 1e-12

#define TOL_REFI 1e-12
#define NM_ITER_REFI 10
#define WITH_ITER_REFI

enum { TAKEOFF_CASE, STICKING_CASE, SLIDING_CASE };

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

static double solve_iterative_refinement3x3(double* restrict A, double* restrict Ainv, double* restrict x, double* restrict b, double* restrict coeffs, size_t max_iter, double tol)
{
  assert(max_iter > 0);
  double res[3] = {b[0], b[1], b[2]};
  double del[3] = {0., 0., 0.};
  double res_l1 = INFINITY;

  size_t iter = 0;

  while(iter < max_iter)
  {
    ++iter;

    /* Solve Ax = b  */
    mv3x3(Ainv, b, x);

    /* Compute r = b - Ax  */
    if (coeffs)
    {
      /*Scale x because Ainv is not A^{-1}, but A^{-1} D */
      diag_scal3(coeffs, x);
    }

    mvm3x3(A, x, res);

    res_l1 = (fabs(res[0]) + fabs(res[1]) + fabs(res[2]))/3.;

    if (res_l1 < tol)
    {
      return res_l1;
    }

    /* Solve Ad = res  */
    mv3x3(Ainv, res, del);

    /* Incremental update x += delta  */
    add3(del, x);

    /* reset residual value*/
    res[0] = b[0];
    res[1] = b[1];
    res[2] = b[2];
  }

  return res_l1;
}

static double solve_iterative_refinement3x3_t(double* restrict A, double* restrict Ainv, double* restrict x, double* restrict b, double* restrict coeffs, size_t max_iter, double tol)
{
  assert(max_iter > 0);
  double res[3] = {b[0], b[1], b[2]};
  double del[3] = {0., 0., 0.};
  double res_l1 = INFINITY;

  size_t iter = 0;

  while(iter < max_iter)
  {
    ++iter;

    /* Solve Ax = b  */
    mtv3x3(Ainv, b, x);
    /* Compute r = b - Ax  */

    if (coeffs)
    {
      /*Scale x because Ainv is not A^{-1}, but A^{-1}D */
      diag_scal3(coeffs, x);
    }

    mtvm3x3(A, x, res);

    res_l1 = (fabs(res[0]) + fabs(res[1]) + fabs(res[2]))/3.;

    if (res_l1 < tol)
    {
      return res_l1;
    }

    /* Solve Ad = res  */
    mv3x3(Ainv, res, del);
    /* Incremental update x += delta  */
    add3(del, x);

    /* reset residual value*/
    res[0] = b[0];
    res[1] = b[1];
    res[2] = b[2];
  }

  return res_l1;
}

static int FC3D_gams_inner_loop_condensed(unsigned iter, idxHandle_t Xptr, gamsxHandle_t Gptr, optHandle_t Optr, gmoHandle_t gmoPtr, char* sysdir, char* model, const char* base_name, double* restrict slack_r, double* restrict slack_y, double* restrict tmpq, double* restrict lambda_r, double* restrict lambda_y, NumericsMatrix* tildeW, double* restrict tilde_omega, double* restrict tilde_omegat, NumericsMatrix* tildeWt, NumericsMatrix* Emat, NumericsMatrix* Akmat)
{

  char msg[GMS_SSSIZE];
  int status;
  unsigned size = (unsigned)tildeW->size0;
  double infos[] = {0., 0.};
  /* Create objects */
  DEBUG_PRINT("FC3D_LCP_GAMS :: creating gamsx object\n");
  if (! gamsxCreateD (&Gptr, sysdir, msg, sizeof(msg))) {
    printf("Could not create gamsx object: %s\n", msg);
    return 1;
  }

  DEBUG_PRINT("FC3D_LCP_GAMS :: creating gdx object\n");
  if (! idxCreateD (&Xptr, sysdir, msg, sizeof(msg))) {
    printf("Could not create gdx object: %s\n", msg);
    return 1;
  }

  DEBUG_PRINT("FC3D_LCP_GAMS :: creating gmo object\n");
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
  DEBUG_PRINT("FC3D_LCP_GAMS :: fc_lcp-condensed.gdx opened\n");

  if ((status=NM_to_GDX(Xptr, "W", "W matrix", tildeW))) {
    printf("Model data not written\n");
    infos[1] = (double)-ETERMINATE;
    goto fail;
  }
  DEBUG_PRINT("FC3D_LCP_GAMS :: W matrix written\n");


  if ((status=NM_to_GDX(Xptr, "Wt", "Wt matrix", tildeWt))) {
    printf("Model data not written\n");
    infos[1] = (double)-ETERMINATE;
    goto fail;
  }
  DEBUG_PRINT("FC3D_LCP_GAMS :: Wt matrix written\n");

  if ((status=NM_to_GDX(Xptr, "E", "E matrix", Emat))) {
    printf("Model data not written\n");
    infos[1] = (double)-ETERMINATE;
    goto fail;
  }
  DEBUG_PRINT("FC3D_LCP_GAMS :: E matrix written\n");

  if ((status=NM_to_GDX(Xptr, "Ak", "Ak matrix", Akmat))) {
    printf("Model data not written\n");
    infos[1] = (double)-ETERMINATE;
    goto fail;
  }
  DEBUG_PRINT("FC3D_LCP_GAMS :: Ak matrix written\n");

  if ((status=NV_to_GDX(Xptr, "q", "q vector", tilde_omega, size))) {
    printf("Model data not written\n");
    infos[1] = (double)-ETERMINATE;
    goto fail;
  }
  DEBUG_PRINT("FC3D_LCP_GAMS :: q vector written\n");

  if ((status=NV_to_GDX(Xptr, "qt", "qt vector", tilde_omegat, size))) {
    printf("Model data not written\n");
    infos[1] = (double)-ETERMINATE;
    goto fail;
  }
  DEBUG_PRINT("FC3D_LCP_GAMS :: qt vector written\n");

/*  if ((status=NV_to_GDX(Xptr, "guess_r", "guess for r", reaction, size))) {
    printf("Model data not written\n");
    infos[1] = (double)-ETERMINATE;
    goto fail;
  }
  DEBUG_PRINT("FC3D_LCP_GAMS :: guess_r vector written\n");

  if ((status=NV_to_GDX(Xptr, "guess_y", "guess for y", velocity, size))) {
    printf("Model data not written\n");
    infos[1] = (double)-ETERMINATE;
    goto fail;
  }
  DEBUG_PRINT("FC3D_LCP_GAMS :: guess_y vector written\n");

  if ((status=NV_to_GDX(Xptr, "guess_lambda_r", "guess for lambda_r", lambda_r, Akmat->size0))) {
    printf("Model data not written\n");
    infos[1] = (double)-ETERMINATE;
    goto fail;
  }
  DEBUG_PRINT("FC3D_LCP_GAMS :: lambda_r vector written\n");

  if ((status=NV_to_GDX(Xptr, "guess_lambda_y", "guess for lambda_y", lambda_y, Akmat->size0))) {
    printf("Model data not written\n");
    infos[1] = (double)-ETERMINATE;
    goto fail;
  }
  DEBUG_PRINT("FC3D_LCP_GAMS :: lambda_y vector written\n");
*/
  if (idxClose(Xptr))
    idxerrorR(idxGetLastError(Xptr), "idxClose");

   cp(gdxFileName, "fc3d_lcp-condensed.gdx");

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
  memset(slack_r, 0, size*sizeof(double));
  if ((status=GDX_to_NV(Xptr, "sr", slack_r, size))) {
    printf("Model data not read\n");
    infos[1] = (double)-ETERMINATE;
    goto fail;
  }

  memset(slack_y, 0, size*sizeof(double));
  if ((status=GDX_to_NV(Xptr, "sy", slack_y, size))) {
    printf("Model data not read\n");
    infos[1] = (double)-ETERMINATE;
    goto fail;
  }

  if ((status=GDX_to_NV(Xptr, "infos", infos, 2))) {
    printf("Model data not read\n");
    infos[1] = (double)-ETERMINATE;
    goto fail;
  }

  if (idxClose(Xptr))
    idxerrorR(idxGetLastError(Xptr), "idxClose");

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

/* 
  size_t nb_slice;
  if (fabs(delta_angle) < 2*M_PI)
  {
    nb_slice = nb_approx - 1;
  }
  else
  {
    nb_slice = nb_approx;
  }
  double slice_angle = delta_angle/(nb_slice);
*/

static size_t fc3d_lcp_data_generation_one_reaction_force(double mu, size_t nb_angles, size_t nb_extra, size_t i3, size_t offset_row, unsigned char* indx_basis, double* restrict angles, double* restrict extra_hyperplanes, double* restrict omega_i, double* restrict tilde_omega_i, double* restrict tilde_omegat_i, double* restrict inv_change_basis_ii, CSparseMatrix* restrict Ak_triplet, CSparseMatrix* restrict E_triplet)
{
  assert(nb_angles >= 3);
  assert(omega_i);
  assert(tilde_omega_i);
  assert(tilde_omegat_i);
  assert(inv_change_basis_ii);

  double theta1 = angles[indx_basis[0]];
  double theta2 = angles[indx_basis[1]];
  double theta3 = angles[indx_basis[2]];

  int exp_factor;
  bool scaling;

  if (indx_basis[2] <= nb_angles-1)
  {
    double sin12 = sin(theta1 - theta2);
    double sin23 = sin(theta2 - theta3);
    double sin13 = sin(theta1 - theta3);

    double cos1;
    double cos2;
    double cos3;

    double sin1;
    double sin2;
    double sin3;

    /*  XXX be careful with the value of the determinant! */
    double det = sin12 + sin23 - sin13;
    if (det > 0.)
    {
      cos1 = cos(theta1);
      cos2 = cos(theta2);
      cos3 = cos(theta3);
      sin1 = sin(theta1);
      sin2 = sin(theta2);
      sin3 = sin(theta3);
    }
    else
    {
      cos1 = -cos(theta1);
      cos2 = -cos(theta2);
      cos3 = -cos(theta3);
      sin1 = -sin(theta1);
      sin2 = -sin(theta2);
      sin3 = -sin(theta3);

      /* Fastest way to flip signs  */
      sin12 = -sin12;
      sin23 = -sin23;
      sin13 = -sin13;
    }
    /* This is A_B^{-1} in col-major  */
    /* get some factor in here */

    inv_change_basis_ii[0] = sin23/mu;
    inv_change_basis_ii[1] = -sin2 + sin3;
    inv_change_basis_ii[2] = cos2 - cos3;

    inv_change_basis_ii[3] = -sin13/mu;
    inv_change_basis_ii[4] = sin1 - sin3;
    inv_change_basis_ii[5] = cos3 - cos1;

    inv_change_basis_ii[6] = sin12/mu;
    inv_change_basis_ii[7] = sin2 - sin1;
    inv_change_basis_ii[8] = cos1 - cos2;

    frexp(theta1 - theta2, &exp_factor);
    exp_factor = exp_factor + exp_factor/3;
    printf("DEBUG SCALING exp = %d\n", exp_factor);
    for (size_t j = 0; j < 9; ++j)
    {
      inv_change_basis_ii[j] = ldexp(inv_change_basis_ii[j], -exp_factor);
    }

    scaling = true;
  }
  else
  {
    double cos13 = -cos(theta1 - theta3);
    double cos12 = -cos(theta1 - theta2);
    double cos2 = cos(theta2);
    double cos3 = cos(theta3);
    double sin2 = sin(theta2);
    double sin3 = sin(theta3);
    if (sin(theta2 - theta3) < 0)
    {
      cos13 = -cos13;
      cos12 = -cos12;
      cos2  = -cos2;
      cos3  = -cos3;
      sin2  = -sin2;
      sin3  = -sin3;
    }
    inv_change_basis_ii[0] = 1.;
    inv_change_basis_ii[1] = 0.;
    inv_change_basis_ii[2] = 0.;

    inv_change_basis_ii[3] = cos13/mu;
    inv_change_basis_ii[4] = cos3;
    inv_change_basis_ii[5] = sin3;

    inv_change_basis_ii[6] = cos12/mu;
    inv_change_basis_ii[7] = cos2;
    inv_change_basis_ii[8] = sin2;

    scaling = false;
    exp_factor = 0;
  }

  /* E_ii = AbT [* * *]
   *            [0 0 0]
   *            [0 0 0]
   */
  for (size_t j = 0; j < 3; ++j)
  {
    double factor = scaling ? inv_change_basis_ii[3*j] : 1.;
    printf("DEBUG EMAT factor = %e\n", factor);
    cs_entry(E_triplet, i3,     i3 + j, factor*inv_change_basis_ii[0]);
    cs_entry(E_triplet, i3 + 1, i3 + j, factor*inv_change_basis_ii[3]);
    cs_entry(E_triplet, i3 + 2, i3 + j, factor*inv_change_basis_ii[6]);
  }

  /*  Compute the transformed \tilde{ω} = A_B⁻ᵀ ω  */
  mtv3x3(inv_change_basis_ii, omega_i, tilde_omega_i);
  double omegat_i[3] = {0., omega_i[1], omega_i[2]};
  mtv3x3(inv_change_basis_ii, omegat_i, tilde_omegat_i);

  /*  Compute A_{\bar{B}} A^{-1}_B */
  unsigned j = 0;
  for (unsigned i = 0; i < nb_angles; ++i)
  {
    if (i == indx_basis[j])
    {
#ifndef NDBEUG
      double angle = angles[i];
      double p[3] = {mu, cos(angle), sin(angle)};
      p[0] = ldexp(p[0], exp_factor);
      p[1] = ldexp(p[1], exp_factor);
      p[2] = ldexp(p[2], exp_factor);
      double r[3];
      mtv3x3(inv_change_basis_ii, p, r);
      printf("DEBUG INV MAT: r = [%e; %e; %e]; j = %d\n", r[0], r[1], r[2], j);
      for (unsigned k = 0; k < 3; k++)
      {
        if (k != j)
        {
          if (fabs(r[k]) > 1e-10) { printf("r[%d] = %e but should be 0 (j = %d)\n", k, r[k], j); }
          //          assert(fabs(r[k]) < 1e-10);
        }
      }
#endif /*  NDEBUG */
      if (j < 2) { ++j; }
    }
    else
    {
      double angle = angles[i];
      double p[3] = {mu, cos(angle), sin(angle)};
      double r[3];
      p[0] = ldexp(p[0], -exp_factor);
      p[1] = ldexp(p[1], -exp_factor);
      p[2] = ldexp(p[2], -exp_factor);
      mtv3x3(inv_change_basis_ii, p, r);
      cs_entry(Ak_triplet, offset_row, i3,     r[0]);
      cs_entry(Ak_triplet, offset_row, i3 + 1, r[1]);
      cs_entry(Ak_triplet, offset_row, i3 + 2, r[2]);
      ++offset_row;
    }
  }

  /*  Otherwise we have a row of 0 ... */
  for (unsigned i = 0; i < nb_extra; ++i)
  {
    int flip = 1;
    if (i + nb_angles == indx_basis[j])
    {
#ifndef NDBEUG
      flip = -flip;
      double angle = angles[nb_angles + i];
      double p[3] = {0, flip*sin(angle), -flip*cos(angle)};
      p[0] = ldexp(p[0], exp_factor);
      p[1] = ldexp(p[1], exp_factor);
      p[2] = ldexp(p[2], exp_factor);
      double r[3];
      mtv3x3(inv_change_basis_ii, p, r);
      printf("DEBUG INV MAT: r = [%e; %e; %e]; j = %d\n", r[0], r[1], r[2], j);
      for (unsigned k = 0; k < 3; k++)
      {
        if (k != j)
        {
          if (fabs(r[k]) > 1e-10) { printf("r[%d] = %e but should be 0 (j = %d)\n", k, r[k], j); }
          //          assert(fabs(r[k]) < 1e-10);
        }
      }
#endif /*  NDEBUG */
      if (j < 2) { ++j; }
    }
    else
    {
      double r[3];
      double p[3];
      int exp_factor_extra;
      if (scaling)
      {
        frexp(angles[0] - angles[nb_angles-1], &exp_factor_extra);
        exp_factor_extra -= 1;
      }
      else
      {
        exp_factor_extra = 0;
      }
      p[0] = ldexp(extra_hyperplanes[3*i], -exp_factor_extra);
      p[1] = ldexp(extra_hyperplanes[3*i+1], -exp_factor_extra);
      p[2] = ldexp(extra_hyperplanes[3*i+2], -exp_factor_extra);
      /*    p[0] = extra_hyperplanes[3*i];
            p[1] = extra_hyperplanes[3*i+1];
            p[2] = extra_hyperplanes[3*i+2];*/
      mtv3x3(inv_change_basis_ii, p, r);
      cs_entry(Ak_triplet, offset_row, i3,     r[0]);
      cs_entry(Ak_triplet, offset_row, i3 + 1, r[1]);
      cs_entry(Ak_triplet, offset_row, i3 + 2, r[2]);
      ++offset_row;
    }
  }

  return offset_row;
}


static void FC3D_gams_generate_first_constraints(NumericsMatrix* Akmat, NumericsMatrix* Emat, double* restrict mus, double* restrict omega, double* restrict tilde_omega, double* restrict tilde_omegat, double* restrict changeBasis)
{
  unsigned nb_contacts = (unsigned)Akmat->size1/3;
  assert(nb_contacts*3 == (unsigned)Akmat->size1);
  unsigned nb_approx = (unsigned)Akmat->size0/nb_contacts;
  assert(nb_approx*nb_contacts == (unsigned)Akmat->size0);
  unsigned offset_row = 0;

  double slice_angle = 2*M_PI/(NB_APPROX);
  DEBUG_PRINTF("angle: %g\n", slice_angle);

  double angles[NB_APPROX-1];
  for (size_t i = 0; i < NB_APPROX-1; ++i)
  {
    angles[i] = i*slice_angle;
  }

  /*  Select the first 3 constraints as our basis var */
  unsigned char indx_basis[3] = {0, 1, 2};

  CSparseMatrix* Ak_triplet = NM_triplet(Akmat);
  CSparseMatrix* E_triplet = NM_triplet(Emat);

  for (unsigned j = 0, i3 = 0, indxMat = 0; j < nb_contacts; ++j, i3 += 3, indxMat += 9)
  {
    /* TODO add a real starting angle  */
    double *inv_change_basis_ii = &changeBasis[indxMat];
    double *tilde_omegat_i = &tilde_omegat[i3];
    double *tilde_omega_i = &tilde_omega[i3];
    double *omega_i = &omega[i3];
    offset_row = fc3d_lcp_data_generation_one_reaction_force(mus[j], NB_APPROX-1, 0, i3, offset_row, indx_basis, angles, NULL, omega_i, tilde_omega_i, tilde_omegat_i, inv_change_basis_ii, Ak_triplet, E_triplet);
  }
}

static int fc3d_lcp_gams_base(FrictionContactProblem* problem, double *reaction, double *velocity, SolverOptions* options, const char* solverName)
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

  int status;
  char sysdir[GMS_SSSIZE], model[GMS_SSSIZE], msg[GMS_SSSIZE], template_filename[GMS_SSSIZE], hdf5_filename[GMS_SSSIZE];
  const char defModel[] = SPACE_CONC(GAMS_MODELS_SHARE_DIR, "/fc_lcp-condensed.gms");
  const char defGAMSdir[] = GAMS_DIR;

  unsigned size = (unsigned) problem->dimension*problem->numberOfContacts;

  NumericsMatrix Wtmat;
  NumericsMatrix Emat;
  NumericsMatrix Akmat;

  NM_null(&Wtmat);
  NM_null(&Emat);
  NM_null(&Akmat);

  DEBUG_PRINT("FC3D_LCP_GAMS :: seeting basic directories\n");
  SN_Gams_set_dirs(options->solverParameters, defModel, defGAMSdir, model, sysdir, "/fc_lcp-condensed.gms");

  const char* filename = GAMSP_get_filename(options->solverParameters);

  DEBUG_PRINT("FC3D_LCP_GAMS :: creating opt object\n");
  if (! optCreateD (&Optr, sysdir, msg, sizeof(msg))) {
    printf("Could not create opt object: %s\n", msg);
    return 1;
  }

  DEBUG_PRINT("FC3D_LCP_GAMS :: creating solveropt object\n");
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
  if (strcmp(solverName, "path"))
  {
    optSetStrStr(Optr, "emp", solverName);
    optSetStrStr(solverOptPtr, "avi_start", "ray_first");
    optSetStrStr(solverOptPtr, "ratio_tester", "expand");
  }
  else // only for path
  {
//    optSetDblStr(solverOptPtr, "convergence_tolerance", options->dparam[0]);
//    optSetIntStr(solverOptPtr, "output_linear_model", 1);
    optSetDblStr(solverOptPtr, "proximal_perturbation", 0.);
    optSetStrStr(solverOptPtr, "crash_method", "none");
    optSetIntStr(solverOptPtr, "crash_perturb", 0);
//    optSetIntStr(solverOptPtr, "output_minor_iterations_frequency", 1);
//    optSetIntStr(solverOptPtr, "output_linear_model", 1);

  }
  optSetIntStr(solverOptPtr, "minor_iteration_limit", 100000);
  optSetDblStr(solverOptPtr, "expand_delta", 1e-10);
//  optSetDblStr(solverOptPtr, "convergence_tolerance", 1e-12);
  optSetDblStr(solverOptPtr, "convergence_tolerance", options->dparam[0]);
  optWriteParameterFile(solverOptPtr, msg);

  optSetIntStr(Optr, "Keep", 0);


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

  unsigned size_l = (NB_APPROX+3)*problem->numberOfContacts;
  /* Create a csc matric for Ab  */
  NumericsMatrix Ab;
  NM_null(&Ab);
  Ab.storageType = NM_SPARSE;
  Ab.size0 = 3*problem->numberOfContacts;
  Ab.size1 = 3*problem->numberOfContacts;
  NM_csc_alloc(&Ab, 3*Ab.size0);

  /* Let the dark magic begin: we create the indices and colptr  */
  CS_INT* Ab_rowindx = NM_sparse(&Ab)->csc->i;
  CS_INT* Ab_colptr = NM_sparse(&Ab)->csc->p;

  for (unsigned i = 0, i3 = 0, data_indx = 0; i3 < 3*problem->numberOfContacts; i3 += 3, i += 9)
  {
    Ab_rowindx[i] = Ab_rowindx[i+3] = Ab_rowindx[i+6] = i3;
    Ab_rowindx[i+1] = Ab_rowindx[i+4] = Ab_rowindx[i+7] = i3+1;
    Ab_rowindx[i+2] = Ab_rowindx[i+5] = Ab_rowindx[i+8] = i3+2;

    Ab_colptr[i3] = data_indx;
    data_indx += 3;
    Ab_colptr[i3+1] = data_indx;
    data_indx += 3;
    Ab_colptr[i3+2] = data_indx;
    data_indx += 3;
  }

  Ab_colptr[Ab.size1] = 3*Ab.size0;

  /* csc matrix for \tilde{W} */
  NumericsMatrix tildeW;
  NM_null(&tildeW);
  tildeW.storageType = NM_SPARSE;

  NumericsMatrix tildeWt;
  NM_null(&tildeWt);
  tildeWt.storageType = NM_SPARSE;

  double* tilde_omega = (double*) malloc(3*problem->numberOfContacts*sizeof(double));
  double* tilde_omegat = (double*) malloc(3*problem->numberOfContacts*sizeof(double));
  /* Fills Akmat and Ab.x  */
  FC3D_gams_generate_first_constraints(&Akmat, &Emat, problem->mu, problem->q, tilde_omega, tilde_omegat, NM_csc(&Ab)->x);

  NumericsMatrix* tmpWmat = NM_create(NM_SPARSE, 3*problem->numberOfContacts, 3*problem->numberOfContacts);

  double* tmpq = (double*)malloc(size * sizeof(double));
  double* lambda_r = (double*)calloc(size_l, sizeof(double));
  double* lambda_y = (double*)calloc(size_l, sizeof(double));
  double* reaction_old = (double*)calloc(size, sizeof(double));
  double* velocity_old = (double*)calloc(size, sizeof(double));

  double* slack_r = (double*)calloc(size, sizeof(double));
  double* slack_y = (double*)calloc(size, sizeof(double));

  double* coeffs = (double*)calloc(size, sizeof(double));

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

  bool done = false;
  double total_residual = 0.;
  double old_residual = 1e20;

  double current_nb_approx = NB_APPROX;

  unsigned iter = 0;
  unsigned maxiter = 20;

  /* Logger starting  */
  if (filename)
  {
    strncpy(hdf5_filename, filename, sizeof(hdf5_filename));
  }
  else
  {
    strncpy(hdf5_filename, "logger", sizeof(hdf5_filename));
  }

  strncat(hdf5_filename, ".hdf5", sizeof(hdf5_filename) - strlen(hdf5_filename) - 1);

  SN_logh5* logger_s = SN_logh5_init(hdf5_filename, maxiter);

  while (!done && (iter < maxiter))
  {
    iter++;
    total_residual = 0.;
    filename_datafiles(iter, filename, sizeof(template_filename), template_filename);

    SN_logh5_new_iter(iter, logger_s);

    SN_logh5_scalar_double(old_residual, "old_residual", logger_s->group);
    SN_logh5_vec_double(size, reaction, "reaction_guess", logger_s->group);
    SN_logh5_vec_double(size, velocity, "velocity_guess", logger_s->group);

    SN_logh5_vec_double(Akmat.size0, lambda_r, "lambda_r", logger_s->group);
    SN_logh5_vec_double(Akmat.size0, lambda_y, "lambda_y", logger_s->group);

    /* Compute tildeW and tildeWt*/
    NM_gemm(1., problem->M, &Ab, 0., tmpWmat);

    CSparseMatrix* AbT = cs_transpose(NM_sparse(&Ab)->csc, 1);
    NM_sparse(&tildeW)->csc = cs_multiply(AbT, NM_csc(tmpWmat));
    NM_update_size(&tildeW);

    /* Now compute tildeWt */
    cs_fkeep(NM_csc(tmpWmat), &SN_rm_normal_part, NULL);
    NM_sparse(&tildeWt)->csc = cs_multiply(AbT, NM_csc(tmpWmat));
    NM_update_size(&tildeWt);

#ifndef  NDEBUG
    NumericsMatrix* NM_AbT = NM_create(NM_SPARSE, AbT->m,  AbT->n);
    NM_sparse(NM_AbT)->csc = AbT;
    SN_logh5_NM(NM_AbT, "AbT", logger_s);
    NM_free(NM_AbT);
#else
    cs_spfree(AbT);
#endif

    DEBUG_PRINT("FC3D_LCP_GAMS :: tildeWt matrix constructed\n");


    int solverStat = FC3D_gams_inner_loop_condensed(iter, Xptr, Gptr, Optr, gmoPtr, sysdir, model, template_filename, slack_r, slack_y, tmpq, lambda_r, lambda_y, &tildeW, tilde_omega, tilde_omegat, &tildeWt, &Emat, &Akmat);
    SN_logh5_vec_double(size, slack_r, "reaction_slack", logger_s->group);
    SN_logh5_vec_double(size, slack_y, "velocity_slack", logger_s->group);

    double* change_basis = NM_csc(&Ab)->x;
/*     double* change_basis_inv = NM_csc(&Ab_real)->x;
    for (unsigned i3 = 0, indxMat = 0; i3 < size; i3 += 3, indxMat += 9)
    {
#ifdef WITH_ITER_REFI
      solve_iterative_refinement3x3(&change_basis[indxMat], &change_basis_inv[indxMat], &reaction[i3], &slack_r[i3], coeffs[i3], NB_ITER_REFI, TOL_REFI);
#else
      mv3x3(&change_basis[indxMat], &slack_r[i3], &reaction[i3]);
#endif
    }
    */
    //DEBUG_PRINT_VEC(reaction, size);
    //DEBUG_PRINT_VEC(velocity, size);

    // CSC form should be OK
    SN_logh5_NM(&Akmat, "Akmat", logger_s);
    SN_logh5_NM(&Ab, "Ab", logger_s);
    SN_logh5_NM(&Emat, "Emat", logger_s);
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
      double mu = problem->mu[i];
      /* Step 1. project r on the cone */

//      DEBUG_PRINTF("contact %d, before projection: theta_r = %.*e\n", i, DECIMAL_DIG, rad2deg(atan2(reaction[i3+2], reaction[i3+1])));
      projectionOnCone(&reaction[i3], mu);
//      DEBUG_PRINTF("contact %d, after projection:  theta_r = %.*e\n", i, DECIMAL_DIG, rad2deg(atan2(reaction[i3+2], reaction[i3+1])));
    }

    memset(lambda_r, 0, size_l * sizeof(double));
    memset(lambda_y, 0, size_l * sizeof(double));
    memset(tilde_omega, 0, 3*problem->numberOfContacts*sizeof(double));
    memset(tilde_omegat, 0, 3*problem->numberOfContacts*sizeof(double));

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
    NM_clearSparseStorage(&Emat);
    /* This is now a upper bound ...  */
    Akmat.matrix2->triplet = cs_spalloc(NB_APPROX*problem->numberOfContacts, size, NB_APPROX*problem->numberOfContacts*3, 1, 1);
    Emat.matrix2->triplet = cs_spalloc(3*problem->numberOfContacts, size, NB_APPROX*problem->numberOfContacts, 1, 1);
    Ak_triplet = NM_triplet(&Akmat);
    CSparseMatrix* E_triplet = NM_triplet(&Emat);
    /************************************************
     * Compute the error on each contact point + 
     ************************************************/
    unsigned offset_row = 0;
    double* xtmp = (double*)calloc(size, sizeof(double));
    for (unsigned i3 = 0, i = 0, indxMat = 0; i3 < size; ++i, i3 += 3, indxMat += 9)
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
      fc3d_unitary_compute_and_add_error(ri, ui, mu, &res);
      residual_contact[i] = sqrt(res);
      DEBUG_EXPR_WE(if (res > old_residual) { printf("Contact %d, res = %g > %g = old_residual\n", i, sqrt(res), old_residual); });
      total_residual += res;
      unsigned p = NB_APPROX;
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
      else if (ri[0] > TOL_RN) // if r= 0 :(
      {
        //double delta_angle = atan2(-ri[1]*ui[2] + ui[1]*ri[2], ri[1]*ri[2] + ui[1]*ui[2]);
        double minus_r_angle = atan2(ri[2], ri[1])+ M_PI;
        double delta_angle = atan2(ui[2], ui[1]) - minus_r_angle;
        delta_angles[i] = rad2deg(delta_angle);
        real_angles[i] = rad2deg(-minus_r_angle);
        if ((fabs(delta_angle) > M_PI/2))
        {
          if (fabs(delta_angle+2*M_PI) > M_PI/2)
          {
            printf("Contact %d, something bad happened, angle value is %g (rad) or %g (deg)\n", i, delta_angle, rad2deg(delta_angle));
            printf("r = [%g; %g; %g]\tu = [%g; %g; %g]\tres = %g\n", ri[0], ri[1], ri[2], ui[0], ui[1], ui[2], sqrt(res));
            if (((ri[1]*ri[1] + ri[2]*ri[2]) < mu*mu * ri[0]*ri[0]))
            {
              printf("r is the in the interior of the cone ... |r_r| = %g < %g = r_n*mu\n", sqrt((ri[1]*ri[1] + ri[2]*ri[2])), sqrt(mu*mu * ri[0]*ri[0]));
            }
            goto bad_angle;
          }
          else
          {
            delta_angle = (delta_angle + 2*M_PI);
          }
        }

        /* Ok so here we support that we are in the sliding case */
        type_contact[i] = SLIDING_CASE;
        DEBUG_PRINTF("contact %d, delta angle = %g, theta_r = %.*e, theta_u = %.*e\n", i, rad2deg(delta_angle), DECIMAL_DIG, rad2deg(atan2(ri[2],ri[1])), 
            DECIMAL_DIG, rad2deg(atan2(ui[2], ui[1])));
        if (fabs(delta_angle) < MIN_DELTA_ANGLE) { printf("Contact %d, delta_angle too small %g; set to 1e-12", i, delta_angle); delta_angle = copysign(MIN_DELTA_ANGLE, delta_angle);}

        /* now compute minus the angle, since we want to compute the constraints  */
        double slice_angle = delta_angle/(p-1);
        DEBUG_PRINTF("contact %d, slice_angle = %g\n", i, rad2deg(slice_angle));

        /* Generate the supporting hyperplanes  */
        double* angles = malloc((p+2) * sizeof(double));
        angles[0] = minus_r_angle;
        unsigned offset_row_bck = offset_row;
        for (unsigned j = 1; j < p; ++j)
        {
          angles[j] = angles[j-1] + slice_angle;
          DEBUG_PRINTF("contact %d, row entry %d, picking a point at the angle %g\n", i, offset_row + j, rad2deg(angles[j]));
        }
        //assert(fabs(angle - minus_r_angle - angle) < 1e-12);
        if(fabs(angles[p-1] - minus_r_angle - delta_angle) > 1e-12)
        {
          printf("warning, big difference betwwen original angle and result: %g !\n", angles[p-1] - minus_r_angle - delta_angle);
        }

        /* Add the last constraint  */
        /* XXX we already have computed those ...  --xhub */
        double middle_point[] = {(cos(minus_r_angle) + cos(angles[p-1]))/2., (sin(minus_r_angle) + sin(angles[p-1]))/2.};

        unsigned nb_closing_hyperplanes;
#ifdef SMALL_APPROX
        double closing_hyperplanes[3];
        closing_hyperplanes[0] = -mu*(hypot(middle_point[0], middle_point[1]));
        closing_hyperplanes[1] = -cos((minus_r_angle+angles[p-1])/2);
        closing_hyperplanes[2] = -sin((minus_r_angle+angles[p-1])/2);
        /* Update the row index */
        nb_closing_hyperplanes = 1;
#else

        double closing_hyperplanes[6];
        closing_hyperplanes[0] = 0.;
        closing_hyperplanes[3] = 0.;

        if (delta_angle > 0) /* We need to rotate - pi/2 the original angle */
        {
          closing_hyperplanes[1] = cos(minus_r_angle - M_PI/2); /* XXX there are formulae for this ... */
          closing_hyperplanes[2] = sin(minus_r_angle - M_PI/2);
          closing_hyperplanes[4] = cos(angles[p-1] + M_PI/2);
          closing_hyperplanes[5] = sin(angles[p-1] + M_PI/2);
          angles[p] = minus_r_angle;
          angles[p+1] = angles[p-1];
        }
        else /* We need to rotate of pi/2 */
        {
          closing_hyperplanes[1] = cos(minus_r_angle + M_PI/2);
          closing_hyperplanes[2] = sin(minus_r_angle + M_PI/2);
          closing_hyperplanes[4] = cos(angles[p-1] - M_PI/2);
          closing_hyperplanes[5] = sin(angles[p-1] - M_PI/2);
          angles[p] = angles[p-1];
          angles[p+1] = minus_r_angle;
        }

        nb_closing_hyperplanes = 2;
#endif

        unsigned char middle = (p-1)/2;
//        unsigned char indx_basis[3] = {middle-1, middle, middle+1};
        unsigned char indx_basis[3] = {middle, p, p+1};
        double *inv_change_basis_ii = &(NM_csc(&Ab)->x[indxMat]);

        offset_row = fc3d_lcp_data_generation_one_reaction_force(mu, p, nb_closing_hyperplanes, i3, offset_row, indx_basis, angles, closing_hyperplanes, &problem->q[i3], &tilde_omega[i3], &tilde_omegat[i3], inv_change_basis_ii, Ak_triplet, E_triplet);

        free(angles);
        angles = NULL;

        /* update ri =   */
        unsigned pos = indx_basis[1];
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

      }
      else // r = 0, or r in int(cone) but other interactions moved u
bad_angle:
      {
        double slice_angle = 2*M_PI/(NB_APPROX + 1);
        DEBUG_PRINTF("angle: %g\n", slice_angle);
        if (ri[0] < TOL_RN)
        {
          type_contact[i] = TAKEOFF_CASE;
        }
        else if ((ri[1]*ri[1] + ri[2]*ri[2]) < (1.+1e-10)*mu*mu * ri[0]*ri[0]) /* We should have r \in int K and u = y = 0 = dual(y) */
        {
          type_contact[i] = STICKING_CASE;
          /* TODO? update r based on the changes in the contact forces  */
          /* lambda_r and lambda_y is already set to 0 */
        }
        else
        {
          printf("bad_angle case\n");
        }

        double angles[NB_APPROX-1];
        for (size_t i = 0; i < NB_APPROX-1; ++i)
        {
          angles[i] = i*slice_angle;
        }

        /*  Select the first 3 constraints as our basis var */
        unsigned char indx_basis[3] = {0, 1, 2};

          /* TODO add a real starting angle  */
        double *inv_change_basis_ii = &(NM_csc(&Ab)->x[indxMat]);
        offset_row = fc3d_lcp_data_generation_one_reaction_force(mu, NB_APPROX-1, 0, i3, offset_row, indx_basis, angles, NULL, &problem->q[i3], &tilde_omega[i3], &tilde_omegat[i3], inv_change_basis_ii, Ak_triplet, E_triplet);

      }

    }
    /* Update the dimension of Ak */
    Akmat.size0 = offset_row;
    Akmat.size1 = size;
    Ak_triplet->m = offset_row;

    cblas_dcopy(size, reaction, 1, reaction_old, 1);
    cblas_dcopy(size, velocity, 1, velocity_old, 1);

    total_residual = sqrt(total_residual);
    optSetDblStr(solverOptPtr, "expand_delta", fmax(1e-13, fmin(1e-10, total_residual*1e-7)));
    optWriteParameterFile(solverOptPtr, msg);
//    optSetDblStr(solverOptPtr, "convergence_tolerance", 1e-12);
    printf("FrictionContact3D_LCP_gams :: residual = %g\n", total_residual);
    DEBUG_PRINTF("FrictionContact3D_LCP_gams :: residual = %g\n", total_residual);
//    done = (total_residual < options->dparam[0]);
    done = (total_residual < 1e-8);
    if (total_residual > 10*old_residual)
    {
      printf("FrictionContact3D_LCP_gams :: failure, new residual %g is bigger than old one %g\n", total_residual, old_residual);
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
//      optSetStrStr(solverOptPtr, "avi_start", "regular");
    }

    SN_logh5_vec_double(problem->numberOfContacts, delta_angles, "delta_angles", logger_s->group);
    SN_logh5_vec_double(problem->numberOfContacts, real_angles, "real_angles", logger_s->group);
    SN_logh5_vec_double(problem->numberOfContacts, residual_contact, "residual_contact", logger_s->group);
    /*  XXX buggy here, implement a SN_logh5_vec_integer */
    SN_logh5_vec_int64(problem->numberOfContacts, type_contact, "type_contact", logger_s->group);
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
  SN_logh5_NM(problem->M, "W", logger_s);

  /*  Truly, truly, this is the end */
  SN_logh5_end(logger_s);


  optFree(&Optr);
  optFree(&solverOptPtr);
  NM_free(&Wtmat);
  NM_free(&Emat);
  NM_free(&Akmat);
  NM_free(tmpWmat);

  free(reaction_old);
  free(velocity_old);
  free(tmpq);
  free(lambda_r);
  free(lambda_y);
  free(slack_r);
  free(slack_y);
  free(predicted_angles);
  free(delta_angles); //status = FrictionContact3D_compute_error(problem, reaction, velocity, options->dparam[0], options, &(options->dparam[1]));
  free(real_angles);
  free(residual_contact);
  free(type_contact);
  free(tilde_omega);
  free(tilde_omegat);
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
  return status;
}

void fc3d_lcp_gams_path(FrictionContactProblem* problem, double* reaction, double* velocity, int *info, SolverOptions* options)
{
  *info = fc3d_lcp_gams_base(problem, reaction, velocity, options, "path");
}

void fc3d_lcp_gams_pathvi(FrictionContactProblem* problem, double* reaction, double* velocity, int *info, SolverOptions* options)
{
  *info = fc3d_lcp_gams_base(problem, reaction, velocity, options, "pathvi");
}

#else

void fc3d_lcp_gams_path(FrictionContactProblem* problem, double *reaction, double *velocity, int *info, SolverOptions* options)
{
  printf("fc3d_gams :: gams was not enabled at compile time!\n");
  exit(EXIT_FAILURE);
}

void fc3d_lcp_gams_pathvi(FrictionContactProblem* problem, double *reaction, double *velocity, int *info, SolverOptions* options)
{
  printf("fc3d_gams :: gams was not enabled at compile time!\n");
  exit(EXIT_FAILURE);
}
#endif
