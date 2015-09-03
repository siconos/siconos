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

#include "NumericsMatrix.h"
#include "FrictionContactProblem.h"
#include "FrictionContact3D_Solvers.h"
#include "FrictionContact3D_compute_error.h"
#include "projectionOnCone.h"

#ifdef HAVE_GAMS_C_API

#include "GAMSlink.h"

#include <math.h>

#include "sanitizer.h"

#define DEBUG_NOCOLOR
#define DEBUG_STDOUT
#define DEBUG_MESSAGES
#include "debug.h"

#define NB_APPROX 10

#define TOL_RN 1e-12

#define TOL2 1e-20

#define ETERMINATE 4242

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

static int FC3D_gams_inner_loop_condensed(unsigned iter, idxHandle_t Xptr, gamsxHandle_t Gptr, optHandle_t Optr, gmoHandle_t gmoPtr, char* sysdir, char* model, double* restrict reaction, double* restrict velocity, NumericsMatrix* W, double* restrict q, NumericsMatrix* Wtmat, NumericsMatrix* Emat, NumericsMatrix* Akmat)
{

  char msg[GMS_SSSIZE];
  int status;
  unsigned size = (unsigned)W->size0;
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

  char gdxFileName[GMS_SSSIZE];
  char iterStr[4];
  sprintf(iterStr, "%d", iter);
  strncpy(gdxFileName, "fc3d_avi-condensed", sizeof(gdxFileName));
  strncat(gdxFileName, iterStr, sizeof(iterStr));
  strncat(gdxFileName, ".gdx", sizeof(gdxFileName));

  idxOpenWrite(Xptr, "fc3d_avi-condensed.gdx", "Siconos/Numerics NM_to_GDX", &status);
  if (status)
    idxerrorR(status, "idxOpenWrite");
  DEBUG_PRINT("FC3D_AVI_GAMS :: fc3d_avi-condensed.gdx opened\n");

  if ((status=NM_to_GDX(Xptr, "W", "W matrix", W))) {
    printf("Model data not written\n");
    return -ETERMINATE;
  }
  DEBUG_PRINT("FC3D_AVI_GAMS :: W matrix written\n");


  cblas_dcopy(W->size0, q, 1, reaction, 1);
  for (unsigned i = 0; i < size; i += 3)
  {
    reaction[i] = 0.;
  }


  if ((status=NM_to_GDX(Xptr, "Wt", "Wt matrix", Wtmat))) {
    printf("Model data not written\n");
    return -ETERMINATE;
  }
  DEBUG_PRINT("FC3D_AVI_GAMS :: Wt matrix written\n");

  if ((status=NM_to_GDX(Xptr, "E", "E matrix", Emat))) {
    printf("Model data not written\n");
    return -ETERMINATE;
  }
  DEBUG_PRINT("FC3D_AVI_GAMS :: E matrix written\n");

  if ((status=NM_to_GDX(Xptr, "Ak", "Ak matrix", Akmat))) {
    printf("Model data not written\n");
    return -ETERMINATE;
  }
  DEBUG_PRINT("FC3D_AVI_GAMS :: Ak matrix written\n");

  if ((status=NV_to_GDX(Xptr, "q", "q vector", q, size))) {
    printf("Model data not written\n");
    return -ETERMINATE;
  }
  DEBUG_PRINT("FC3D_AVI_GAMS :: q vector written\n");

  if ((status=NV_to_GDX(Xptr, "qt", "qt vector", reaction, size))) {
    printf("Model data not written\n");
    return -ETERMINATE;
  }
  DEBUG_PRINT("FC3D_AVI_GAMS :: qt vector written\n");

  if (idxClose(Xptr))
    idxerrorR(idxGetLastError(Xptr), "idxClose");

   cp(gdxFileName, "fc3d_avi-condensed.gdx");

  if ((status=CallGams(Gptr, Optr, sysdir, model))) {
    printf("Call to GAMS failed\n");
    return -ETERMINATE;
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
    return -ETERMINATE;
  }

  memset(velocity, 0, size*sizeof(double));
  if ((status=GDX_to_NV(Xptr, "velocity", velocity, size))) {
    printf("Model data not read\n");
    return -ETERMINATE;
  }

  double infos[] = {0., 0.};
  if ((status=GDX_to_NV(Xptr, "infos", infos, 2))) {
    printf("Model data not read\n");
    return -ETERMINATE;
  }

  if (idxClose(Xptr))
    idxerrorR(idxGetLastError(Xptr), "idxClose");

  printf("SolveStat = %d, ModelStat = %d\n", (int)infos[1], (int)infos[0]);
  gmoGetModelStatusTxt(gmoPtr, (int)infos[0], msg);
  printf(msg);
  gmoGetSolveStatusTxt(gmoPtr, (int)infos[1], msg);
  printf(msg);

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
  gmoHandle_t gmoPtr = NULL;

  int status;
  char sysdir[GMS_SSSIZE], model[GMS_SSSIZE], msg[GMS_SSSIZE];
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
  SN_Gams_set_dirs(options->solverParameters, defModel, defGAMSdir, model, sysdir, "/fc_vi-condensed.gms");


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
  strncat(msg, ".opt", sizeof(msg));

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
  else // only for path
  {
//    optSetDblStr(solverOptPtr, "convergence_tolerance", options->dparam[0]);

  }
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

  bool done = false;
  double total_residual = 0.;
  double old_residual = 1e20;

  double current_nb_approx = NB_APPROX;

  unsigned iter = 0;
  unsigned maxiter = 20;
  while (!done && (iter < maxiter))
  {
    iter++;
    total_residual = 0.;
    int solverStat = FC3D_gams_inner_loop_condensed(iter, Xptr, Gptr, Optr, gmoPtr, sysdir, model, reaction, velocity, problem->M, problem->q, &Wtmat, &Emat, &Akmat);
    DEBUG_PRINT_VEC(reaction, size);
    DEBUG_PRINT_VEC(velocity, size);

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
      DEBUG_PRINTF("contact %d, before projection: theta_r = %g\n", i, rad2deg(atan2(reaction[i3+2],reaction[i3+1])));
      projectionOnCone(&reaction[i3], mu);
      DEBUG_PRINTF("contact %d, after projection: theta_r = %g\n", i, rad2deg(atan2(reaction[i3+2],reaction[i3+1])));
    }
    /************************************************
     * (Re)compute the local velocities
     ************************************************/
    cblas_dcopy(size, problem->q, 1, velocity, 1);
    NM_gemv(1., problem->M, reaction, 1., velocity);

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
    for (unsigned i3 = 0, i = 0; i3 < size; ++i, i3 += 3)
    {
      double res = 0.;
      /* Step 2. recompute the local velocities and  */
      double* ri = &reaction[i3];
      double* ui = &velocity[i3];
      assert(i < (unsigned)problem->numberOfContacts);
      double mu = problem->mu[i];
      FrictionContact3D_unitary_compute_and_add_error(ri, ui, mu, &res);
      total_residual += res;
      unsigned p = 10;
      /* TODO we may want to revisit this, since err < TOL2 should be enough to
       * really reduce the number of constraints ...*/
      /* Well we do not want to mess with the sliding case ( both r and u on
       * the boundaries)*/
      if ((res < TOL2) && ((ri[0] < TOL_RN) || ((ri[1]*ri[1] + ri[2]*ri[2]) < mu*mu * ri[0]*ri[0])))
      {
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
      }
      else if (ri[0] > TOL_RN) // if r= 0 :(
      {
        //double delta_angle = atan2(-ri[1]*ui[2] + ui[1]*ri[2], ri[1]*ri[2] + ui[1]*ui[2]);
        double minus_r_angle = atan2(ri[2], ri[1])+ M_PI;
        double delta_angle = atan2(ui[2], ui[1]) - minus_r_angle;
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
        DEBUG_PRINTF("contact %d, delta angle = %g, theta_r = %g, theta_u = %g\n", i, rad2deg(delta_angle), rad2deg(atan2(ri[2],ri[1])), rad2deg(atan2(ui[2], ui[1])));
        /* now compute minus the angle, since we want to compute the constraints  */
        double slice_angle = delta_angle/(p-1);
        DEBUG_PRINTF("contact %d, slice_angle = %g\n", i, rad2deg(slice_angle));

        /* Generate the supporting hyperplanes  */
        double angle = minus_r_angle;
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
        DEBUG_PRINTF("contact %d, row entry %d, last_entry: angle = %g, norm = %g; coeff = %g, %g, %g\n", i, p + offset_row, rad2deg(atan2(middle_point[1], middle_point[0])), hypot(middle_point[0], middle_point[1]), -mu*hypot(middle_point[0], middle_point[1]), -middle_point[0], -middle_point[1]);
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
        display(&Akmat);
        DEBUG_PRINTF("contact %d, checking feasibility\n", i);
        free(xtmp2);
        free(pts);*/


      }
      else // r = 0, or r in int(cone) but other interactions moved u
bad_angle:
      {
        double angle = 2*M_PI/(NB_APPROX + 1);
        DEBUG_PRINTF("angle: %g\n", angle);

        for (unsigned k = 0; k < NB_APPROX; ++k)
        {
          cs_entry(Ak_triplet, k + offset_row, i3, mu);
          cs_entry(Ak_triplet, k + offset_row, i3 + 1, cos(k*angle));
          cs_entry(Ak_triplet, k + offset_row, i3 + 2, sin(k*angle));
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
    total_residual = sqrt(total_residual);
    DEBUG_PRINTF("FrictionContact3D_AVI_gams :: residual = %g\n", total_residual);
//    done = (total_residual < options->dparam[0]);
    done = (total_residual < 1e-8);
    if (total_residual > 10*old_residual)
    {
      printf("FrictionContact3D_AVI_gams :: failure, new residual %g is bigger than old one %g\n", total_residual, old_residual);
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
  }

  /********************************************************
   * Compute the residual and update the local velocities u
   ********************************************************/

  /********************************************************
   * Generate new angles
   ********************************************************/

TERMINATE:
  optFree(&Optr);
  optFree(&solverOptPtr);
  freeNumericsMatrix(&Wtmat);
  freeNumericsMatrix(&Emat);
  freeNumericsMatrix(&Akmat);

  //status = FrictionContact3D_compute_error(problem, reaction, velocity, options->dparam[0], options, &(options->dparam[1]));
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
