/* GAMS stuff */

#define _XOPEN_SOURCE 700

#if (__linux ||  __APPLE__)
#elif _MSC_VER
#define strdup _strdup
#else
static inline char* strdup(const char* src)
{
  size_t len = strlen(src) + 1;
  char* dest = (char*)malloc(len * sizeof(char));
  strncpy(dest, src, len);
  return dest;
}
#endif


#ifdef HAVE_GAMS_C_API
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>

#include "CSparseMatrix.h"
#include "NumericsMatrix.h"
#include "FrictionContactProblem.h"
#include "SolverOptions.h"

#include "GAMSlink.h"

#include <math.h>

#include "sanitizer.h"

#define DEBUG_NOCOLOR
//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "debug.h"

#define ETERMINATE 4242

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

void filename_datafiles(const int iter, const int solverId, const char* base_name, unsigned len, char* template_name, char* log_filename)
{
  char iterStr[40];
  snprintf(iterStr, sizeof(iterStr), "-i%d-%s", iter, solver_options_id_to_name(solverId));
  if (base_name)
  {
    strncpy(template_name, base_name, len);
    strncpy(log_filename, base_name, len);
  }
  else
  {
    strncpy(template_name, "fc3d_avi-condensed", len);
    strncpy(log_filename, "fc3d_avi-condense-log", len);
  }

  strncat(template_name, iterStr, len - strlen(template_name) - 1);
  strncat(log_filename, iterStr, len - strlen(log_filename) - 1);
  strncat(log_filename, ".log", len - strlen(log_filename) - 1);
}



int NM_to_GDX(idxHandle_t Xptr, const char* name, const char* descr, NumericsMatrix* M)
{
  char msg[GMS_SSSIZE];

  int dims[2];
  dims[0] = M->size0;
  dims[1] = M->size1;
  if (idxDataWriteStart(Xptr, name, descr, 2, dims, msg, GMS_SSSIZE) == 0)
    idxerrorR(idxGetLastError(Xptr), "idxDataWriteStart");

  switch (M->storageType)
  {
  case NM_DENSE:
  {
    assert(M->matrix0);
    idxDataWriteDenseColMajor(Xptr, 2, M->matrix0);
    break;
  }
  case NM_SPARSE_BLOCK: /* Perform a conversion to sparse storage */
  case NM_SPARSE:
  {
    CSparseMatrix* cs = NM_csc(M);
    assert(cs->p);
    assert(cs->i);
    assert(cs->x);
    int* p_int = (int*)malloc((cs->n+1) * sizeof(int));
    int* i_int = (int*)malloc(cs->nzmax * sizeof(int));
    assert(cs->n == M->size1);
    assert(cs->m == M->size0);
    for (unsigned i = 0; i < cs->n+1; ++i)
    {
      p_int[i] = (int) cs->p[i];
    }

    for (unsigned i = 0; i < cs->nzmax; ++i)
    {
      i_int[i] = (int) cs->i[i];
    }

    idxDataWriteSparseColMajor(Xptr, p_int, i_int, cs->x);

    free(p_int);
    free(i_int);
    break;
  }
  default:
  {
    printf("NM_to_GDX :: unsupported matrix storage\n");
    exit(EXIT_FAILURE);
  }
  }


  if (0==idxDataWriteDone(Xptr))
    idxerrorR(idxGetLastError(Xptr), "idxDataWriteDone");

  return 0;
}int SN_gams_solve(unsigned iter, optHandle_t Optr, char* sysdir, char* model, const char* base_name, SolverOptions* options, SN_GAMS_gdx* gdx_data)
{
  assert(gdx_data);
  SN_GAMS_NM_gdx* mat_for_gdx = gdx_data->mat_for_gdx;
  SN_GAMS_NV_gdx* vec_for_gdx = gdx_data->vec_for_gdx;
  SN_GAMS_NV_gdx* vec_from_gdx = gdx_data->vec_from_gdx;

  char msg[GMS_SSSIZE];
  int status;
  gamsxHandle_t Gptr = NULL;
  idxHandle_t Xptr = NULL;
  gmoHandle_t gmoPtr = NULL;
  double infos[4] = {0.};
  /* Create objects */

  DEBUG_PRINT("FC3D_AVI_GAMS :: creating gamsx object\n");
  if (! gamsxCreateD (&Gptr, sysdir, msg, sizeof(msg))) {
    fprintf(stderr, "Could not create gamsx object: %s\n", msg);
    return 1;
  }

  DEBUG_PRINT("FC3D_AVI_GAMS :: creating gdx object\n");
  if (! idxCreateD (&Xptr, sysdir, msg, sizeof(msg))) {
    fprintf(stderr, "Could not create gdx object: %s\n", msg);
    return 1;
  }

  DEBUG_PRINT("FC3D_AVI_GAMS :: creating gmo object\n");
  if (! gmoCreateD (&gmoPtr, sysdir, msg, sizeof(msg))) {
    fprintf(stderr, "Could not create gmo object: %s\n", msg);
    return 1;
  }

  /* create input and output gdx names*/
  char gdxFileName[GMS_SSSIZE];
  char solFileName[GMS_SSSIZE];
//  char paramFileName[GMS_SSSIZE];

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

   while (mat_for_gdx)
   {
     char mat_descr[30];
     assert(mat_for_gdx->name);
     assert(mat_for_gdx->mat);
     snprintf(mat_descr, sizeof(mat_descr), "%s matrix", mat_for_gdx->name);
     if ((status=NM_to_GDX(Xptr, mat_for_gdx->name, mat_descr, mat_for_gdx->mat))) {
       fprintf(stderr, "Model data for matrix %s not written\n", mat_for_gdx->name);
       infos[1] = (double)-ETERMINATE;
       goto fail;
     }
     DEBUG_PRINTF("GAMSlink :: %s matrix written\n", mat_for_gdx->name);
     mat_for_gdx = mat_for_gdx->next;
   }

   while (vec_for_gdx)
   {
     char vec_descr[30];
     assert(vec_for_gdx->name);
     assert(vec_for_gdx->vec);
     assert(vec_for_gdx->size > 0);
     snprintf(vec_descr, sizeof(vec_descr), "%s vector", vec_for_gdx->name);

     if ((status=NV_to_GDX(Xptr, vec_for_gdx->name, vec_descr, vec_for_gdx->vec, vec_for_gdx->size))) {
       fprintf(stderr, "Model data for vector %s not written\n", vec_for_gdx->name);
       infos[1] = (double)-ETERMINATE;
       goto fail;
     }
     DEBUG_PRINTF("FC3D_AVI_GAMS :: %s vector written\n", vec_for_gdx->name);
     vec_for_gdx = vec_for_gdx->next;

   }

  if (idxClose(Xptr))
    idxerrorR(idxGetLastError(Xptr), "idxClose");


//   cp(gdxFileName, "fc3d_avi-condensed.gdx");


  if ((status=CallGams(Gptr, Optr, sysdir, model))) {
    fprintf(stderr, "Call to GAMS failed\n");
    infos[1] = (double)-ETERMINATE;
    goto fail;
  }


  /************************************************
   * Read back solution
   ************************************************/
  idxOpenRead(Xptr, solFileName, &status);
  if (status)
    idxerrorR(status, "idxOpenRead");

  while (vec_from_gdx)
  {
    assert(vec_from_gdx->name);
    assert(vec_from_gdx->vec);
    assert(vec_from_gdx->size > 0);
    double* data = vec_from_gdx->vec;
    unsigned size = vec_from_gdx->size;
    /* GAMS does not set a value to 0 ... --xhub */
    memset(data, 0, size*sizeof(double));
    if ((status=GDX_to_NV(Xptr, vec_from_gdx->name, data, size))) {
      fprintf(stderr, "Model data %s could not be read\n", vec_from_gdx->name);
      infos[1] = (double)-ETERMINATE;
      goto fail;
    }
    vec_from_gdx = vec_from_gdx->next;
  }

  if ((status=GDX_to_NV(Xptr, "infos", infos, sizeof(infos)/sizeof(double)))) {
    fprintf(stderr, "infos could not be read\n");
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

#define GAMS_ADD_OPT(GAMSP_OPT_L, GAMSP_OPT_T) \
GAMSP_OPT_T* next_opt = GAMSP_OPT_L; \
GAMSP_OPT_T* new_opt; \
if (next_opt) \
{ \
  while (next_opt->next_opt) \
  { \
    next_opt = next_opt->next_opt; \
  } \
  next_opt->next_opt = (GAMSP_OPT_T*)malloc(sizeof(GAMSP_OPT_T)); \
  new_opt = next_opt->next_opt; \
} \
else \
{ \
  GAMSP_OPT_L = (GAMSP_OPT_T*)malloc(sizeof(GAMSP_OPT_T)); \
  new_opt = GAMSP_OPT_L; \
} \
new_opt->name = lname; \
new_opt->value = value; \
new_opt->type = type; \
new_opt->next_opt = NULL; \

#define GAMS_ADD_PREP(GP, name) \
assert(GP); \
assert(name); \
char* lname = strdup(name);

void add_GAMS_opt_str(SN_GAMSparams* GP, char* name, char* value_orig, unsigned type)
{
  GAMS_ADD_PREP(GP, name);
  assert(value_orig);
  char* value = strdup(value_orig);
  GAMS_ADD_OPT(GP->opt_str_list, GAMS_opt_str);
}

void add_GAMS_opt_bool(SN_GAMSparams* GP, char* name, bool value, unsigned type)
{
  GAMS_ADD_PREP(GP, name);
  GAMS_ADD_OPT(GP->opt_bool_list, GAMS_opt_bool);
}

void add_GAMS_opt_int(SN_GAMSparams* GP, char* name, int value, unsigned type)
{
  GAMS_ADD_PREP(GP, name);
  GAMS_ADD_OPT(GP->opt_int_list, GAMS_opt_int);
}

void add_GAMS_opt_double(SN_GAMSparams* GP, char* name, double value, unsigned type)
{
  GAMS_ADD_PREP(GP, name);
  GAMS_ADD_OPT(GP->opt_double_list, GAMS_opt_double);
}

SN_GAMSparams* createGAMSparams(char* model_dir, char* gams_dir)
{
  SN_GAMSparams* GP = (SN_GAMSparams*) malloc(sizeof(SN_GAMSparams));

  GP->model_dir = strdup(model_dir);
  GP->gams_dir = strdup(gams_dir);
  assert(GP->model_dir);
  assert(GP->gams_dir);
  GP->filename = NULL;
  GP->filename_suffix = NULL;
  GP->opt_str_list = NULL;
  GP->opt_bool_list = NULL;
  GP->opt_int_list = NULL;
  GP->opt_double_list = NULL;

  return GP;
}

void deleteGAMSparams(SN_GAMSparams* GP)
{
  if (GP->model_dir)
  {
    free(GP->model_dir);
    GP->model_dir = NULL;
  }

  if (GP->gams_dir)
  {
    free(GP->gams_dir);
    GP->gams_dir = NULL;
  }

  if (GP->opt_str_list)
  {
    GAMS_opt_str* next_opt = GP->opt_str_list;
    do 
    {
      GAMS_opt_str* str_opt = next_opt;
      next_opt = str_opt->next_opt;
      assert(str_opt->name);
      free(str_opt->name);
      str_opt->name = NULL;
      assert(str_opt->value);
      free(str_opt->value);
      str_opt->value = NULL;
      str_opt->next_opt = NULL;
      free(str_opt);
    }
    while (next_opt);
    GP->opt_str_list = NULL;
  }
  if (GP->opt_bool_list)
  {
    GAMS_opt_bool* next_opt = GP->opt_bool_list;
    do 
    {
      GAMS_opt_bool* bool_opt = next_opt;
      next_opt = bool_opt->next_opt;
      bool_opt->name = NULL;
      bool_opt->value = false;
      bool_opt->next_opt = NULL;
      free(bool_opt);
    }
    while (next_opt);
    GP->opt_bool_list = NULL;
  }
  if (GP->opt_int_list)
  {
    GAMS_opt_int* next_opt = GP->opt_int_list;
    do 
    {
      GAMS_opt_int* int_opt = next_opt;
      next_opt = int_opt->next_opt;
      int_opt->name = NULL;
      int_opt->value = 0;
      int_opt->next_opt = NULL;
      free(int_opt);
    }
    while (next_opt);
    GP->opt_int_list = NULL;
  }
  if (GP->opt_double_list)
  {
    GAMS_opt_double* next_opt = GP->opt_double_list;
    do 
    {
      GAMS_opt_double* double_opt = next_opt;
      next_opt = double_opt->next_opt;
      double_opt->name = NULL;
      double_opt->value = 0.;
      double_opt->next_opt = NULL;
      free(double_opt);
    }
    while (next_opt);
    GP->opt_double_list = NULL;
  }
  free(GP);
}

/*
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
*/
#endif
