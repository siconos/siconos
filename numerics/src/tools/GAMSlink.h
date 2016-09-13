/* GAMS stuff */

#ifndef GAMSLINK_H
#define GAMSLINK_H

#include <stdbool.h>
#include "NumericsMatrix.h"
#include <assert.h>


//#include "SolverOptions.h"

/** Simply linked list of bool option for GAMS
 */

enum { GAMS_OPT_GENERAL, GAMS_OPT_SOLVER };

typedef struct GAMS_opt_bool_ {
  char* name; /**< Name of the option */
  bool value; /**< Value of the option */
  unsigned type; /**< Type of option (general or solver-specific)*/
  struct GAMS_opt_bool_* next_opt; /**< Link to the next option*/
} GAMS_opt_bool;

/** Simply linked list of integer option for GAMS
 */
typedef struct GAMS_opt_int_ {
  char* name; /**< Name of the option */
  int value; /**< Value of the option */
  unsigned type; /**< Type of option (general or solver-specific)*/
  struct GAMS_opt_int_* next_opt; /**< Link to the next option*/
} GAMS_opt_int;

/** Simply linked list of double option for GAMS
 */
typedef struct GAMS_opt_double_ {
  char* name; /**< Name of the option */
  double value; /**< Value of the option */
  unsigned type; /**< Type of option (general or solver-specific)*/
  struct GAMS_opt_double_* next_opt; /**< Link to the next option*/
} GAMS_opt_double;

/** Simply linked list of string option for GAMS
 */
typedef struct GAMS_opt_str_ {
  char* name; /**< Name of the option */
  char* value; /**< Value of the option */
  unsigned type; /**< Type of option (general or solver-specific)*/
  struct GAMS_opt_str_* next_opt; /**< Link to the next option*/
} GAMS_opt_str;

/** Parameters for GAMS */
typedef struct {
  char* model_dir; /**<  Directory where the GAMS model are */
  char* gams_dir;  /**<  GAMS directory */
  char* filename; /**< name of the problem (used as a gdx filename) */
  char* filename_suffix; /**< suffix for the filename. Useful when solving the same problem with the same solver, but different options */
  GAMS_opt_str* opt_str_list; /**< list of string options */
  GAMS_opt_bool* opt_bool_list; /**< list of boolean options */
  GAMS_opt_int* opt_int_list; /**< list of integer options */
  GAMS_opt_double* opt_double_list; /**< list of double options */
} SN_GAMSparams;

static inline SN_GAMSparams* createGAMSparams(char* model_dir, char* gams_dir)
{
  SN_GAMSparams* GP = (SN_GAMSparams*) malloc(sizeof(SN_GAMSparams));

  GP->model_dir = model_dir;
  GP->gams_dir = gams_dir;
  GP->filename = NULL;
  GP->filename_suffix = NULL;
  GP->opt_str_list = NULL;
  GP->opt_bool_list = NULL;
  GP->opt_int_list = NULL;
  GP->opt_double_list = NULL;

  return GP;
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
new_opt->name = name; \
new_opt->value = value; \
new_opt->type = type; \
new_opt->next_opt = NULL; \

static inline void add_GAMS_opt_str(SN_GAMSparams* GP, char* name, char* value, unsigned type)
{
   GAMS_ADD_OPT(GP->opt_str_list, GAMS_opt_str);
}

static inline void add_GAMS_opt_bool(SN_GAMSparams* GP, char* name, bool value, unsigned type)
{
   GAMS_ADD_OPT(GP->opt_bool_list, GAMS_opt_bool);
}

static inline void add_GAMS_opt_int(SN_GAMSparams* GP, char* name, int value, unsigned type)
{
   GAMS_ADD_OPT(GP->opt_int_list, GAMS_opt_int);
}

static inline void add_GAMS_opt_double(SN_GAMSparams* GP, char* name, double value, unsigned type)
{
   GAMS_ADD_OPT(GP->opt_double_list, GAMS_opt_double);
}

static inline void deleteGAMSparams(SN_GAMSparams* GP)
{
  if (GP->opt_str_list)
  {
    GAMS_opt_str* next_opt = GP->opt_str_list;
    do 
    {
      GAMS_opt_str* str_opt = next_opt;
      next_opt = str_opt->next_opt;
      str_opt->name = NULL;
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

static inline const char* GAMSP_get_filename(const void* GP)
{
  return ((SN_GAMSparams*) GP)->filename;
}

static inline void GAMSP_set_filename(void* GP, char* filename)
{
  ((SN_GAMSparams*) GP)->filename = filename;
}

static inline const char* GAMSP_get_filename_suffix(const void* GP)
{
  return ((SN_GAMSparams*) GP)->filename_suffix;
}

static inline void GAMSP_set_filename_suffix(void* GP, char* filename_suffix)
{
  ((SN_GAMSparams*) GP)->filename_suffix = filename_suffix;
}


typedef struct SN_GAMS_NM_gdx_
{
  NumericsMatrix* mat;
  char* name;
  struct SN_GAMS_NM_gdx_* next;
} SN_GAMS_NM_gdx;

typedef struct SN_GAMS_NV_gdx_
{
  double* vec;
  char* name;
  unsigned size;
  struct SN_GAMS_NV_gdx_* next;
} SN_GAMS_NV_gdx;

typedef struct
{
  SN_GAMS_NM_gdx* mat_for_gdx;
  SN_GAMS_NV_gdx* vec_for_gdx;
  SN_GAMS_NV_gdx* vec_from_gdx;
} SN_GAMS_gdx;

#define SN_FREE_TILL_NEXT(X, T, ELT) \
  while(X) { T* next = X->next; X->ELT = NULL; X->name = NULL; X->next = NULL;\
    free(X); X = next; }
  

static inline void SN_free_SN_GAMS_gdx(SN_GAMS_gdx* gdx_data)
{
  assert(gdx_data);

  SN_GAMS_NM_gdx* mat_for_gdx = gdx_data->mat_for_gdx;
  SN_GAMS_NV_gdx* vec_for_gdx = gdx_data->vec_for_gdx;
  SN_GAMS_NV_gdx* vec_from_gdx = gdx_data->vec_from_gdx;
  SN_FREE_TILL_NEXT(mat_for_gdx, SN_GAMS_NM_gdx, mat);
  SN_FREE_TILL_NEXT(vec_for_gdx, SN_GAMS_NV_gdx, vec);
  SN_FREE_TILL_NEXT(vec_from_gdx, SN_GAMS_NV_gdx, vec);

}


static inline void SN_GAMS_add_NM_to_gdx(SN_GAMS_gdx* gdx_data, NumericsMatrix* M, char* name)
{
  assert(gdx_data);
  assert(M);
  assert(name);

  SN_GAMS_NM_gdx* mat_for_gdx;

  if (gdx_data->mat_for_gdx)
  {
    mat_for_gdx = gdx_data->mat_for_gdx;
    while (mat_for_gdx->next)
    {
      mat_for_gdx = mat_for_gdx->next;
    }
    mat_for_gdx->next = (SN_GAMS_NM_gdx*)malloc(sizeof(SN_GAMS_NM_gdx));

    mat_for_gdx = mat_for_gdx->next;
  }
  else
  {
    gdx_data->mat_for_gdx = (SN_GAMS_NM_gdx*)malloc(sizeof(SN_GAMS_NM_gdx));
    mat_for_gdx = gdx_data->mat_for_gdx;
  }

  mat_for_gdx->mat = M;
  mat_for_gdx->name = name;
  mat_for_gdx->next = NULL;
}

#define SN_GAMS_ADD_GDX(X) \
  SN_GAMS_NV_gdx* X; \
 \
  if (gdx_data->X) \
  { \
    X = gdx_data->X; \
    while (X->next) \
    { \
      X = X->next; \
    } \
    X->next = (SN_GAMS_NV_gdx*)malloc(sizeof(SN_GAMS_NV_gdx)); \
 \
    X = X->next; \
  } \
  else \
  { \
    gdx_data->X = (SN_GAMS_NV_gdx*)malloc(sizeof(SN_GAMS_NV_gdx)); \
    X = gdx_data->X; \
  } \
 \
  X->vec = vec; \
  X->name = name; \
  X->size = size; \
  X->next = NULL;

static inline void SN_GAMS_add_NV_to_gdx(SN_GAMS_gdx* gdx_data, double* vec, char* name, unsigned size)
{
  assert(gdx_data);
  assert(vec);
  assert(name);

  SN_GAMS_ADD_GDX(vec_for_gdx)
}

static inline void SN_GAMS_add_NV_from_gdx(SN_GAMS_gdx* gdx_data, double* vec, char* name, unsigned size)
{
  assert(gdx_data);
  assert(vec);
  assert(name);

  SN_GAMS_ADD_GDX(vec_from_gdx)
}

#ifdef HAVE_GAMS_C_API

#include "gclgms.h"
#include "gamsxcc.h"
#include "idxcc.h"
#include "optcc.h"
#include "gevmcc.h"
#include "gmomcc.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define idxerror(i, s) { idxErrorStr(Xptr, i, msg, GMS_SSSIZE); \
  printf("%s failed: %s\n",s,msg); return; }

#define idxerrorR(i, s) { idxErrorStr(Xptr, i, msg, GMS_SSSIZE); \
  printf("%s failed: %s\n",s,msg); return 1; }


#define STR_VALUE(arg)      #arg
#define SPACE_CONC(str1,str2) str1 " " str2

static inline void SN_Gams_set_dirs(const SN_GAMSparams* solverParameters, const char* defModel, const char* defGAMSdir, char* model, char* GAMSdir, char* model_file)
{
  assert(solverParameters);
  assert(defModel);
  assert(defGAMSdir);
  assert(model);
  assert(GAMSdir);
  assert(model_file);

  if (solverParameters->model_dir)
  {
    size_t len1 = strlen(solverParameters->model_dir);
    strncpy(model, solverParameters->model_dir, len1);
    strncpy(&model[len1], model_file, GMS_SSSIZE-len1-2);
    model[GMS_SSSIZE-1] = '\0';
  }
  else
  {
    strncpy(model, defModel, GMS_SSSIZE-2);
    model[GMS_SSSIZE-1] = '\0';
  }

  if (solverParameters->gams_dir)
  {
    strncpy(GAMSdir, solverParameters->gams_dir, GMS_SSSIZE-2);
    GAMSdir[GMS_SSSIZE-1] = '\0';
  }
  else
  {
    strncpy(GAMSdir, defGAMSdir, GMS_SSSIZE-2);
    GAMSdir[GMS_SSSIZE-1] = '\0';
  }
}

#define WALK_GAMSP_OPTS(GAMSP_OPT_L, GAMSP_OPT_T, GAMS_OPT_FUN) \
  if (GAMSP_OPT_L) \
  { \
    GAMSP_OPT_T* next_opt = GAMSP_OPT_L; \
    do \
    { \
      GAMSP_OPT_T* opt = next_opt; \
      next_opt = opt->next_opt; \
      GAMS_OPT_FUN(Opts[opt->type], opt->name, opt->value); \
    } \
    while (next_opt); \
  }

static inline void SN_Gams_set_options(const SN_GAMSparams* GP, optHandle_t* Opts)
{
  assert(GP);
  assert(Opts);
  WALK_GAMSP_OPTS(GP->opt_str_list, GAMS_opt_str, optSetStrStr);
  WALK_GAMSP_OPTS(GP->opt_bool_list, GAMS_opt_bool, optSetIntStr);
  WALK_GAMSP_OPTS(GP->opt_int_list, GAMS_opt_int, optSetIntStr);
  WALK_GAMSP_OPTS(GP->opt_double_list, GAMS_opt_double, optSetDblStr);
}


static inline int getGamsSolverOpt(const optHandle_t Optr, const char* sysdir, const char* solverDefName)
{
  assert(Optr);
  assert(sysdir);
  assert(solverDefName);

  char deffile[GMS_SSSIZE];
  char msg[GMS_SSSIZE];
  strncpy(deffile, sysdir, sizeof(deffile));
  strncat(deffile, "/opt", sizeof(deffile) - strlen(deffile) - 1);
  strncat(deffile, solverDefName, sizeof(deffile) - strlen(deffile) - 1);
  strncat(deffile, ".def", sizeof(deffile) - strlen(deffile) - 1);

  if (optReadDefinition(Optr,deffile)) {
    int itype;
    for (int i=1; i<=optMessageCount(Optr); ++i) {
      optGetMessage(Optr, i, msg, &itype);
      printf("%s\n", msg);
    }
    return 1;
  }
  return 0;
}
static inline int getGamsOpt(const optHandle_t Optr, const char *sysdir)
{
  char msg[GMS_SSSIZE];
  char deffile[GMS_SSSIZE];
  strncpy(deffile, sysdir, sizeof(deffile));
  strncat(deffile, "/optgams.def", sizeof(deffile) - strlen(deffile) - 1);

  if (optReadDefinition(Optr,deffile)) {
    int itype;
    for (int i=1; i<=optMessageCount(Optr); ++i) {
      optGetMessage(Optr, i, msg, &itype);
      printf("%s\n", msg);
    }
    return 1;
  }
  optSetStrStr(Optr, "sysdir", sysdir);

  return 0;
}

static inline int CallGams(const gamsxHandle_t Gptr, const optHandle_t Optr, const char *sysdir, const char *model)
{
  char msg[GMS_SSSIZE];

  assert(Gptr); assert(Optr);

  optSetStrStr(Optr, "input", model);
  optSetIntStr(Optr, "logoption", 2);
//  optSetIntStr(Optr, "keep", 1);
  optSetIntStr(Optr, "optfile", 1);
//  optSetDblStr(Optr,"OptCA", 1e-12);

  if (gamsxRunExecDLL(Gptr, Optr, sysdir, 1, msg)) {
    printf ("Could not execute RunExecDLL: %s\n", msg);
    return 1;
  }

  return 0;
}

static inline int iparam_to_GDX(idxHandle_t Xptr, const char* name, const char* descr, double param)
{
  char msg[GMS_SSSIZE];
  int dim = 1;

  if (idxDataWriteStart(Xptr, name, descr, 0, &dim, msg, GMS_SSSIZE) == 0)
    idxerrorR(idxGetLastError(Xptr), "idxDataWriteStart");

  idxDataWrite(Xptr, 0, param);

  if (0==idxDataWriteDone(Xptr))
    idxerrorR(idxGetLastError(Xptr), "idxDataWriteDone");

  return 0;
}


static inline int NV_to_GDX(idxHandle_t Xptr, const char* name, const char* descr, const double* vector, unsigned size)
{
  char msg[GMS_SSSIZE];

  int dim = size;
  if (idxDataWriteStart(Xptr, name, descr, 1, &dim, msg, GMS_SSSIZE) == 0)
    idxerrorR(idxGetLastError(Xptr), "idxDataWriteStart");

  idxDataWriteDenseColMajor(Xptr, 1, vector);

  if (0==idxDataWriteDone(Xptr))
    idxerrorR(idxGetLastError(Xptr), "idxDataWriteDone");

  return 0;
}

static inline int NM_to_GDX(idxHandle_t Xptr, const char* name, const char* descr, NumericsMatrix* M)
{
  char msg[GMS_SSSIZE];

  int dims[2];
  dims[0] = M->size0;
  dims[1] = M->size1;
  if (idxDataWriteStart(Xptr, name, descr, 2, dims, msg, GMS_SSSIZE) == 0)
    idxerrorR(idxGetLastError(Xptr), "idxDataWriteStart");

  if (M->storageType == NM_DENSE)
  {
    assert(M->matrix0);
    idxDataWriteDenseColMajor(Xptr, 2, M->matrix0);
  }
  else
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
  }


  if (0==idxDataWriteDone(Xptr))
    idxerrorR(idxGetLastError(Xptr), "idxDataWriteDone");

  return 0;
}

static inline int GDX_to_NV(idxHandle_t Xptr, const char* name, double* vector, unsigned size)
{
  char msg[GMS_SSSIZE];

  int nbdims, nbelts;
  int dims[GLOBAL_MAX_INDEX_DIM];
  if (idxDataReadStart(Xptr, name, &nbdims, dims, &nbelts, msg, GMS_SSSIZE) == 0)
    idxerrorR(idxGetLastError(Xptr), "idxDataReadStart");

  if (nbdims != 1 || dims[0] != (int)size)
  {
    printf("GDX_to_NV :: inconsistency between expected size and actual one, variable %s\n", name);
    printf("expected dimension: %d; actual one: %d\n", 1, nbdims);
    printf("expected number of elements: %d; actual one: %d\n", size, nbelts);
  }
  idxDataReadDenseColMajor(Xptr, vector);

  if (0==idxDataReadDone(Xptr))
    idxerrorR(idxGetLastError(Xptr), "idxDataWriteDone");

  return 0;
}


#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  int SN_gams_solve(unsigned iter, optHandle_t Optr, char* sysdir, char* model, const char* base_name, SolverOptions* options, SN_GAMS_gdx* gdx_data);

  void filename_datafiles(const int iter, const int solverId, const char* base_name, unsigned len, char* template_name, char* log_filename);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#else

#endif

#endif /* GAMSLINK_H  */
