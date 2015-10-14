/* GAMS stuff */

#include <stdbool.h>

/** Simply linked list of bool option for GAMS
 */

typedef struct GAMS_opt_bool_ {
  char* name; /**< Name of the option */
  bool value; /**< Value of the option */
  struct GAMS_opt_bool_* next_opt; /**< Link to the next option*/
} GAMS_opt_bool;

/** Simply linked list of integer option for GAMS
 */
typedef struct GAMS_opt_int_ {
  char* name; /**< Name of the option */
  int value; /**< Value of the option */
  struct GAMS_opt_int_* next_opt; /**< Link to the next option*/
} GAMS_opt_int;

/** Simply linked list of double option for GAMS
 */
typedef struct GAMS_opt_double_ {
  char* name; /**< Name of the option */
  double value; /**< Value of the option */
  struct GAMS_opt_double_* next_opt; /**< Link to the next option*/
} GAMS_opt_double;

/** Simply linked list of string option for GAMS
 */
typedef struct GAMS_opt_str_ {
  char* name; /**< Name of the option */
  char* value; /**< Value of the option */
  struct GAMS_opt_str_* next_opt; /**< Link to the next option*/
} GAMS_opt_str;

/** Parameters for GAMS */
typedef struct {
  char* model_dir; /**<  Directory where the GAMS model are */
  char* gams_dir;  /**<  GAMS directory */
  char* filename; /**< name of the problem (used as a gdx filename) */
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
  GP->opt_str_list = NULL;
  GP->opt_bool_list = NULL;
  GP->opt_int_list = NULL;
  GP->opt_double_list = NULL;

  return GP;
}

static inline void add_GAMS_opt_str(SN_GAMSparams* GP, char* name, char* value)
{
   GAMS_opt_str* next_opt = GP->opt_str_list;
   GAMS_opt_str* new_opt;
   if (next_opt)
   {
     do
     {
       next_opt = next_opt->next_opt;
     }
     while (next_opt->next_opt);
     next_opt->next_opt = (GAMS_opt_str*)malloc(sizeof(GAMS_opt_str));
     new_opt = next_opt->next_opt;
   }
   else
   {
     next_opt = (GAMS_opt_str*)malloc(sizeof(GAMS_opt_str));
     new_opt = next_opt;
   }
   new_opt->name = name;
   new_opt->value = value;
   new_opt->next_opt = NULL;
}

static inline void add_GAMS_opt_bool(SN_GAMSparams* GP, char* name, bool value)
{
   GAMS_opt_bool* next_opt = GP->opt_bool_list;
   GAMS_opt_bool* new_opt;
   if (next_opt)
   {
     do
     {
       next_opt = next_opt->next_opt;
     }
     while (next_opt->next_opt);
     next_opt->next_opt = (GAMS_opt_bool*)malloc(sizeof(GAMS_opt_bool));
     new_opt = next_opt->next_opt;
   }
   else
   {
     next_opt = (GAMS_opt_bool*)malloc(sizeof(GAMS_opt_bool));
     new_opt = next_opt;
   }
   new_opt->name = name;
   new_opt->value = value;
   new_opt->next_opt = NULL;
}

static inline void add_GAMS_opt_int(SN_GAMSparams* GP, char* name, int value)
{
   GAMS_opt_int* next_opt = GP->opt_int_list;
   GAMS_opt_int* new_opt;
   if (next_opt)
   {
     do
     {
       next_opt = next_opt->next_opt;
     }
     while (next_opt->next_opt);
     next_opt->next_opt = (GAMS_opt_int*)malloc(sizeof(GAMS_opt_int));
     new_opt = next_opt->next_opt;
   }
   else
   {
     next_opt = (GAMS_opt_int*)malloc(sizeof(GAMS_opt_int));
     new_opt = next_opt;
   }
   new_opt->name = name;
   new_opt->value = value;
   new_opt->next_opt = NULL;
}

static inline void add_GAMS_opt_double(SN_GAMSparams* GP, char* name, double value)
{
   GAMS_opt_double* next_opt = GP->opt_double_list;
   GAMS_opt_double* new_opt;
   if (next_opt)
   {
     do
     {
       next_opt = next_opt->next_opt;
     }
     while (next_opt->next_opt);
     next_opt->next_opt = (GAMS_opt_double*)malloc(sizeof(GAMS_opt_double));
     new_opt = next_opt->next_opt;
   }
   else
   {
     next_opt = (GAMS_opt_double*)malloc(sizeof(GAMS_opt_double));
     new_opt = next_opt;
   }
   new_opt->name = name;
   new_opt->value = value;
   new_opt->next_opt = NULL;
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
  }
  free(GP);
}

static inline const char* GAMSP_get_filename(const SN_GAMSparams* GP)
{
  return GP->filename;
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
#include <assert.h>
#include <stdlib.h>

#include "NumericsMatrix.h"

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

  if (M->storageType == 0)
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


#else

#endif
