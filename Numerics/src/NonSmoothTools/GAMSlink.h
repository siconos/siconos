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

/** Parameters for GAMS */
typedef struct {
  char* model_dir; /**<  Directory where the GAMS model are */
  char* gams_dir;  /**<  GAMS directory */
} SN_GAMSparams;



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
  strncat(deffile, "/opt", sizeof(deffile));
  strncat(deffile, solverDefName, sizeof(deffile));
  strncat(deffile, ".def", sizeof(deffile));

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
  strncat(deffile, "/optgams.def", sizeof(deffile));

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
  optSetIntStr(Optr, "logoption", 4);
  optSetIntStr(Optr, "keep", 1);
  optSetIntStr(Optr, "optfile", 1);
//  optSetDblStr(Optr,"OptCA", 1e-12);

  if (gamsxRunExecDLL(Gptr, Optr, sysdir, 1, msg)) {
    printf ("Could not execute RunExecDLL: %s", msg);
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
    idxDataWriteSparseColMajor(Xptr, cs->p, cs->i, cs->x);
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
