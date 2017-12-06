%{
#include "fc3d_Solvers.h"
#include "Friction_cst.h"
#include "FrictionContactProblem.h"
#include "fc3d_compute_error.h"
#ifdef WITH_FCLIB
// avoid a conflict with old csparse.h in case fclib.h includes it
#define _CS_H
#include "fclib_interface.h"
#endif
#include "fc2d_Solvers.h"
#include "gfc3d_Solvers.h"
#include "fc3d_AlartCurnier_functions.h"
#include "AlartCurnierGenerated.h"
#include "FischerBurmeisterGenerated.h"
#include "NaturalMapGenerated.h"
%}

%include "FrictionContactProblem.h"
#ifdef WITH_FCLIB
// avoid a conflict with old csparse.h in case fclib.h includes it
#define _CS_H
%include fclib_interface.h
#endif

%include "fc3d_AlartCurnier_functions.h"
%include "fc3d_nonsmooth_Newton_AlartCurnier.h"
%include "fc3d_nonsmooth_Newton_FischerBurmeister.h"
%include "fc3d_nonsmooth_Newton_natural_map.h"
%include "AlartCurnierGenerated.h"
%include "FischerBurmeisterGenerated.h"
%include "NaturalMapGenerated.h"
%include "fclib_interface.h"
%include "fc3d_Solvers.h"
%include "fc3d_unitary_enumerative.h"
%include "fc2d_Solvers.h"
%include "Friction_cst.h"
%include "fc3d_compute_error.h"

%extend FrictionContactProblem
{
  FrictionContactProblem()
  {
    return newFCP();
  }


  /* copy constructor */
  FrictionContactProblem(SN_OBJ_TYPE *o)
  {
    FrictionContactProblem* fcp;
    FrictionContactProblem* FCP;

    %SN_INPUT_CHECK_RETURN(o, fcp, FrictionContactProblem);

    FCP = newFCP();
    FCP->dimension = fcp->dimension;

    FCP->M = NM_create(fcp->M->storageType, fcp->M->size0, fcp->M->size1);

    NM_copy(FCP->M, fcp->M);

    FCP->numberOfContacts = fcp->numberOfContacts;
    size_t size =  fcp->dimension * fcp->numberOfContacts;

    FCP->q = (double*) malloc(size * sizeof(double));
    memcpy(FCP->q, fcp->q, size * sizeof(double));

    FCP->q = (double*) malloc(size * sizeof(double));
    memcpy(FCP->mu, fcp->mu, size * sizeof(double));

    return FCP;
  }

  /* */
  FrictionContactProblem(SN_OBJ_TYPE *dim, SN_OBJ_TYPE *numberOfContacts, SN_OBJ_TYPE *M, SN_OBJ_TYPE *q, SN_OBJ_TYPE *mu)
  {
    FrictionContactProblem * FC = newFCP();
    SWIG_AsVal_int(dim, &FC->dimension);
    SWIG_AsVal_int(numberOfContacts, &FC->numberOfContacts);

    %NM_convert_from_target(M, (&FC->M), return NULL);

    int is_new_objectq=0;
    SN_ARRAY_TYPE* vector = obj_to_sn_vector(q, &is_new_objectq);

    sn_check_size_mat_vec(FC->M->size0, vector, return NULL);
    set_vec_from_target(FC->q, vector, , return NULL);

    int is_new_objectmu = 0;
    SN_ARRAY_TYPE* vectormu = obj_to_sn_vector(mu, &is_new_objectmu);

    sn_check_size_mat_vec(FC->numberOfContacts, vectormu, return NULL);
    set_vec_from_target(FC->mu, vectormu, , return NULL);

    target_mem_mgmt(is_new_objectq, vector);
    target_mem_mgmt(is_new_objectmu, vectormu);

    return FC;

  }

  FrictionContactProblem(SN_OBJ_TYPE *dim, SN_OBJ_TYPE * M, SN_OBJ_TYPE * q, SN_OBJ_TYPE * mu)
  {
    FrictionContactProblem * FC = newFCP();

    SWIG_AsVal_int(dim, &FC->dimension);

    %NM_convert_from_target(M, (&FC->M), return NULL);

    FC->numberOfContacts = FC->M->size0 / FC->dimension;

    int is_new_objectq=0;
    SN_ARRAY_TYPE* vector = obj_to_sn_vector(q, &is_new_objectq);

    sn_check_size_mat_vec(FC->M->size0, vector, return NULL);
    set_vec_from_target(FC->q, vector, , return NULL);

    int is_new_objectmu = 0;
    SN_ARRAY_TYPE* vectormu = obj_to_sn_vector(mu, &is_new_objectmu);

    sn_check_size_mat_vec(FC->numberOfContacts, vectormu, return NULL);
    set_vec_from_target(FC->mu, vectormu, , return NULL);

    target_mem_mgmt(is_new_objectq, vector);
    target_mem_mgmt(is_new_objectmu, vectormu);


    return FC;
  }

  ~FrictionContactProblem()
  {
    freeFrictionContactProblem($self);
  }

};

%inline %{

#include <stdio.h>
  static FrictionContactProblem* frictionContactProblemFromFile
    (const char * filename)
  {
    FILE * finput = fopen(filename, "r");
    if (finput)
    {
      FrictionContactProblem* problem =
        (FrictionContactProblem *) malloc(sizeof(FrictionContactProblem));
      if (frictionContact_newFromFile(problem,finput))
      {
      char msg[1024];
      snprintf(msg, sizeof(msg), "frictionContactProblemFromFile: cannot load %s\n",filename);
      SWIG_Error(SWIG_RuntimeError, msg);
      free(problem);
      fclose(finput);
      return NULL;
      }
      else
      {
        fclose(finput);
        return problem;
      }
    }
    else
    {
      char msg[1024];
      snprintf(msg, sizeof(msg), "frictionContactProblemFromFile: cannot open %s\n",filename);
      SWIG_Error(SWIG_RuntimeError, msg);
      return NULL;
    }

  }

%}
