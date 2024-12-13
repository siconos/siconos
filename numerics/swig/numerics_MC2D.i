%{
#include "mc2d_solvers.h"
#include "Plasticity_cst.h"
#include "MohrCoulomb2DProblem.h"
#include "mc2d_compute_error.h"
#ifdef WITH_FCLIB
// avoid a conflict with old csparse.h in case fclib.h includes it
#define _CS_H
#include "fclib_interface.h"
#endif
#include "mc2d_solvers.h"
%}

%include "MohrCoulomb2DProblem.h"
#ifdef WITH_FCLIB
// avoid a conflict with old csparse.h in case fclib.h includes it
#define _CS_H
%include fclib_interface.h
#endif

%include "mc2d_solvers.h"
%include "Plasticity_cst.h"
%include "mc2d_compute_error.h"

%extend MohrCoulomb2DProblem
{
  MohrCoulomb2DProblem()
  {
    return mohrCoulomb2DProblem_new();
  }


  /* copy constructor */
  MohrCoulomb2DProblem(SN_OBJ_TYPE *o)
  {
    MohrCoulomb2DProblem* mc2d;
    MohrCoulomb2DProblem* MC2D;

    %SN_INPUT_CHECK_RETURN(o, mc2d, MohrCoulomb2DProblem);

    MC2D = mohrCoulomb2DProblem_new();
    MC2D->dimension = mc2d->dimension;

    MC2D->M = NM_create(mc2d->M->storageType, mc2d->M->size0, mc2d->M->size1);

    NM_copy(MC2D->M, mc2d->M);

    MC2D->numberOfCones = mc2d->numberOfCones;
    size_t size =  mc2d->dimension * mc2d->numberOfCones;

    MC2D->q = (double*) malloc(size * sizeof(double));
    memcpy(MC2D->q, mc2d->q, size * sizeof(double));

    MC2D->eta = (double*) malloc(size * sizeof(double));
    memcpy(MC2D->theta, mc2d->eta, size * sizeof(double));
    
    MC2D->theta = (double*) malloc(size * sizeof(double));
    memcpy(MC2D->theta, mc2d->theta, size * sizeof(double));

    return MC2D;
  }

  /* */
  MohrCoulomb2DProblem(SN_OBJ_TYPE *dim, SN_OBJ_TYPE *numberOfCones, SN_OBJ_TYPE *M, SN_OBJ_TYPE *q, SN_OBJ_TYPE *eta, SN_OBJ_TYPE *theta)
  {
    MohrCoulomb2DProblem * MC2D = mohrCoulomb2DProblem_new();
    SWIG_AsVal_int(dim, &MC2D->dimension);
    SWIG_AsVal_int(numberOfCones, &MC2D->numberOfCones);

    %NM_convert_from_target(M, (&MC2D->M), return NULL);

    int is_new_objectq=0;
    SN_ARRAY_TYPE* vector = obj_to_sn_vector(q, &is_new_objectq);

    sn_check_size_mat_vec(MC2D->M->size0, vector, return NULL);
    set_vec_from_target(MC2D->q, vector, , return NULL);

    int is_new_objecteta = 0;
    SN_ARRAY_TYPE* vectoreta = obj_to_sn_vector(eta, &is_new_objecteta);

    sn_check_size_mat_vec(MC2D->numberOfCones, vectoreta, return NULL);
    set_vec_from_target(MC2D->eta, vectoreta, , return NULL);
    
    int is_new_objecttheta = 0;
    SN_ARRAY_TYPE* vectortheta = obj_to_sn_vector(theta, &is_new_objecttheta);

    sn_check_size_mat_vec(MC2D->numberOfCones, vectortheta, return NULL);
    set_vec_from_target(MC2D->theta, vectortheta, , return NULL);

    target_mem_mgmt(is_new_objectq, vector);
    target_mem_mgmt(is_new_objecteta, vectoreta);
    target_mem_mgmt(is_new_objecttheta, vectortheta);

    return MC2D;

  }

  MohrCoulomb2DProblem(SN_OBJ_TYPE *dim, SN_OBJ_TYPE * M, SN_OBJ_TYPE * q, SN_OBJ_TYPE *eta, SN_OBJ_TYPE * theta)
  {
    MohrCoulomb2DProblem * MC2D = mohrCoulomb2DProblem_new();

    SWIG_AsVal_int(dim, &MC2D->dimension);

    %NM_convert_from_target(M, (&MC2D->M), return NULL);

    MC2D->numberOfCones = MC2D->M->size0 / MC2D->dimension;

    int is_new_objectq=0;
    SN_ARRAY_TYPE* vector = obj_to_sn_vector(q, &is_new_objectq);

    sn_check_size_mat_vec(MC2D->M->size0, vector, return NULL);
    set_vec_from_target(MC2D->q, vector, , return NULL);

    int is_new_objecteta = 0;
    SN_ARRAY_TYPE* vectoreta = obj_to_sn_vector(eta, &is_new_objecteta);

    sn_check_size_mat_vec(MC2D->numberOfCones, vectoreta, return NULL);
    set_vec_from_target(MC2D->eta, vectoreta, , return NULL);
    
    int is_new_objecttheta = 0;
    SN_ARRAY_TYPE* vectortheta = obj_to_sn_vector(theta, &is_new_objecttheta);

    sn_check_size_mat_vec(MC2D->numberOfCones, vectortheta, return NULL);
    set_vec_from_target(MC2D->theta, vectortheta, , return NULL);

    target_mem_mgmt(is_new_objectq, vector);
    target_mem_mgmt(is_new_objecteta, vectoreta);
    target_mem_mgmt(is_new_objecttheta, vectortheta);

    
    return MC2D;
  }

  ~MohrCoulomb2DProblem()
  {
    mohrCoulomb2DProblem_free($self);
  }

};

%inline %{

#include <stdio.h>
  static MohrCoulomb2DProblem* mohrCoulomb2DProblemFromFile
    (const char * filename)
  {
    return mohrCoulomb2D_new_from_filename(filename);
  }

%}
