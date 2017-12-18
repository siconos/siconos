
%apply (double *b) { (double *b_bck_gfc3d) };

// for b in GlobalFrictionContact
%typemap(out) (double* b) {

  if (!arg1->H) { SWIG_exception_fail(SWIG_RuntimeError, "H is not present, don't known the size"); }

  if ($1)
  {
    SN_OBJ_TYPE *obj;
    C_to_target_lang1(obj, arg1->H->size1, $1, SWIG_fail);
    $result = obj;
  }
  else
    SWIG_fail;
}

// vectors of size col(H)
%typemap(memberin) (double *b) {
  // Still some dark magic :( --xhub
 char msg[1024];
  assert(arg1);
  if (!arg1->H)
  {
    SWIG_exception_fail(SWIG_RuntimeError, "H is not initialized, it sould be done first!");
  }

  int size = arg1->H->size1;
  if (size !=  array_size((SN_ARRAY_TYPE*)array2, 0))
  {
    snprintf(msg, sizeof(msg), "Size of b is %ld, but the size of H is %d! Both should be equal!\n", array_size((SN_ARRAY_TYPE*)array2, 0), size);
    SWIG_exception_fail(SWIG_RuntimeError, msg);
  }

  if (!$1) { $1 = (double*)malloc(size * sizeof(double)); }
  memcpy($1, $input, size * sizeof(double));

}

%{
#include "GlobalFrictionContactProblem.h"
#include "GlobalFrictionContactProblem.h"
#include "gfc3d_Solvers.h"
#include "gfc3d_compute_error.h"
  
%}

%include "GlobalFrictionContactProblem.h"
%include "gfc3d_Solvers.h"
%include "gfc3d_compute_error.h"

%extend GlobalFrictionContactProblem
{

  GlobalFrictionContactProblem()
    {

      GlobalFrictionContactProblem *FC;
      // return pointer : free by std swig destructor
      FC = (GlobalFrictionContactProblem *) malloc(sizeof(GlobalFrictionContactProblem));

      globalFrictionContact_null(FC);

      return FC;
    }

  GlobalFrictionContactProblem(SN_OBJ_TYPE *dim)
    {

      GlobalFrictionContactProblem *FC;
      // return pointer : free by std swig destructor
      FC = (GlobalFrictionContactProblem *) malloc(sizeof(GlobalFrictionContactProblem));

      globalFrictionContact_null(FC);
      SWIG_AsVal_int(dim, &FC->dimension);

      return FC;
    }

  GlobalFrictionContactProblem(SN_OBJ_TYPE *dim, SN_OBJ_TYPE *o1, SN_OBJ_TYPE *o2, SN_OBJ_TYPE *o3)
    {

      int is_new_object1=0;
      int is_new_object2=0;
      int is_new_object3=0;

      SN_ARRAY_TYPE* array = obj_to_sn_array(o1, &is_new_object1);
      SN_ARRAY_TYPE* vector = obj_to_sn_vector(o2, &is_new_object2);
      SN_ARRAY_TYPE* mu_vector = obj_to_sn_vector(o3, &is_new_object3);
      GlobalFrictionContactProblem * FC = (GlobalFrictionContactProblem *) malloc(sizeof(GlobalFrictionContactProblem));
      globalFrictionContact_null(FC);

      size_t size0 = array_size((SN_ARRAY_TYPE*)array,0);
      size_t size1 = array_size((SN_ARRAY_TYPE*)array,1);
      FC->M = NM_create(NM_DENSE, size0, size1);

      memcpy(FC->M->matrix0,array_data(array),size0*size1*sizeof(double));
      SWIG_AsVal_int(dim, &FC->dimension);
      FC->numberOfContacts = size0 / FC->dimension;
      FC->q = (double *) malloc(size0*sizeof(double));
      memcpy(FC->q,array_data(vector),size0*sizeof(double));
      FC->mu = (double *) malloc(FC->numberOfContacts*sizeof(double));
      memcpy(FC->mu,array_data(mu_vector),FC->numberOfContacts*sizeof(double));

      // python mem management
      target_mem_mgmt(is_new_object1, array)
      target_mem_mgmt(is_new_object2, vector)
      target_mem_mgmt(is_new_object3, mu_vector)

      return FC;
    }

  ~GlobalFrictionContactProblem()
  {
    freeGlobalFrictionContactProblem($self);
  }

};


// remove typemaps for b
%clear double* b;
%apply (double *b_bck_gfc3d) { (double *b) };

