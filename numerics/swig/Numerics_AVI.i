%{
#include "AffineVariationalInequalities.h"
#include "AVI_Solvers.h"
#include "AVI_cst.h"
#include "VariationalInequality_Solvers.h"
%}


%include "AffineVariationalInequalities.h"
%include "AVI_Solvers.h"
%include "AVI_cst.h"

%extend AffineVariationalInequalities
{
  AffineVariationalInequalities()
   {
     AffineVariationalInequalities* avi = newAVI();
     return avi;
   }


  AffineVariationalInequalities(SN_OBJ_TYPE* mat, SN_OBJ_TYPE* vec)
  {


     AffineVariationalInequalities* avi = newAVI();

     int is_new_object2=0;
     SN_ARRAY_TYPE* vector = obj_to_sn_vector(vec, &is_new_object2); 

     %NM_convert_from_target(mat, (&avi->M), return NULL);

     if (avi->M->size0 != avi->M->size1)
     {
       SWIG_Error(SWIG_ValueError, "A non square matrix has been given");
       freeAVI(avi);
       return NULL;
     }

     sn_check_array_type(vector, return NULL);
     sn_check_size_mat_vec(avi->M->size0, vector, return NULL);

     avi->size = avi->M->size0;
     set_vec_from_target(avi->q, vector, , return NULL);

      // python mem management
      target_mem_mgmt(is_new_object2, vector);

     return avi;
   }

   void set_polyhedron(SN_OBJ_TYPE* H_mat, SN_OBJ_TYPE* K_vec)
   {
     $self->poly.split = (polyhedron*) malloc(sizeof(polyhedron));
      int is_new_object2=0;
      %NM_convert_from_target(H_mat, (&$self->poly.split->H), TARGET_ERROR_VERBOSE);

      if ($self->poly.split->H->size1 != $self->size)
      {
        SWIG_Error(SWIG_TypeError, "The matrix does not have the right number of column");
        TARGET_ERROR_VERBOSE;
      }

      SN_ARRAY_TYPE* vector = obj_to_sn_vector(K_vec, &is_new_object2); 
      sn_check_array_type(vector, TARGET_ERROR_VERBOSE);
      sn_check_size_mat_vec($self->poly.split->H->size0, vector, TARGET_ERROR_VERBOSE);

      set_vec_from_target($self->poly.split->K, vector, , TARGET_ERROR_VERBOSE);

      $self->poly.split->size_ineq = $self->poly.split->H->size0;
      $self->poly.split->size_eq = 0;
      $self->poly.split->Heq = NULL;
      $self->poly.split->Keq = NULL;

   }

  ~AffineVariationalInequalities()
  {
    freeAVI($self);
  }
};


