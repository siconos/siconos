// LCP

%ignore lcp_compute_error_only;

%{
  #include "LinearComplementarityProblem.h"
  #include "LCP_Solvers.h"
  #include "lcp_cst.h"
  %}

%include "LinearComplementarityProblem.h"
%include "LCP_Solvers.h"
%include "lcp_cst.h"

%extend LinearComplementarityProblem
{
  LinearComplementarityProblem(SN_OBJ_TYPE *M, SN_OBJ_TYPE *q)
    {

      int is_new_object2=0;
      SN_ARRAY_TYPE* vector = obj_to_sn_vector(q, &is_new_object2);

      LinearComplementarityProblem* LC = newLCP();
      %NM_convert_from_target(M, (&LC->M), return NULL);
      int size0 = LC->M->size0;
      sn_check_array_type(vector, return NULL);
      sn_check_size_mat_vec(LC->M->size0, vector, return NULL);

      LC->size = size0;
      set_vec_from_target(LC->q, vector, NULL, return NULL);

      // python mem management
      target_mem_mgmt(is_new_object2, vector);

      return LC;

    }


  ~LinearComplementarityProblem()
  {
    freeLinearComplementarityProblem($self);
  }

};


