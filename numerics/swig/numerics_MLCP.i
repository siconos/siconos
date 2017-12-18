// MLCP
%apply (double *q) { (double *qbck) };

// redefine typemap on q for MLCP
%typemap(out) (double* q) {

  if ($1)
  {
    SN_OBJ_TYPE *obj;
    C_to_target_lang1(obj, arg1->m + arg1->n, $1, SWIG_fail);
    $result = obj;
  }
  else
    SWIG_fail;
 }

%{
  
  #include "mlcp_cst.h"
  #include "MixedLinearComplementarityProblem.h"
  #include "MLCP_Solvers.h"
  #include "SiconosCompat.h"
  
  %}

%include "mlcp_cst.h"
%include "MixedLinearComplementarityProblem.h"
%include "MLCP_Solvers.h"

%extend MixedLinearComplementarityProblem
{
  MixedLinearComplementarityProblem()
   {
     MixedLinearComplementarityProblem* MLCP = newMLCP();
     return MLCP;
   }

  MixedLinearComplementarityProblem(SN_OBJ_TYPE *dim, SN_OBJ_TYPE *o1, SN_OBJ_TYPE *o2)
    {

      int is_new_object2=0;
      SN_ARRAY_TYPE* vector = obj_to_sn_vector(o2, &is_new_object2);

      MixedLinearComplementarityProblem *MLCP = newMLCP();
      // return pointer : free by std swig destructor

      %NM_convert_from_target(o1, (&MLCP->M), return NULL);

      if (MLCP->M->size0 != MLCP->M->size1)
      {
        SWIG_Error(SWIG_ValueError, "A non square matrix has been given");
        freeMixedLinearComplementarityProblem(MLCP);
        return NULL;
      }

      sn_check_array_type(vector, return NULL);
      sn_check_size_mat_vec(MLCP->M->size0, vector, return NULL);


      SWIG_AsVal_int(dim, &MLCP->n);
      MLCP->m = MLCP->M->size0 - MLCP->n;
      MLCP->blocksRows = (int *) malloc(3*sizeof(int));
      MLCP->blocksIsComp = (int *) malloc(2*sizeof(int));


      MLCP->blocksRows[0]=0;
      MLCP->blocksRows[1]=MLCP->n;
      MLCP->blocksRows[2]=MLCP->n+MLCP->m;
      MLCP->blocksIsComp[0]=0;
      MLCP->blocksIsComp[1]=1;


      MLCP->isStorageType1 = 1;
      MLCP->isStorageType2 = 0;
      MLCP->A = NULL;
      MLCP->B = NULL;
      MLCP->C = NULL;
      MLCP->D = NULL;
      MLCP->a = NULL;
      MLCP->b = NULL;

      if (MLCP->M->size0 !=  array_size(vector,0))
      {
        SWIG_Error(SWIG_ValueError, "Matrix and vector of incompatible lengths");
        return NULL;
      }
      set_vec_from_target(MLCP->q, vector, , return NULL);

      // python mem management
      target_mem_mgmt(is_new_object2, vector);

      return MLCP;

    }



  ~MixedLinearComplementarityProblem()
  {
    freeMixedLinearComplementarityProblem($self);
  }

  // MixedLinearComplementarityProblem * newFromFilename(SN_OBJ_TYPE * o1)
  // {
  //   int result;
  //   MixedLinearComplementarityProblem *MLCP;
  //   // return pointer : free by std swig destructor
  //   MLCP = (MixedLinearComplementarityProblem *) malloc(sizeof(MixedLinearComplementarityProblem));

  //   char *arg1 = (char *) 0 ;
  //   int res1 ;
  //   char *buf1 = 0 ;
  //   int alloc1 = 0 ;

  //   res1 = SWIG_AsCharPtrAndSize(o1, &buf1, NULL, &alloc1);
  //   // if (!SWIG_IsOK(res1)) {
  //   //   SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "MixedLinearComplementarity_newFromFilename" "', argument " "1"" of type '" "char *""'");
  //   // }
  //   arg1 = reinterpret_cast< char * >(buf1);
  //   {
  //     try
  //     {
  //       result = (int)mixedLinearComplementarity_newFromFilename(MLCP,arg1);
  //     }
  //     catch (const std::invalid_argument& e)
  //     {
  //       // SWIG_exception(SWIG_ValueError, e.what());
  //     }
  //   }

  //   return MLCP;

  // }

};


%clear double* q;
%apply (double *qbck) { (double *q) };
