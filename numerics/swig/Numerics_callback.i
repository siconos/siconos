

%define %TARGET_MATRIX_FROM_CALL(mat, target_mat)
{
  switch(mat->storageType)
  {
  case NM_DENSE:
  {
    set_existing_dense_mat_from_target(mat->matrix0, (SN_ARRAY_TYPE*)target_mat, mat->size0, mat->size1, TARGET_ERROR_VERBOSE);
    break;
  }
  case NM_SPARSE_BLOCK:
  case NM_SPARSE:
  {
    // XXX Quite a bit of overhead here
    %NM_convert_from_target(target_mat, &(mat), TARGET_ERROR_VERBOSE);
    break;
  }
  default:
  {
    SWIG_Error(SWIG_TypeError, "TARGET_MATRIX_FROM_CALL :: unsupported storage type");
    TARGET_ERROR_VERBOSE;
  }
  }
}
%enddef

%inline %{

enum {ENV_IS_C_STRUCT = -1, ENV_IS_UNKNOWN, ENV_IS_PYTHON_CLASS, ENV_IS_PYTHON_FUNCTIONS, ENV_IS_PYTHON_FUNCTIONS_WITH_PROJ, ENV_IS_MATLAB_FUNCTION_HANDLES, ENV_IS_MATLAB_FUNCTION_NAMES};

#define CALL_COMPUTE_FN_NAME(FN_NAME, PNAME) FN_NAME ## PNAME
%}

#define CALL_COMPUTE_F(PNAME, ENV_FROM_PROBLEM) \
static void call_compute_F (void *problem, int n, double* z, double* F) \
  { \
    SN_OBJ_TYPE* target_z; \
    C_to_target_lang1(target_z, n, z, TARGET_ERROR_VERBOSE); \
 \
    SN_OBJ_TYPE* target_F; \
    TARGET_VECTOR_TO_CALL(F, target_F, n); \
 \
    SN_OBJ_TYPE* target_n = SWIG_From_int(n); \
 \
    void* env =  ENV_FROM_PROBLEM(problem); \
 \
    TARGET_CALL(env, "compute_F", env_compute_function, target_F, target_n, target_z); \
 \
    TARGET_VECTOR_FROM_CALL(F, (SN_ARRAY_TYPE*)target_F, n); \
 \
    target_mem_mgmt_instr(target_z); \
    target_mem_mgmt_instr(target_F); \
    target_mem_mgmt_instr(target_n); \
  };

#define CALL_COMPUTE_NABLA_F(PNAME, ENV_FROM_PROBLEM) \
  static void call_compute_nabla_F (void * problem, int n, double* z, NumericsMatrix* nabla_F) \
  { \
    SN_OBJ_TYPE* target_z; \
    C_to_target_lang1(target_z, n, z, TARGET_ERROR_VERBOSE); \
 \
    SN_OBJ_TYPE* target_nabla_F; \
    TARGET_MATRIX_TO_CALL(nabla_F, target_nabla_F); \
 \
    SN_OBJ_TYPE* target_n = SWIG_From_int(n); \
 \
    void* env = ENV_FROM_PROBLEM(problem); \
 \
    TARGET_CALL(env, "compute_nabla_F", env_compute_jacobian, target_nabla_F, target_n, target_z); \
 \
    %TARGET_MATRIX_FROM_CALL(nabla_F, target_nabla_F); \
 \
    target_mem_mgmt_instr(target_z); \
    target_mem_mgmt_instr(target_nabla_F); \
    target_mem_mgmt_instr(target_n); \
  };
