#ifdef SWIGPYTHON

#define SN_ARRAY_TYPE PyArrayObject
#define SN_OBJ_TYPE PyObject
#define TARGET_ERROR_VERBOSE PyErr_PrintEx(0)
#define SN_SWIG_ERROR_CODE NULL
#define SN_ARRAY_INT_TYPE npy_intp

%inline %{

#define SN_ARRAY_TYPE PyArrayObject
#define SN_OBJ_TYPE PyObject
#define TARGET_ERROR_VERBOSE PyErr_PrintEx(0)

#define C_to_target_lang1_detail(out, dim, cobj, FAILED_ACTION, CLASSID) \
{ \
  npy_intp pydims[] = { dim }; \
  out = FPyArray_SimpleNewFromData(1, pydims, CLASSID, cobj); \
  if (!out) { SWIG_Error(SWIG_RuntimeError, "Could not create an array from C object. Please file a bug"); FAILED_ACTION; } \
}

#define C_to_target_lang1(out, dim, cobj, FAILED_ACTION) C_to_target_lang1_detail(out, dim, cobj, FAILED_ACTION, NPY_DOUBLE)
#define C_to_target_lang1_int(out, dim, cobj, FAILED_ACTION) C_to_target_lang1_detail(out, dim, cobj, FAILED_ACTION, NPY_INT)

#define C_to_target_lang1_alloc(out, outarray, dim, FAILED_ACTION) \
{ \
  npy_intp pydims[] = { dim }; \
  outarray = PyArray_SimpleNew(1, pydims, NPY_DOUBLE); \
  if (!outarray) { SWIG_Error(SWIG_RuntimeError, "Could not create an array"); FAILED_ACTION; } \
  out = (double*) array_data(outarray); \
}

#define C_to_target_lang2_alloc(out, outarray, dim0, dim1, FAILED_ACTION) \
{ \
  npy_intp pydims[] = { dim0, dim1 }; \
  outarray = PyArray_SimpleNew(2, pydims, NPY_DOUBLE); \
  if (!outarray) { SWIG_Error(SWIG_RuntimeError, "Could not create an array"); FAILED_ACTION; } \
  out = (double*) array_data(outarray); \
}

#define C_to_target_lang2(out, dim0, dim1, cobj, FAILED_ACTION) \
{ \
  npy_intp pydims[] = { dim0, dim1 }; \
  out = FPyArray_SimpleNewFromData(2, pydims, NPY_DOUBLE, cobj); \
  if (!out) { SWIG_Error(SWIG_RuntimeError, "Could not create an array from C object. Please file a bug"); FAILED_ACTION; } \
}

#define sn_check_array_type(X, FAILED_ACTION) if (!X) { SWIG_Error(SWIG_RuntimeError, "Could not create an array from C object. Please file a bug"); FAILED_ACTION; }
#define sn_check_size_mat_vec(size, vec, FAILED_ACTION) if ((array_numdims(vec) == 1 && size != array_size(vec, 0)) || (array_numdims(vec) == 2 && size != array_size(vec, 0)*array_size(vec, 1) && (array_size(vec, 0)==1 || array_size(vec, 1)==1))) { SWIG_Error(SWIG_RuntimeError, "Matrix and vector given as arguments have different sizes"); FAILED_ACTION; }
#define sn_check_size_mat(dim0, dim1, mat, FAILED_ACTION) if (!CHECK_ARRAY_SIZE(dim0, mat, 0) || !CHECK_ARRAY_SIZE(dim1, mat, 1)) { SWIG_Error(SWIG_RuntimeError, "Matrix does not have the right size"); FAILED_ACTION; }

#define set_vec_from_target(dest, src, CTRL_ACTION, FAILED_ACTION) \
{ \
  size_t sizeB = array_size(src, 0) * sizeof(double); \
  if (dest) { SWIG_Error(SWIG_RuntimeError, "destination vector is already allocated"); FAILED_ACTION; } \
  dest = (double*)malloc(sizeB); \
  memcpy(dest, array_data(src), sizeB); \
}

%}

%define check_save_target_fn(py_fn, env, env_field, callback, problem_field, arg_nb_str)
  if (PyCallable_Check(py_fn))
  {
    problem_field = callback;
    if (!env) env = malloc(sizeof(functions_env_python));
    functions_env_python* vi_env_python = (functions_env_python*) env;
    vi_env_python->id = ENV_IS_PYTHON_FUNCTIONS;
    vi_env_python->env_field = py_fn;
  }
  else
  {
    SWIG_Error(SWIG_TypeError, "argument " #arg_nb_str " must be callable");
    TARGET_ERROR_VERBOSE;
  }
%enddef

#define target_mem_mgmt_instr(obj) Py_DECREF(obj)
#define target_mem_mgmtX_instr(obj) Py_XDECREF(obj)
#define target_mem_mgmt(cond, obj) if(cond && obj) { Py_DECREF(obj); }
#define target_mem_mgmtX(cond, obj) if(cond && obj) { Py_XDECREF(obj); }

#define TARGET_CALL PY_CALL_METHOD_OR_FUNCTION

// We stole a reference in Python
#define TARGET_VECTOR_TO_CALL(vec, target_vec, n) C_to_target_lang1(target_vec, n, vec, TARGET_ERROR_VERBOSE);

#define TARGET_VECTOR_FROM_CALL(vec, target_vec, n) sn_check_size_mat_vec(n, target_vec, TARGET_ERROR_VERBOSE)

#define TARGET_MATRIX_TO_CALL(mat, target_mat) target_mat = NM_to_python(mat);

#endif /* SWIGPYTHON */

#ifdef SWIGMATLAB
#define SN_ARRAY_TYPE mxArray
#define SN_OBJ_TYPE mxArray
#define TARGET_ERROR_VERBOSE mexErrMsgIdAndTxt(SWIG_ErrorType(SWIG_lasterror_code), SWIG_lasterror_msg)
#define SN_SWIG_ERROR_CODE 1
#define SN_ARRAY_INT_TYPE int

// XXX really unsure  ???
#define SWIG_POINTER_NEW 1

%inline %{
// XXX HACK --xhub
#pragma GCC diagnostic ignored "-Wunused-variable"

#define SN_ARRAY_TYPE mxArray
#define SN_OBJ_TYPE mxArray
#define TARGET_ERROR_VERBOSE mexErrMsgIdAndTxt(SWIG_ErrorType(SWIG_lasterror_code), SWIG_lasterror_msg)

#define C_to_target_lang1_detail(out, dim, cobj, FAILED_ACTION, CLASSID, C_TYPE) \
  out = mxCreateUninitNumericMatrix(dim, 1, CLASSID, mxREAL);\
  if (!out) { SWIG_Error(SWIG_RuntimeError, "Could not create an array from C object. Check your available memory and file a bug"); FAILED_ACTION; }\
  memcpy(mxGetData(out), cobj, dim * sizeof(C_TYPE));

#define C_to_target_lang1(out, dim, cobj, FAILED_ACTION) C_to_target_lang1_detail(out, dim, cobj, FAILED_ACTION, mxDOUBLE_CLASS, double);
#define C_to_target_lang1_int(out, dim, cobj, FAILED_ACTION) C_to_target_lang1_detail(out, dim, cobj, FAILED_ACTION, sizeof(int) == 8 ? mxINT64_CLASS : mxINT32_CLASS, int);

#define C_to_target_lang1_alloc(out, outarray, dim, FAILED_ACTION) \
  outarray = mxCreateUninitNumericMatrix(dim, 1, mxDOUBLE_CLASS, mxREAL);\
  if (!out) { SWIG_Error(SWIG_RuntimeError, "Could not create an array from C object. Check your available memory"); FAILED_ACTION; }\
  out = (double*)  array_data((SN_ARRAY_TYPE*)outarray);

#define C_to_target_lang2_alloc(out, outarray, dim0, dim1, FAILED_ACTION) \
  outarray = mxCreateUninitNumericMatrix(dim0, dim1, mxDOUBLE_CLASS, mxREAL); \
  if (!out) { SWIG_Error(SWIG_RuntimeError, "Could not create an array from C object. Check your available memory"); FAILED_ACTION; }\
  out = (double*)  array_data((SN_ARRAY_TYPE*)outarray);

#define C_to_target_lang2(out, dim0, dim1, cobj, FAILED_ACTION) \
  out = mxCreateUninitNumericMatrix(dim0, dim1, mxDOUBLE_CLASS, mxREAL);\
  if (!out) { SWIG_Error(SWIG_RuntimeError, "Could not create an array from C object. Check your available memory and file a bug"); FAILED_ACTION; }\
  memcpy(mxGetData(out), cobj, dim0 * dim1 * sizeof(double));

#define sn_check_array_type(X, FAILED_ACTION) if (!X) { SWIG_Error(SWIG_RuntimeError, "Could not create an array from C object. Please file a bug"); FAILED_ACTION; }
#define sn_check_size_mat_vec(size, vec, FAILED_ACTION) if (mxGetNumberOfDimensions(vec) != 2 || !((array_size(vec, 0) == 1 && (array_size(vec, 1) == size)) || (array_size(vec, 0) == size && (array_size(vec, 1) == 1)))) { SWIG_Error(SWIG_RuntimeError, "Matrix and vector given as arguments have different sizes"); FAILED_ACTION; }
#define sn_check_size_mat(dim0, dim1, mat, FAILED_ACTION) if (!CHECK_ARRAY_SIZE(dim0, mat, 0) || !CHECK_ARRAY_SIZE(dim1, mat, 1)) { SWIG_Error(SWIG_RuntimeError, "Matrix does not have the right size"); FAILED_ACTION; }


#define set_vec_from_target(dest, src, CTRL_ACTION, FAILED_ACTION) \
{ \
  size_t sizeB = array_size(src, 0) > 1 ? array_size(src, 0) * sizeof(double) : array_size(src, 1) * sizeof(double); \
  if (dest) { SWIG_Error(SWIG_RuntimeError, "destination vector is already allocated"); FAILED_ACTION; } \
  dest = (double*)malloc(sizeB); \
  memcpy(dest, array_data(src), sizeB); \
}

#define STR_CONCAT(A, B) A ## _ ## B

#define TARGET_CALL(env, fn_str, field_fn, output, ...) \
{ \
  functions_env_matlab* env_matlab = (functions_env_matlab*)env; \
  switch (env_matlab->id) \
  { \
  case ENV_IS_MATLAB_FUNCTION_HANDLES: \
  { \
    SN_OBJ_TYPE* plhs[] = {NULL};\
    SN_OBJ_TYPE* _arr_[] = {env_matlab->field_fn, __VA_ARGS__}; \
    int nlhs = 1; \
    int nrhs = sizeof(_arr_)/sizeof(SN_OBJ_TYPE*); \
\
    mexCallMATLAB(nlhs, plhs, nrhs, _arr_, "feval");\
\
    if (!plhs[0])\
    {\
      SWIG_Error(SWIG_RuntimeError, "matlab function should return one array.");\
      TARGET_ERROR_VERBOSE;\
    }\
    output = plhs[0]; \
    break;\
  }\
  case ENV_IS_MATLAB_FUNCTION_NAMES:\
  {\
    SN_OBJ_TYPE* plhs[] = {NULL};\
    SN_OBJ_TYPE* _arr_[] = {NULL, __VA_ARGS__}; \
    int nlhs = 1; \
    char* vstr = STR_CONCAT(env_matlab->field_fn,str); \
    if (vstr) { \
      if (!env_matlab->field_fn) { env_matlab->field_fn = mexGetVariablePtr("base", vstr); } \
      if (!env_matlab->field_fn) { env_matlab->field_fn = mexGetVariablePtr("caller", vstr); } \
      if (!env_matlab->field_fn) { env_matlab->field_fn = mexGetVariablePtr("global", vstr); } \
      if (!env_matlab->field_fn) { SWIG_Error(SWIG_RuntimeError, "null pointer to %s", vstr); TARGET_ERROR_VERBOSE; } \
      _arr_[0] =  env_matlab->field_fn; \
      mexCallMATLAB(nlhs, plhs, sizeof(_arr_)/sizeof(SN_OBJ_TYPE*), _arr_, "feval");\
      output = plhs[0]; \
      break;\
    } \
    int nrhs = sizeof(_arr_)/sizeof(SN_OBJ_TYPE*)-1; \
\
    mexCallMATLAB(nlhs, plhs, nrhs, nrhs > 0 ? &_arr_[1] : (SN_OBJ_TYPE**)NULL, (char*)env_matlab->field_fn);\
\
    if (!plhs[0])\
    {\
      SWIG_Error(SWIG_RuntimeError, format_msg_concat("matlab function named should return 1 array. Faulty function: ", (char*)env_matlab->field_fn));\
      TARGET_ERROR_VERBOSE;\
    }\
    output = plhs[0]; \
    break; \
  }\
  default:\
  {\
    SWIG_Error(SWIG_RuntimeError, "unsupported env type %d in TARGET_CALL for MATLAB. File a bug report", env_matlab->id);\
    TARGET_ERROR_VERBOSE;\
  }\
  } \
}

%}

%define check_save_target_fn(matlab_fn, env, env_field, callback, problem_field, arg_nb_str)
{
  switch(mxGetClassID(matlab_fn))
  {
  case mxFUNCTION_CLASS:
  {
    SWIG_Error(SWIG_RuntimeError, "Unsupported function handle argument");
    TARGET_ERROR_VERBOSE;

    problem_field = callback;
    if (!env) env = malloc(sizeof(functions_env_matlab));
    functions_env_matlab* vi_env_matlab = (functions_env_matlab*) env;
    vi_env_matlab->id = ENV_IS_MATLAB_FUNCTION_HANDLES;
    vi_env_matlab->env_field = (void*)matlab_fn;
    break;
  }
  case mxCHAR_CLASS:
  {
    SN_OBJ_TYPE* plhs[] = {NULL};
    SN_OBJ_TYPE* prhs[] = {matlab_fn};
    if (mexCallMATLAB(1, plhs, 1, prhs, "exist"))
    {
      SWIG_Error(SWIG_RuntimeError, "Call to exists failed");
      TARGET_ERROR_VERBOSE;
    }
    size_t m = mxGetM(matlab_fn);
    size_t n = mxGetN(matlab_fn);
    size_t len_str = (m > n ? m : n) + 1;
    if (!env) env = malloc(sizeof(functions_env_matlab));
    functions_env_matlab* vi_env_matlab = (functions_env_matlab*) env;
    vi_env_matlab->env_field = (char*) malloc(len_str * sizeof(char));
    vi_env_matlab->STR_CONCAT(env_field,str) = NULL;
    if (mxGetString(matlab_fn, vi_env_matlab->env_field, len_str))
    {
      free(vi_env_matlab->env_field);
      vi_env_matlab->env_field = NULL;
      SWIG_Error(SWIG_RuntimeError, "Impossible to get a C string representation");
      TARGET_ERROR_VERBOSE;
    }
    if (mxGetScalar(plhs[0]) == 0)
    {
     vi_env_matlab->STR_CONCAT(env_field,str) = vi_env_matlab->env_field;
     vi_env_matlab->env_field = NULL;
    }
    mxDestroyArray(plhs[0]);
    problem_field = callback;
    vi_env_matlab->id = ENV_IS_MATLAB_FUNCTION_NAMES;
    break;
  }
  default:
  {
    SWIG_Error(SWIG_TypeError, "argument " #arg_nb_str " has to be a function name (string) or a handle to a function");
    TARGET_ERROR_VERBOSE;
  }
  }
}
%enddef

#define target_mem_mgmt_instr(obj) mxDestroyArray(obj)
#define target_mem_mgmtX_instr(obj) if(obj) { mxDestroyArray (obj); }
#define target_mem_mgmt(cond, obj)
#define target_mem_mgmtX(cond, obj)

#define TARGET_VECTOR_TO_CALL(vec, target_vec, n)

#define TARGET_VECTOR_FROM_CALL(vec, target_vec, n) set_existing_vec_from_target(vec, target_vec, n, , TARGET_ERROR_VERBOSE)

#define TARGET_MATRIX_TO_CALL(mat, target_mat)

#endif /* SWIGMATLAB */

//common stuff

#define set_existing_vec_from_target(dest, src, size, CTRL_ACTION, FAILED_ACTION) \
{ \
  if (CHECK_ARRAY_VECTOR((SN_ARRAY_TYPE*)src)) { SWIG_Error(SWIG_RuntimeError, "object is not a vector"); FAILED_ACTION; } \
  sn_check_size_mat_vec(size, (SN_ARRAY_TYPE*)src, FAILED_ACTION); \
  if (!dest) { SWIG_Error(SWIG_RuntimeError, "destination vector is not allocated"); FAILED_ACTION; } \
  memcpy(dest, array_data((SN_ARRAY_TYPE*)src), size * sizeof(double)); \
}

#define set_existing_dense_mat_from_target(dest, src, dim0, dim1, FAILED_ACTION) \
{ \
  if (CHECK_ARRAY_MATRIX((SN_ARRAY_TYPE*)src)) { SWIG_Error(SWIG_RuntimeError, "object is not a matrix"); FAILED_ACTION; } \
  sn_check_size_mat(dim0, dim1, (SN_ARRAY_TYPE*)src, FAILED_ACTION); \
  if (!dest) { SWIG_Error(SWIG_RuntimeError, "destination (dense) matrix is not allocated"); FAILED_ACTION; } \
  memcpy(dest, array_data((SN_ARRAY_TYPE*)src), dim0 * dim1 * sizeof(double)); \
}
