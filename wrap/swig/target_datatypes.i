#ifdef SWIGPYTHON

#define SN_ARRAY_TYPE PyArrayObject
#define SN_OBJ_TYPE PyObject
#define TARGET_ERROR_VERBOSE PyErr_PrintEx(0)
#define SN_SWIG_ERROR_CODE NULL
#define SN_ARRAY_INT_TYPE npy_intp

%inline %{

#define SN_ARRAY_TYPE PyArrayObject

#define C_to_target_lang1_detail(out, dim, cobj, FAILED_ACTION, CLASSID) \
  npy_intp pydims[] = { dim }; \
  out = FPyArray_SimpleNewFromData(1, pydims, CLASSID, cobj); \
  if (!out) { SWIG_Error(SWIG_RuntimeError, "Could not create an array from C object. Please file a bug"); FAILED_ACTION; }

#define C_to_target_lang1(out, dim, cobj, FAILED_ACTION) C_to_target_lang1_detail(out, dim, cobj, FAILED_ACTION, NPY_DOUBLE)
#define C_to_target_lang1_int(out, dim, cobj, FAILED_ACTION) C_to_target_lang1_detail(out, dim, cobj, FAILED_ACTION, NPY_INT)

#define C_to_target_lang1_alloc(out, outarray, dim, FAILED_ACTION) \
  npy_intp pydims[] = { dim }; \
  outarray = PyArray_SimpleNew(1, pydims, NPY_DOUBLE); \
  if (!outarray) { SWIG_Error(SWIG_RuntimeError, "Could not create an array"); FAILED_ACTION; } \
  out = (double*) array_data(outarray);

#define C_to_target_lang2_alloc(out, outarray, dim0, dim1, FAILED_ACTION) \
    npy_intp pydims[] = { dim0, dim1 }; \
    outarray = PyArray_SimpleNew(2, pydims, NPY_DOUBLE); \
    if (!outarray) { SWIG_Error(SWIG_RuntimeError, "Could not create an array"); FAILED_ACTION; } \
    out = (double*) array_data(outarray);

#define C_to_target_lang2(out, dim0, dim1, cobj, FAILED_ACTION) \
  npy_intp pydims[] = { dim0, dim1 }; \
  out = FPyArray_SimpleNewFromData(2, pydims, NPY_DOUBLE, cobj); \
  if (!out) { SWIG_Error(SWIG_RuntimeError, "Could not create an array from C object. Please file a bug"); FAILED_ACTION; }

#define sn_check_array_type(X, FAILED_ACTION) if (!X) { SWIG_Error(SWIG_RuntimeError, "Could not create an array from C object. Please file a bug"); FAILED_ACTION; }
#define sn_check_size_mat_vec(size, vec, FAILED_ACTION) if (array_numdims(vec) != 1 || size != array_size(vec, 0)) { SWIG_Error(SWIG_RuntimeError, "Matrix and vector given as arguments have different sizes"); FAILED_ACTION; }

#define set_vec_from_target(dest, src, CTRL_ACTION, FAILED_ACTION) \
{ \
  size_t sizeB = array_size(src, 0) * sizeof(double); \
  if (dest) { SWIG_Error(SWIG_RuntimeError, "destination vector is already allocated"); FAILED_ACTION; } \
  dest = (double*)malloc(sizeB); \
  memcpy(dest, array_data(src), sizeB); \
}

%}

#define target_mem_mgmt_instr(obj) Py_DECREF(obj)
#define target_mem_mgmtX_instr(obj) Py_XDECREF(obj)
#define target_mem_mgmt(cond, obj) if(cond && obj) { Py_DECREF(obj); }
#define target_mem_mgmtX(cond, obj) if(cond && obj) { Py_XDECREF(obj); }
#endif /* SWIGPYTHON */

#ifdef SWIGMATLAB
#define SN_ARRAY_TYPE mxArray
#define SN_OBJ_TYPE mxArray
#define TARGET_ERROR_VERBOSE
#define SN_SWIG_ERROR_CODE 1
#define SN_ARRAY_INT_TYPE int

// XXX really unsure  ???
#define SWIG_POINTER_NEW 1

%inline %{
// XXX HACK --xhub
#pragma GCC diagnostic ignored "-Wunused-variable"

#define SN_ARRAY_TYPE mxArray

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

#define set_vec_from_target(dest, src, CTRL_ACTION, FAILED_ACTION) \
{ \
  size_t sizeB = array_size(src, 0) > 1 ? array_size(src, 0) * sizeof(double) : array_size(src, 1) * sizeof(double); \
  if (dest) { SWIG_Error(SWIG_RuntimeError, "destination vector is already allocated"); FAILED_ACTION; } \
  dest = (double*)malloc(sizeB); \
  memcpy(dest, array_data(src), sizeB); \
}
%}

#define target_mem_mgmt_instr(obj)
#define target_mem_mgmtX_instr(obj)
#define target_mem_mgmt(cond, obj)
#define target_mem_mgmtX(cond, obj)

#endif /* SWIGMATLAB */

