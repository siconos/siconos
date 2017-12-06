//#define GET_INT(OBJ,VAR)                                                \
//  if(!PyInt_Check(OBJ))                                                 \
//  {                                                                     \
//    throw(std::invalid_argument(#OBJ ": expecting an int")); \
//  }                                                                     \
//  VAR = PyInt_AsLong(OBJ)

// need a PySequence_Check
#define GET_INTS(OBJ,INDEX,VAR)                                          \
  PyObject *_TEMP##VAR = PySequence_GetItem(OBJ,INDEX);                 \
  if (!PyInt_Check(_TEMP##VAR))                                         \
  {                                                                     \
    Py_XDECREF(_TEMP##VAR);                                             \
    PyErr_SetString(PyExc_RuntimeError, "expecting an int for " # VAR);       \
    PyObject_Print(OBJ, stderr, 0);                                     \
    return 0;                                                           \
  }                                                                     \
  VAR = PyInt_AsLong(_TEMP##VAR);                                   \
  Py_DECREF(_TEMP##VAR)

%inline %{

#include "SparseMatrix.h"

static inline bool sane_pyseq_check(PyObject *o)
{
  if (PySequence_Check(o) && (PyList_Check(o) || PyTuple_Check(o)))
  {
    PyObject* elt = PySequence_GetItem(o, 0);
    if (elt && (PyInt_Check(elt) || PyFloat_Check(elt)))
    {
      Py_DECREF(elt);
      return true;
    }
    Py_XDECREF(elt);
  }
  return false;
}
#define CHECK_PYINT(PYSEQ, INDX, OUT) \
  PyObject *_TEMP##OUT = PySequence_GetItem(PYSEQ, INDX);                 \
  if (!PyInt_Check(_TEMP##OUT))                                         \
  {                                                                     \
    Py_XDECREF(_TEMP##OUT);                                             \
    PyErr_SetString(PyExc_RuntimeError, "expecting an int for " # OUT);       \
    PyObject_Print(PYSEQ, stderr, 0);                                     \
  }                                                                     \
  OUT = PyInt_AsLong(_TEMP##OUT);                                   \
  Py_DECREF(_TEMP##OUT)


#ifndef NDEBUG
static inline void _sn_check_nnz(PyObject** mat, CSparseMatrix *M)
{
  if (!*mat) { return; }
  PyObject *auto_nnz = PyObject_GetAttrString(*mat, "nnz");
  CS_INT nz;
  if (M->nz >= 0) { nz = M->nz; } else { nz = M->nzmax; }
  if(PyInt_AsLong(auto_nnz) != nz) {  PyErr_SetString(PyExc_RuntimeError, "number of nnz is inconsistent"); *mat = NULL; }
  Py_XDECREF(auto_nnz);
}

static inline void _sn_check_shape(PyObject** mat, CSparseMatrix *M)
{
  if (!*mat) { return; }
  PyObject *shape_ = PyObject_GetAttrString(*mat, "shape");
  unsigned nrows, ncols;
  CHECK_PYINT(shape_, 0, nrows);
  CHECK_PYINT(shape_, 1, ncols);

  Py_DECREF(shape_);

  if (nrows != M->m) { PyErr_SetString(PyExc_RuntimeError, "number of rows is inconsistent"); *mat = NULL;}
  if (ncols != M->n) { PyErr_SetString(PyExc_RuntimeError, "number of columns is inconsistent"); *mat = NULL;}
}
#else
static inline void _sn_check_nnz(PyObject** mat, CSparseMatrix *M) {};
static inline void _sn_check_shape(PyObject** mat, CSparseMatrix *M) {};
#endif

#include "SiconosConfig.h"

// Work-around for "disappearing" SICONOS_INT64 problem in use of
// %#ifdef below
#ifdef SICONOS_INT64
#define _SICONOS_INT64
#endif

#if defined(SICONOS_INT64) && !defined(SICONOS_FORCE_NPY_INT32)
#define NPY_INT_TYPE NPY_INT64
#else
#define NPY_INT_TYPE NPY_INT32
#endif

#if defined(SICONOS_FORCE_NPY_INT32) && defined(SICONOS_INT64)

#define INT_TO_NPY_INT(dim, intp, out, copy) \
  { \
  int32_t* int32p = (int32_t*) malloc(dim[0] * sizeof(int32_t)); \
  if (!int32p) {  PyErr_SetString(PyExc_RuntimeError, "Could not allocate memory to convert " # intp "to 32 bits"); return NULL; }; \
  for (size_t i = 0; i < (size_t)dim[0]; ++i) { int32p[i] = intp[i]; }        \
  out  = PyArray_SimpleNewFromData(1, dim, NPY_INT_TYPE, int32p); \
  if(!out) {  PyErr_SetString(PyExc_RuntimeError, "Could not extract " # intp); return NULL; }; \
  PyArray_UpdateFlags((PyArrayObject*)out, NPY_ARRAY_OWNDATA); \
  }

#else

#define INT_TO_NPY_INT(dim, intp, out, copy) \
  { \
  CS_INT * int_p; \
  if (copy) \
  { \
    int_p = (CS_INT*)malloc(dim[0] * sizeof(CS_INT)); \
    memcpy(int_p, intp, dim[0] * sizeof(CS_INT)); \
  } \
  else \
  { \
    int_p = intp; \
  } \
  out  = PyArray_SimpleNewFromData(1, dim, NPY_INT_TYPE, int_p); \
  if(!out) {  PyErr_SetString(PyExc_RuntimeError, "Could not extract " # intp); return NULL; }; \
  if (copy) { PyArray_UpdateFlags((PyArrayObject*)out, NPY_ARRAY_OWNDATA); } \
  } \

#endif

#define CS_TO_SCIPY(TYPE, P_LEN, copy) \
  if (!M) \
  { \
    Py_INCREF(Py_None); \
    return Py_None; \
  } \
  else \
  { \
    /* get sys.modules dict */ \
    PyObject* sys_mod_dict = PyImport_GetModuleDict(); \
\
    /* get the csr module object */ \
    PyObject* scipy_mod = PyMapping_GetItemString(sys_mod_dict, (char *)"scipy.sparse." #TYPE);\
\
    if (!scipy_mod) \
    { \
      PyErr_SetString(PyExc_RuntimeError, "Did you import scipy.sparse." #TYPE "?"); \
      return NULL; \
    } \
\
    npy_intp this_M_x_dims[1] = { M->nzmax }; \
    npy_intp this_M_i_dims[1] = { M->nzmax }; \
    npy_intp this_M_p_dims[1] = { P_LEN + 1 }; \
\
    double* data_p; \
    if (copy) \
    { \
      data_p = (double*)malloc(M->nzmax * sizeof(double)); \
      memcpy(data_p, M->x, M->nzmax * sizeof(double)); \
    } \
    else \
    { \
      data_p = M->x; \
    } \
\
    PyObject* out_data = PyArray_SimpleNewFromData(1, this_M_x_dims, NPY_DOUBLE, data_p); \
    if(!out_data) { PyErr_SetString(PyExc_RuntimeError, "Could not extract M->x"); return NULL; }; \
    if (copy) { PyArray_UpdateFlags((PyArrayObject*)out_data, NPY_ARRAY_OWNDATA); } \
\
    PyObject* out_indices; \
    INT_TO_NPY_INT(this_M_i_dims, M->i, out_indices, copy); \
\
    PyObject* out_indptr ; \
    INT_TO_NPY_INT(this_M_p_dims, M->p, out_indptr, copy); \
\
    /* Warning ! m is the number of rows, n the number of columns ! --xhub */ \
    PyObject* out_shape = PyTuple_Pack(2, PyInt_FromLong(M->m), PyInt_FromLong(M->n)); \
    if(!out_shape) {  PyErr_SetString(PyExc_RuntimeError, "Could not extract M->m or M->n"); return NULL; }; \
\
    PyObject* out_nnz = PyInt_FromLong(M->nzmax); \
    if(!out_nnz) {  PyErr_SetString(PyExc_RuntimeError, "Could not extract M->nzmax"); return NULL; }; \
\
    /* call the class inside the csr module */ \
    PyObject* out_mat = PyObject_CallMethodObjArgs(scipy_mod, PyString_FromString((char *) #TYPE "_matrix"), out_shape, NULL); \
\
    if(out_mat) \
    { \
      PyObject_SetAttrString(out_mat,"data", out_data); \
      PyObject_SetAttrString(out_mat,"indices", out_indices); \
      PyObject_SetAttrString(out_mat,"indptr", out_indptr); \
\
      _sn_check_nnz(&out_mat, M); \
      _sn_check_shape(&out_mat, M); \
\
      return out_mat; \
    } \
    else \
    { \
      PyErr_SetString(PyExc_RuntimeError, "Could not create " #TYPE " matrix"); \
      return NULL; \
    } \
  }

static inline bool is_Pyobject_scipy_sparse_matrix(PyObject* o, PyObject* scipy_mod)
 {
    bool ret;
    PyObject* res = PyObject_CallMethodObjArgs(scipy_mod, PyString_FromString("issparse"), o, NULL);

    if (!res) return false;

    ret = (res == Py_True);
    Py_DECREF(res);

    return ret;
}

%}

%define %SAFE_CAST_INT(pyvar, len, dest_array, array_pyvar, indvar, alloc)
{
    int array_pyvartype_ = PyArray_TYPE((PyArrayObject *)pyvar);
    switch (array_pyvartype_)
    {
      case NPY_INT32:
      {
        array_pyvar = obj_to_array_allow_conversion(pyvar, NPY_INT32, indvar);
        if (!array_pyvar) { PyErr_SetString(PyExc_RuntimeError, "Could not get array for variable" #pyvar); PyObject_Print(pyvar, stderr, 0); return 0; }

%#ifdef _SICONOS_INT64
        PyErr_Warn(PyExc_UserWarning, "Performance warning: the vector of indices or pointers is in int32, but siconos has 64-bits integers: we have to perform a conversion. Consider given sparse matrix in the right format");
        dest_array = (int64_t*) malloc(len * sizeof(int64_t));
        if(!dest_array) { PyErr_SetString(PyExc_RuntimeError, "Allocation of i or p failed (triggered by conversion to int32)"); return 0; }
        for(unsigned i = 0; i < len; ++i)
        {
          dest_array[i] = ((int32_t *) array_data(array_pyvar)) [i];
        }
        if (*indvar) Py_DECREF(array_pyvar);
        *indvar = 0;
        alloc = true;
%#else
        dest_array = (int32_t *) array_data(array_pyvar);
%#endif
        break;
      }
      case NPY_INT64:
      {
        array_pyvar = obj_to_array_allow_conversion(pyvar, NPY_INT64, indvar);
        if (!array_pyvar) { PyErr_SetString(PyExc_RuntimeError, "Could not get array for variable " #pyvar);  PyObject_Print(pyvar, stderr, 0); return 0; }

%#ifdef _SICONOS_INT64
        dest_array = (int64_t *) array_data(array_pyvar);
%#else
        PyErr_Warn(PyExc_UserWarning, "Performance warning: the vector of indices or pointers is in int64, but siconos has 32-bits integers: we have to perform a conversion. Consider given sparse matrix in the right format");
        dest_array = (int32_t*) malloc(len * sizeof(int32_t));
        if(!dest_array) { PyErr_SetString(PyExc_RuntimeError, "Allocation of i or p failed (triggered by conversion to int64)"); return 0; }
        for(unsigned i = 0; i < len; ++i)
        {
          dest_array[i] = ((int64_t *) array_data(array_pyvar)) [i];
        }
        if (*indvar) Py_DECREF(array_pyvar);
        *indvar = 0;
        alloc = true;
%#endif
        break;
      }
      default:
      {
        PyObject *errmsg;
        errmsg = PyUString_FromString("Unknown type ");
        PyUString_ConcatAndDel(&errmsg, PyObject_Repr((PyObject *)PyArray_DESCR((PyArrayObject *)pyvar)));
        PyUString_ConcatAndDel(&errmsg, PyUString_FromFormat(" for variable " #pyvar));
        PyErr_SetObject(PyExc_TypeError, errmsg);
        Py_DECREF(errmsg);
        return 0;
      }
    }
}
%enddef


%fragment("NumericsMatrix", "header", fragment="NumPy_Fragments")
{
  static int cs_convert_from_scipy_sparse(PyObject* obj, CSparseMatrix** m, PyArrayObject** array_data_, int* array_data_ctrl_, PyArrayObject** array_i_, int* array_i_ctrl_, PyArrayObject** array_p_, int* array_p_ctrl_, int* alloc_ctrl)
  {

  assert(m);
  /* get sys.modules dict */
  PyObject* sys_mod_dict = PyImport_GetModuleDict();
  /* get the scipy module object */ 
  PyObject* scipy_mod = PyMapping_GetItemString(sys_mod_dict, (char *)"scipy.sparse");

  if (!scipy_mod) 
  { 
    PyErr_SetString(PyExc_RuntimeError, "Did you import scipy.sparse ?");
    return 0;
  }

  bool isspmat = is_Pyobject_scipy_sparse_matrix(obj, scipy_mod);
  if (isspmat == false)
  {
    return -1;
  }
  else
  {
    PyObject* shape_ = PyObject_GetAttrString(obj, "shape");

    unsigned nrows, ncols;
    GET_INTS(shape_, 0, nrows);
    GET_INTS(shape_, 1, ncols);

    Py_DECREF(shape_);


    PyObject* res;

    res = PyObject_CallMethodObjArgs(scipy_mod, PyString_FromString("isspmatrix_csc"), obj, NULL);
    bool is_csc = (res == Py_True);
    Py_DECREF(res);

    if (is_csc)
    {
      // csc

      PyObject* nnz_ = PyObject_GetAttrString(obj, "nnz");
      size_t nzmax = PyInt_AsLong(nnz_);
      Py_DECREF(nnz_);

      CSparseMatrix* M = (CSparseMatrix*) calloc(1, sizeof(CSparseMatrix));
      if(!M) { PyErr_SetString(PyExc_RuntimeError, "Failed to allocate a cs_sparse"); return 0; }

      M->nz = -1;
      M->m = nrows;
      M->n = ncols;
      M->nzmax = nzmax;
      *m = M;
      if(!M) { PyErr_SetString(PyExc_RuntimeError, "Allocation of the csc matrix failed"); return 0; };

      PyObject* data_ = PyObject_GetAttrString(obj, "data");
      PyObject* indices_ = PyObject_GetAttrString(obj, "indices");
      PyObject* indptr_ = PyObject_GetAttrString(obj, "indptr");

      *array_data_ = obj_to_array_allow_conversion(data_, NPY_DOUBLE, array_data_ctrl_);
      if (!*array_data_) { PyErr_SetString(PyExc_RuntimeError, "Could not get a pointer to the data array");  PyObject_Print(data_, stderr, 0); return 0; }

      M->x = (double*)array_data(*array_data_);

      bool alloc_p = false;
      %SAFE_CAST_INT(indptr_, (M->n + 1), M->p, *array_p_, array_p_ctrl_, alloc_p);
      if (alloc_p) { *alloc_ctrl |= ALLOC_CTRL_P; };
      bool alloc_i = false;
      %SAFE_CAST_INT(indices_, nzmax, M->i, *array_i_, array_i_ctrl_, alloc_i);
      if (alloc_i) { *alloc_ctrl |= ALLOC_CTRL_I; };

      return 1;
    }

#ifdef WITH_CSR
    res = PyObject_CallMethodObjArgs(scipy_mod, PyString_FromString("isspmatrix_csr"), obj, NULL);
    bool is_csr = (res == Py_True);
    Py_DECREF(res);

    if (is_csr)
    {
      // csr

      PyObject* nnz_ = PyObject_GetAttrString(obj, "nnz");
      size_t nzmax = PyInt_AsLong(nnz_);
      Py_DECREF(nnz_);

      CSparseMatrix* M = (CSparseMatrix*) calloc(1, sizeof(CSparseMatrix));
      if(!M) { PyErr_SetString(PyExc_RuntimeError, "Failed to allocate a cs_sparse"); return 0; }

      M->nz = -2;
      M->nzmax = nzmax;
      M->m = nrows;
      M->n = ncols;
      *m = M;

      if(!M) { PyErr_SetString(PyExc_RuntimeError, "Allocation of the M matrix failed"); return 0; };

      PyObject* data_ = PyObject_GetAttrString(obj, "data");
      PyObject* indices_ = PyObject_GetAttrString(obj, "indices");
      PyObject* indptr_ = PyObject_GetAttrString(obj, "indptr");

      *array_data_ = obj_to_array_allow_conversion(data_, NPY_DOUBLE, array_data_ctrl_);
      if (!*array_data_) { PyErr_SetString(PyExc_RuntimeError, "Could not get a pointer to the data array");  PyObject_Print(data_, stderr, 0); return 0; }

      M->x = (double*)array_data(*array_data_);

      bool alloc_p = false;
      %SAFE_CAST_INT(indptr_, (nrows + 1), M->p, *array_p_, array_p_ctrl_, alloc_p);
      if (alloc_p) { *alloc_ctrl |= ALLOC_CTRL_P; };
      bool alloc_i = false;
      %SAFE_CAST_INT(indices_, nzmax, M->i, *array_i_, array_i_ctrl_, alloc_i);
      if (alloc_i) { *alloc_ctrl |= ALLOC_CTRL_I; };

      return 1;
    }
#endif /* WITH_CSR */

    res = PyObject_CallMethodObjArgs(scipy_mod, PyString_FromString("isspmatrix_coo"), obj, NULL);
    bool is_coo = (res == Py_True);
    Py_DECREF(res);

    PyObject* coo;
    int coo_new_alloc;
    if (!is_coo)
    {
      PyErr_Warn(PyExc_UserWarning, "Performance warning: the given sparse matrix is neither csc or coo, we have to perform a conversion to coo");
      coo = PyObject_CallMethodObjArgs(scipy_mod, PyString_FromString("coo_matrix"), obj, NULL);
      if (!coo) { if (!PyErr_Occurred()) { PyErr_SetString(PyExc_RuntimeError, "Conversion to coo failed!"); }; return 0; }
      coo_new_alloc = 1;
    }
    else
    {
      coo = obj;
      coo_new_alloc = 0;
    }

    // triplet
    PyObject* nnz_ = PyObject_GetAttrString(coo, "nnz");
    size_t nnz = PyInt_AsLong(nnz_);
    Py_DECREF(nnz_);


    CSparseMatrix* M = (CSparseMatrix*) calloc(1, sizeof(CSparseMatrix));
    if(!M) { PyErr_SetString(PyExc_RuntimeError, "Failed to allocate a cs_sparse"); return 0; }

    M->m = nrows;
    M->n = ncols;
    M->nzmax = nnz;
    *m = M;
    M->nz = nnz;

    if(!M) { PyErr_SetString(PyExc_RuntimeError, "Allocation of the triplet matrix failed"); return 0; }

    PyObject* data_ = PyObject_GetAttrString(coo, "data");
    PyObject* row_ = PyObject_GetAttrString(coo, "row");
    PyObject* col_ = PyObject_GetAttrString(coo, "col");

    *array_data_ = obj_to_array_allow_conversion(data_, NPY_DOUBLE, array_data_ctrl_);
    if (!*array_data_) { PyErr_SetString(PyExc_RuntimeError, "Could not get a pointer to the data array");  PyObject_Print(data_, stderr, 0); return 0; }

    M->x = (double*)array_data(*array_data_);

    bool alloc_p = false;
    %SAFE_CAST_INT(col_, nnz, M->p, *array_p_, array_p_ctrl_, alloc_p);
    if (alloc_p) { *alloc_ctrl |= ALLOC_CTRL_P; };
    bool alloc_i = false;
    %SAFE_CAST_INT(row_, nnz, M->i, *array_i_, array_i_ctrl_, alloc_i);
    if (alloc_i) { *alloc_ctrl |= ALLOC_CTRL_I; };

    if (coo_new_alloc)
    {
      Py_DECREF(coo);
    }

    return 1;
  }
  }

  static int NM_convert_from_scipy_sparse(PyObject* obj, NumericsMatrix* m, PyArrayObject** array_data_, int* array_data_ctrl_, PyArrayObject** array_i_, int* array_i_ctrl_, PyArrayObject** array_p_, int* array_p_ctrl_, int* alloc_ctrl)
  {
    CSparseMatrix* csm = NULL;
    int res = cs_convert_from_scipy_sparse(obj, &csm, array_data_, array_data_ctrl_, array_i_, array_i_ctrl_, array_p_, array_p_ctrl_, alloc_ctrl);
    if (res > 0)
    {
      m->storageType = NM_SPARSE;
      m->matrix2 = newNumericsSparseMatrix();

      if (csm->nz > 0)
      {
        m->matrix2->triplet = csm;
        m->matrix2->origin = NS_TRIPLET;
      }
      else if (csm->nz == -1)
      {
        m->matrix2->csc = csm;
        m->matrix2->origin = NS_CSC;
      }
      else if (csm->nz == -2)
      {
        m->matrix2->csr = csm;
        m->matrix2->origin = NS_CSR;
      }
      else
      {
        PyErr_SetString(PyExc_RuntimeError, "Unknown CSparseMatrix from cs_convert_from_scipy_sparse");
        return 0;
      }

      NM_update_size(m);
    }

    return res;
  }


  static NumericsMatrix* NM_convert_from_python(PyObject* obj, NumericsMatrix** tmpmat, PyArrayObject** array_data_, int* array_ctrl, PyArrayObject** array_i_, int* array_i_ctrl_, PyArrayObject** array_p_, int* array_p_ctrl_, int* alloc_ctrl)
  {
  void* argp = NULL;
  NumericsMatrix* out = NULL;
  int res = SWIG_ConvertPtr(obj, &argp, $descriptor(NumericsMatrix *), %convertptr_flags);
  if (SWIG_IsOK(res))
  {
    out = (NumericsMatrix *)argp;
  }
  else
  {
    *tmpmat = NM_new();
    out = *tmpmat;
    if (is_array(obj) || sane_pyseq_check(obj))
    {
      PyArrayObject* array_data = obj_to_array_fortran_allow_conversion(obj, NPY_DOUBLE, array_ctrl);

      if (!array_data)
      {
        PyErr_SetString(PyExc_TypeError, "Could not get array obj from the python object");
        PyObject_Print(obj, stderr, 0);
        goto fail;
      }

      if (!require_dimensions(array_data, 2) || !require_native(array_data) || !require_fortran(array_data))
      {
        PyErr_SetString(PyExc_TypeError, "The given object does not have the right structure. We expect a 2 dimensional array (or list, tuple, ...)");
        PyObject_Print(obj, stderr, 0);
        goto fail;
      }

      out->storageType = NM_DENSE;
      out->size0 =  array_size(array_data, 0);
      out->size1 =  array_size(array_data, 1);
      out->matrix0 = (double *)array_data(array_data);

      *array_data_ = array_data;
    }
    else
    {
      int sp_conv = NM_convert_from_scipy_sparse(obj, out, array_data_, array_ctrl, array_i_, array_i_ctrl_, array_p_, array_p_ctrl_, alloc_ctrl);
      if (!sp_conv) { goto fail; }
      else if (sp_conv < 0)
      {
        if (SWIG_IsOK(SWIG_ConvertPtr(obj, &argp, $descriptor(SparseBlockStructuredMatrix *), %convertptr_flags)))
        {
          out->matrix1 = (SparseBlockStructuredMatrix *)argp;
          out->storageType = NM_SPARSE_BLOCK;
          NM_update_size(out);
        }
        else
        {
          PyObject_Print(obj, stderr, 0);
          PyErr_SetString(PyExc_TypeError, "Cannot build a NumericsMatrix from the given python object");
          goto fail;
        }
      }
    }
  }

  return out;

fail:
  if (*tmpmat) { free(*tmpmat); *tmpmat = NULL; }
  return NULL;
  }

static PyObject* cs_sparse_to_csr_matrix(CSparseMatrix *M, bool copy)
{
  CS_TO_SCIPY(csr, M->m, copy);
}

static PyObject* cs_sparse_to_csc_matrix(CSparseMatrix *M, bool copy)
{
  CS_TO_SCIPY(csc, M->n, copy);
}

static PyObject* cs_sparse_to_coo_matrix(CSparseMatrix *M, bool copy)
{
  if (!M)
  {
    Py_INCREF(Py_None);
    return Py_None;
  }
  else
  {
    /* get sys.modules dict */
    PyObject* sys_mod_dict = PyImport_GetModuleDict();

    /* get the csr module object */
    PyObject* scipy_mod = PyMapping_GetItemString(sys_mod_dict, (char *)"scipy.sparse.coo");\

    if (!scipy_mod)
    {
      PyErr_SetString(PyExc_RuntimeError, "Did you import scipy.sparse.coo?");
      return NULL;
    }

    npy_intp this_M_x_dims[1] = { M->nz };
    npy_intp this_M_i_dims[1] = { M->nz };
    npy_intp this_M_p_dims[1] = { M->nz };

    double* data_p;
    if (copy)
    {
      data_p = (double*)malloc(M->nz * sizeof(double));
      memcpy(data_p, M->x, M->nz * sizeof(double));
    }
    else
    {
      data_p = M->x;
    }

    PyObject* out_data = PyArray_SimpleNewFromData(1, this_M_x_dims, NPY_DOUBLE, data_p);
    if(!out_data) { PyErr_SetString(PyExc_RuntimeError, "Could not extract M->x"); return NULL; };
    if (copy) { PyArray_UpdateFlags((PyArrayObject*)out_data, NPY_ARRAY_OWNDATA); }

    PyObject* row_indices;
    PyObject* col_indices;

    INT_TO_NPY_INT(this_M_i_dims, M->i, row_indices, copy);
    INT_TO_NPY_INT(this_M_p_dims, M->p, col_indices, copy);

    PyObject* out_indx = PyTuple_Pack(2, row_indices, col_indices);
    if(!out_indx) { PyErr_SetString(PyExc_RuntimeError, "Could not build (row, col)"); return NULL; };
    PyObject* out_all =  PyTuple_Pack(2, out_data, out_indx);
    if(!out_all) { PyErr_SetString(PyExc_RuntimeError, "Could not build (data, (row, col))"); return NULL; };

    /* Warning ! m is the number of rows, n the number of columns ! --xhub */ \
    PyObject* out_shape = PyTuple_Pack(2, PyInt_FromLong(M->m), PyInt_FromLong(M->n)); \
    if(!out_shape) {  PyErr_SetString(PyExc_RuntimeError, "Could not extract M->m or M->n"); return NULL; }; \

    PyObject* out_nnz = PyInt_FromLong(M->nz);
    if(!out_nnz) {  PyErr_SetString(PyExc_RuntimeError, "Could not extract M->nz"); return NULL; };

    /* call the class inside the csr module */
    PyObject* out_mat = PyObject_CallMethodObjArgs(scipy_mod, PyString_FromString((char *) "coo_matrix"), out_all, out_shape, NULL);

   Py_DECREF(out_indx);
   Py_DECREF(out_all);

    if(out_mat)
    {
      _sn_check_nnz(&out_mat, M);
      _sn_check_shape(&out_mat, M);

      return out_mat;
    }
    else
    {
      PyErr_SetString(PyExc_RuntimeError, "Could not create coo matrix");
      return NULL;
    }

  }

}

  static PyObject* NM_to_python(NumericsMatrix* m)
  {
  if (m)
  {
    npy_intp dims[2];
    dims[0] = m->size0;
    dims[1] = m->size1;
    if (m->matrix0)
    {
      PyObject *obj = FPyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, m->matrix0);
      PyArrayObject *array = (PyArrayObject*) obj;
      if (!array) { PyErr_SetString(PyExc_RuntimeError, "Could not create an array from matrix0. Please file a bug"); return NULL; }
      return obj;
    }
    else if(m->matrix1)
    {
      // matrix is sparse : return opaque pointer
      return SWIG_NewPointerObj(SWIG_as_voidptr(m->matrix1), $descriptor(SparseBlockStructuredMatrix *), 0);
    }
    else if(m->matrix2)
    {
      switch(m->matrix2->origin)
      {
      case NS_CSC:
      {
       return cs_sparse_to_csc_matrix(NM_csc(m), false);
      }
      case NS_TRIPLET:
      {
        return cs_sparse_to_coo_matrix(NM_triplet(m), false);
      }
      case NS_CSR:
      {
        return cs_sparse_to_csr_matrix(NM_csr(m), false);
      }
      default:
      {
        PyErr_SetString(PyExc_RuntimeError, "The given sparse matrix has no valid origin. Please file a bug");
        return NULL;
      }
      }
    }
    else
    {
      PyErr_SetString(PyExc_RuntimeError, "The given matrix is of unknown type. Please file a bug");
      return NULL;
    }
  }
  else
  {
     Py_INCREF(Py_None);
     return  Py_None;
  }
  }


} // end fragment NumericsMatrix


