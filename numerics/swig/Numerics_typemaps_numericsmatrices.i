// Numpy array -> NumericsMatrix (dense storage only!)
%typemap(in) 
(NumericsMatrix* A) 
(PyArrayObject* array=NULL, 
 int is_new_object=0,
 int res = 0,
 void* argp = NULL,

 // free in typemap(freearg)
 NumericsMatrix *nummat = NULL)
{
  res = SWIG_ConvertPtr($input, &argp, $descriptor(NumericsMatrix *), %convertptr_flags);
  if (SWIG_IsOK(res))
  {
    $1 = %reinterpret_cast(argp, NumericsMatrix *);
  }
  else
  {
    nummat = newNumericsMatrix();
    array = obj_to_array_fortran_allow_conversion($input, NPY_DOUBLE,&is_new_object);

    if (!array || !require_dimensions(array,2) ||
        !require_native(array) || !require_fortran(array)) SWIG_fail;

    nummat->storageType = NM_DENSE;
    nummat->size0 =  array_size(array,0);
    nummat->size1 =  array_size(array,1);

    nummat->matrix0 = (double *)array_data(array);
    $1 = nummat;
  }
}

%typemap(freearg) (double *z)
{
 
}

%typemap(freearg) (NumericsMatrix* A) {
  // %typemap(freearg) (NumericsMatrix* A)
  if (nummat$argnum) { free(nummat$argnum); }
  if (is_new_object$argnum && array$argnum)
  { Py_DECREF(array$argnum); }
}

%typemap(out) (NumericsMatrix* M) {
  if ($1)
  {
    npy_intp dims[2];
    dims[0] = $1->size0;
    dims[1] = $1->size1;
    if ($1->matrix0)
    {
      PyObject *obj = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, $1->matrix0);
      PyArrayObject *array = (PyArrayObject*) obj;
      if (!array || !require_fortran(array)) SWIG_fail;
      $result = obj;
    }
    else if($1->matrix1)
    {
      // matrix is sparse : return opaque pointer
      $result = SWIG_NewPointerObj(SWIG_as_voidptr($1->matrix1), $descriptor(SparseBlockStructuredMatrix *), 0);
    }
    else
      SWIG_fail;
   }
   else
   {
     Py_INCREF(Py_None);
     $result = Py_None;
   }
  }

%typemap(out) (double* q) {
  npy_intp dims[2];

  dims[0] = arg1->size;
  dims[1] = 1;
  if ($1)
  {
    PyObject *obj = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, $1);
    PyArrayObject *array = (PyArrayObject*) obj;
    if (!array || !require_fortran(array)) SWIG_fail;
    $result = obj;
  }
  else
   {
     Py_INCREF(Py_None);
     $result = Py_None;
   }
 }

%apply (NumericsMatrix *A) { (NumericsMatrix *m) };
%apply (NumericsMatrix *A) { (NumericsMatrix *M) };
%apply (NumericsMatrix *A) { (NumericsMatrix *H) };

// SBM handling

%typemap(in) (const SparseBlockStructuredMatrix* const A, double *denseMat)
(npy_intp dims[2])
{
  int res1=0;
  void* temp = NULL;
  res1 = SWIG_ConvertPtr($input, &temp, $1_descriptor, 0 |  0);
  if (!SWIG_IsOK(res1))
  {
    SWIG_exception_fail(SWIG_ArgError(res1), 
                        "in method '" 
                        " $symname "
                        "', argument " "1"" of type '" "SparseBlockStructuredMatrix *""'");
  };
  SparseBlockStructuredMatrix* A = (SparseBlockStructuredMatrix *) 0;
  A = reinterpret_cast< SparseBlockStructuredMatrix * >(temp);
  assert(A);
  dims[0] = A->blocksize0[A->blocknumber0-1];
  dims[1] = A->blocksize0[A->blocknumber0-1];
  $1 = A;

  $2 = (double *) malloc(dims[0] * dims[1] * sizeof(double));
}

// FIX: do not work
%newobject SBMtoDense(const SparseBlockStructuredMatrix* const A, double *denseMat);
%typemap(newfree) (double *denseMat)
{
  // %typemap(newfree) (double *denseMat)
  if($1)
  {
    free($1);
  }
}


%typemap(argout) (double *denseMat) 
{
  
  PyObject* pyarray = FPyArray_SimpleNewFromData(2,                
                                                 dims1,                  
                                                 NPY_DOUBLE,            
                                                 $1);
  $result = pyarray;
}

// conversion python string -> FILE
%typemap(in) (FILE *file=NULL)
{
  // %typemap(in) (FILE *file)
  %#if PY_MAJOR_VERSION < 3
  $1 = fopen(PyString_AsString($input), "r");
  %#else
  PyObject* tmp_ascii = PyUnicode_AsASCIIString($input);
  $1 = fopen(PyBytes_AsString(tmp_ascii), "r");
  Py_DECREF(tmp_ascii);
  %#endif
  if (!$1)
  {
    SWIG_exception_fail(SWIG_IOError, "in method '" "$symname" "' cannot fopen file");
  }
}

%typemap(freearg) (FILE *file)
{
  // %typemap(freearg) (FILE *file)
  if($1)
  {
    fclose($1);
  }
}

%typemap(in, numinputs=0) (SparseBlockStructuredMatrix* outSBM) 
{
  $1 = (SparseBlockStructuredMatrix*) malloc(sizeof(SparseBlockStructuredMatrix));
  if(!$1) SWIG_fail;

  $1->block = NULL;
  $1->index1_data = NULL;
  $1->index2_data = NULL;

}

%typemap(argout) (SparseBlockStructuredMatrix* outSBM)
{
  if(!$1) SWIG_fail;

  $result = SWIG_Python_AppendOutput($result,
                                     SWIG_NewPointerObj(SWIG_as_voidptr($1), $1_descriptor, 0));
}


#define GET_INT(OBJ,VAR)                                                \
  if(!PyInt_Check(OBJ))                                                 \
  {                                                                     \
    throw(std::invalid_argument(#OBJ ": expecting an int")); \
  }                                                                     \
  VAR = PyInt_AsLong(OBJ)

// need a PySequence_Check
#define GET_INTS(OBJ,INDEX,VAR)                                          \
  PyObject *_TEMP##VAR = PySequence_GetItem(OBJ,INDEX);                 \
  if (!PyInt_Check(_TEMP##VAR))                                         \
  {                                                                     \
    Py_XDECREF(_TEMP##VAR);                                             \
    throw(std::invalid_argument(#OBJ ": expecting an int")); \
  }                                                                     \
  VAR = PyInt_AsLong(_TEMP##VAR);                                   \
  Py_DECREF(_TEMP##VAR)


%typemap(in) (const SparseBlockStructuredMatrix* const A, cs_sparse *outSparseMat)
{
  void *swig_arp;
  int swig_res = SWIG_ConvertPtr($input,&swig_arp,$1_descriptor, 0 | 0);

  if (SWIG_IsOK(swig_res))
  {
    $1 = (SparseBlockStructuredMatrix*) swig_arp;
    $2 = (cs_sparse*) malloc(sizeof(cs_sparse));
    if(!$2) SWIG_fail;

    SBMtoSparseInitMemory($1,$2);
  }
  else
    SWIG_fail;
}


%inline %{

int cs_sparse_to_csc_matrix(cs_sparse *M, PyObject** csc_matrix)
{
  if (!M)
  {
    Py_INCREF(Py_None);
    *csc_matrix = Py_None;
    return 0;
  }
  else
  {
    /* get sys.modules dict */
    PyObject* sys_mod_dict = PyImport_GetModuleDict();

    /* get the csr module object */
    PyObject* csr_mod = PyMapping_GetItemString(sys_mod_dict, (char *)"scipy.sparse.csr");

    if (!csr_mod)
    {
      PyErr_SetString(PyExc_RuntimeError, "Did you import scipy.sparse.csr?");
      return 1;
    }

    npy_intp this_M_x_dims[1];
    this_M_x_dims[0] = M->nzmax;

    npy_intp this_M_i_dims[1];
    this_M_i_dims[0] = M->nzmax;

    npy_intp this_M_p_dims[1];
    this_M_p_dims[0] = M->m+1;

    PyObject* out_data = PyArray_SimpleNewFromData(1,this_M_x_dims,NPY_DOUBLE,M->x);
    if(!out_data) { PyErr_SetString(PyExc_RuntimeError, "Could not extract M->x"); return 1; };

    PyObject* out_indices = PyArray_SimpleNewFromData(1,this_M_i_dims,NPY_INT64,M->i);
    if(!out_indices) {  PyErr_SetString(PyExc_RuntimeError, "Could not extract M->i"); return 1; };

    PyObject* out_indptr = PyArray_SimpleNewFromData(1,this_M_p_dims,NPY_INT64,M->p);
    if(!out_indptr) {  PyErr_SetString(PyExc_RuntimeError, "Could not extract M->p"); return 1; };

    /* Warning ! m is the number of rows, n the number of columns ! --xhub */
    PyObject* out_shape = PyTuple_Pack(2,PyInt_FromLong(M->m),PyInt_FromLong(M->n));
    if(!out_shape) {  PyErr_SetString(PyExc_RuntimeError, "Could not extract M->m or M->n"); return 1; };

    PyObject* out_nnz = PyInt_FromLong(M->nzmax);
    if(!out_nnz) {  PyErr_SetString(PyExc_RuntimeError, "Could not extract M->nzmax"); return 1; };

    /* call the class inside the csr module */
#if PY_MAJOR_VERSION < 3
    PyObject* out_csr = PyObject_CallMethodObjArgs(csr_mod, PyString_FromString((char *)"csr_matrix"), out_shape, NULL);
#else
    PyObject* out_csr = PyObject_CallMethodObjArgs(csr_mod, PyUnicode_FromString((char *)"csr_matrix"), out_shape, NULL);
#endif

    if(out_csr)
    {
      PyObject_SetAttrString(out_csr,"data",out_data);
      PyObject_SetAttrString(out_csr,"indices",out_indices);
      PyObject_SetAttrString(out_csr,"indptr",out_indptr);

#ifndef NDEBUG
      PyObject *auto_nnz = PyObject_GetAttrString(out_csr,"nnz");
      assert(PyInt_AsLong(auto_nnz) == M->nzmax);
      Py_XDECREF(auto_nnz);
#endif

      *csc_matrix = SWIG_Python_AppendOutput(*csc_matrix, out_csr);
      return 0;
    }
    else
    {
      PyErr_SetString(PyExc_RuntimeError, "Could not create csr matrix");
      return 1;
    }
  }

}

%}


%typemap(argout) (cs_sparse *outSparseMat)
{
  if (cs_sparse_to_csc_matrix($1, &$result)) SWIG_fail;

}

%typemap(out) (cs_sparse *)
{
  if (cs_sparse_to_csc_matrix($1, &$result)) SWIG_fail;

}

%define %SAFE_CAST_INT(pyvar, len, indvar, dest_array)
{
    int array_pyvartype_ = PyArray_TYPE((PyArrayObject *)pyvar);
    switch (array_pyvartype_)
    {
      case NPY_INT32:
      {
        PyArrayObject* array_pyvar = obj_to_array_allow_conversion(pyvar, NPY_INT32, &indvar);
        if (!array_pyvar) { SWIG_fail; }

        for(unsigned int i = 0; i < len; i++)
        {
          dest_array[i] = ((int32_t *) array_data(array_pyvar)) [i];
        }
        break;
      }
      case NPY_INT64:
      {
        PyArrayObject* array_pyvar = obj_to_array_allow_conversion(pyvar, NPY_INT64, &indvar);
        if (!array_pyvar) { SWIG_fail; }

        for(unsigned int i = 0; i < len; i++)
        {
          dest_array[i] = ((int64_t *) array_data(array_pyvar)) [i];
        }
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
        SWIG_fail;
      }
    }
}
%enddef

%typemap(in) (cs_sparse* M) 
  (PyObject *shape_ = NULL,
   PyObject *nnz_ = NULL,
   PyObject *data_ = NULL,
   PyObject *indices_ = NULL,
   PyObject *indptr_ = NULL,
   int is_new_object1=0, 
   int is_new_object2=0,
   int is_new_object3=0,
   PyArrayObject *array_data_ = NULL,
   PyArrayObject *array_indices_ = NULL,
   PyArrayObject *array_indptr_ = NULL,
   cs_sparse *M = NULL)
{
  try
  {
    M = (cs_sparse *) malloc(sizeof(cs_sparse));      
    if(!M) SWIG_fail;

    PyObject *obj = $input;
    
    shape_ = PyObject_GetAttrString(obj,"shape");
    nnz_ = PyObject_GetAttrString(obj,"nnz");
    data_ = PyObject_GetAttrString(obj,"data");
    indices_ = PyObject_GetAttrString(obj,"indices");
    indptr_ = PyObject_GetAttrString(obj,"indptr");

    unsigned int dim0 = 0, dim1 = 0, nzmax = 0;
    GET_INTS(shape_,0,dim0);
    GET_INTS(shape_,1,dim1);
//      GET_INT(nnz,nzmax); fail: type is numpy.int32!
    nzmax = PyInt_AsLong(nnz_);
    
    M->m = dim0;
    M->n = dim1;
    
    M->nzmax = nzmax;
    
    M->nz = -2; // csr only for the moment
    
    M->p = (csi *) malloc((M->m+1) * sizeof(csi));
    M->i = (csi *) malloc(M->nzmax * sizeof(csi));
    M->x = (double *) malloc(M->nzmax * sizeof(double));

    // if we return NULL here, the matrix M is going to be freed. It will miserably fail if
    // M has not been "completly" initialized --xhub
    if(!M->p) SWIG_fail;
    if(!M->i) SWIG_fail;
    if(!M->x) SWIG_fail;

    array_data_ = obj_to_array_allow_conversion(data_, NPY_DOUBLE, &is_new_object1);
    if (!array_data_) { SWIG_fail; }

    memcpy(M->x, (double *) array_data(array_data_), M->nzmax * sizeof(double));

    %SAFE_CAST_INT(indptr_, dim0 + 1, is_new_object2, M->p);
    %SAFE_CAST_INT(indices_, nzmax, is_new_object3, M->i);

    $1 = M;
  }
  catch (const std::invalid_argument& e)
  {
    SWIG_exception(SWIG_ValueError, e.what());
  }
}

%typemap(freearg) (cs_sparse* M)
{

  Py_XDECREF(shape_$argnum);
  Py_XDECREF(nnz_$argnum);
  Py_XDECREF(data_$argnum);
  Py_XDECREF(indices_$argnum);
  Py_XDECREF(indptr_$argnum);
  
  if (array_data_$argnum && is_new_object1$argnum)
  {
    Py_DECREF(array_data_$argnum);
  }
  
  if (array_indptr_$argnum && is_new_object2$argnum)
  {
    Py_DECREF(array_indptr_$argnum);
  }
  
  if (array_indices_$argnum && is_new_object3$argnum)
  {
    Py_DECREF(array_indices_$argnum);
  }
  if(M$argnum)
  {
    cs_spfree(M$argnum);
  }
}

%apply (cs_sparse *M) {(const cs_sparse const * m)};

%apply (cs_sparse *M) {(cs_sparse const * m)};

%apply (cs_sparse *M) {(const cs_sparse * m)};

%apply (cs_sparse *M) {(cs_sparse * m)};

%apply (cs_sparse *M) {(const cs_sparse const *sparseMat)};

%apply (cs_sparse *M) {(cs_sparse *sparseMat)};


%inline
%{
  void getSBM(SparseBlockStructuredMatrix* M, SparseBlockStructuredMatrix* outSBM)
  {
    outSBM=M;
  }
%}

// issue with typemap out and is useless for now
// convert matrix to scipy.sparse.csc and do the job there
%ignore SBMRowToDense;
