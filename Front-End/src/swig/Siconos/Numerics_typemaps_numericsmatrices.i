// Numpy array -> NumericsMatrix (dense storage only!)
%typemap(in) 
(NumericsMatrix* A) 
(PyArrayObject* array=NULL, 
 int is_new_object=0,

 // free in typemap(freearg)
 NumericsMatrix *nummat = (NumericsMatrix *) malloc(sizeof(NumericsMatrix)))
{
  array = obj_to_array_fortran_allow_conversion($input, NPY_DOUBLE,&is_new_object);

  if (!array || !require_dimensions(array,2) ||
      !require_native(array) || !require_fortran(array)) SWIG_fail;

  nummat->storageType = 0;
  nummat->size0 =  array_size(array,0);
  nummat->size1 =  array_size(array,1);

  nummat->matrix0 = (double *)array_data(array);
  $1 = nummat;
}

%typemap(freearg) (double *z)
{
 
}

%typemap(freearg) (NumericsMatrix* A) {
  // %typemap(freearg) (NumericsMatrix* A)
  free(nummat$argnum);
  if (is_new_object$argnum && array$argnum)
  { Py_DECREF(array$argnum); }
}

%typemap(out) (NumericsMatrix* M) {
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
    $result = SWIG_NewPointerObj(SWIG_as_voidptr($1->matrix1), SWIGTYPE_p_SparseBlockStructuredMatrix, 0);
  }
  else
    SWIG_fail;
 }

%typemap(out) (double* q) {
  npy_intp dims[2];

  dims[0] = arg1->size;
  dims[1] = 1;
  if (arg1->q)
  {
    PyObject *obj = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, arg1->q);
    PyArrayObject *array = (PyArrayObject*) obj;
    if (!array || !require_fortran(array)) SWIG_fail;
    $result = obj;
  }
  else
    SWIG_fail;
 }

%apply (NumericsMatrix *A) { (NumericsMatrix *m) };
%apply (NumericsMatrix *A) { (NumericsMatrix *M) };

// SBM handling


%typemap(in) (const SparseBlockStructuredMatrix* const A) 
  (npy_intp dims[2])
{
  int res1=0;
  int newmem = 0;
  void* temp = 0;
  res1 = SWIG_ConvertPtr($input, &temp, SWIGTYPE_p_SparseBlockStructuredMatrix, 0 |  0);
  if (!SWIG_IsOK(res1)) 
  {
    SWIG_exception_fail(SWIG_ArgError(res1), 
                        "in method '" 
                        BOOST_PP_STRINGIZE($symname)
                        "', argument " "1"" of type '" "SparseBlockStructuredMatrix *""'");
  };
  SparseBlockStructuredMatrix* A = (SparseBlockStructuredMatrix *) 0;
  A = reinterpret_cast< SparseBlockStructuredMatrix * >(temp);
  assert(A);
  dims[0] = A->blocksize0[A->blocknumber0-1];
  dims[1] = A->blocksize0[A->blocknumber0-1];
  $1 = A;
}

%typemap(in, numinputs=0) (double *denseMat)
{
  // %typemap(in, numinputs=0) (double *denseMat)
  // before %typemap(in) (const SparseBlockStructuredMatrix* const A) 
  // ... but
}

%typemap(check) (double *denseMat) 
{
  // yes! ...  %typemap(check) (double *denseMat) 
  // before %typemap(in) (const SparseBlockStructuredMatrix* const A) 
  // ... what a mess
  $1 = (double *) malloc(dims1[0] * dims1[1] * sizeof(double));
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
    SWIG_exception_fail(SWIG_IOError, 
                        "in method '" 
                        BOOST_PP_STRINGIZE($symname)
                        "cannot fopen file");
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
                                     SWIG_NewPointerObj(SWIG_as_voidptr($1), 
                                                        SWIGTYPE_p_SparseBlockStructuredMatrix, 0));
}


#define GET_INT(OBJ,VAR)                                                \
  if(!PyInt_Check(OBJ))                                                 \
  {                                                                     \
    throw(std::invalid_argument(BOOST_PP_STRINGIZE(OBJ: expecting an int))); \
  }                                                                     \
  VAR = PyInt_AsLong(OBJ)

// need a PySequence_Check
#define GET_INTS(OBJ,INDEX,VAR)                                          \
  PyObject *_TEMP##VAR = PySequence_GetItem(OBJ,INDEX);                 \
  if (!PyInt_Check(_TEMP##VAR))                                         \
  {                                                                     \
    Py_XDECREF(_TEMP##VAR);                                             \
    throw(std::invalid_argument(BOOST_PP_STRINGIZE(OBJ: expecting an int))); \
  }                                                                     \
  VAR = PyInt_AsLong(_TEMP##VAR);                                   \
  Py_DECREF(_TEMP##VAR)


%typemap(in) (const SparseBlockStructuredMatrix* const A, cs_sparse *outSparseMat)
{
  void *swig_arp;
  int swig_res = SWIG_ConvertPtr($input,&swig_arp,SWIGTYPE_p_SparseBlockStructuredMatrix, 0 | 0);

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


%typemap(argout) (cs_sparse *outSparseMat)
{

  cs_sparse *M=$1;

  /* get sys.modules dict */
  PyObject* sys_mod_dict = PyImport_GetModuleDict();
  
  /* get the csr module object */
  PyObject* csr_mod = PyMapping_GetItemString(sys_mod_dict, (char *)"scipy.sparse.csr");
  
  npy_intp this_M_x_dims[1];
  this_M_x_dims[0] = M->nzmax;

  npy_intp this_M_i_dims[1];
  this_M_i_dims[0] = M->nzmax;

  npy_intp this_M_p_dims[1];
  this_M_p_dims[0] = M->m+1;

  PyObject* out_data = PyArray_SimpleNewFromData(1,this_M_x_dims,NPY_DOUBLE,M->x);
  if(!out_data) SWIG_fail;

  PyObject* out_indices = PyArray_SimpleNewFromData(1,this_M_i_dims,NPY_INT,M->i);
  if(!out_indices) SWIG_fail;

  PyObject* out_indptr = PyArray_SimpleNewFromData(1,this_M_p_dims,NPY_INT,M->p);
  if(!out_indptr) SWIG_fail;

  PyObject* out_shape = PyTuple_Pack(2,PyInt_FromLong(M->n),PyInt_FromLong(M->m));
  if(!out_shape) SWIG_fail;

  PyObject* out_nnz = PyInt_FromLong(M->nzmax);
  if(!out_nnz) SWIG_fail;

  /* call the class inside the csr module */
  %#if PY_MAJOR_VERSION < 3
  PyObject* out_csr = PyObject_CallMethodObjArgs(csr_mod, PyString_FromString((char *)"csr_matrix"), out_shape, NULL);
  %#else
  PyObject* out_csr = PyObject_CallMethodObjArgs(csr_mod, PyUnicode_FromString((char *)"csr_matrix"), out_shape, NULL);
  %#endif

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

    $result = SWIG_Python_AppendOutput($result,out_csr);
  }
  else
  {
    SWIG_fail;
  }

}

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

    int dim0, dim1, nzmax;
    GET_INTS(shape_,0,dim0);
    GET_INTS(shape_,1,dim1);
//      GET_INT(nnz,nzmax); fail: type is numpy.int32!
    nzmax = PyInt_AsLong(nnz_);
    
    array_data_ = obj_to_array_allow_conversion(data_, NPY_DOUBLE, &is_new_object1);
    array_indices_ = obj_to_array_allow_conversion(indices_, NPY_INT32, &is_new_object2);
    array_indptr_ = obj_to_array_allow_conversion(indptr_, NPY_INT32, &is_new_object3);
    
    M->m = dim0;
    M->n = dim1;
    
    M->nzmax = nzmax;
    
    M->nz = -2; // csr only for the moment
    
    M->p = (int *) malloc((M->m+1) * sizeof(int));
    if(!M->p) SWIG_fail;

    M->i = (int *) malloc(M->nzmax * sizeof(int));
    if(!M->i) SWIG_fail;

    M->x = (double *) malloc(M->nzmax * sizeof(double));
    if(!M->x) SWIG_fail;

    for(unsigned int i = 0; i < (M->m+1); i++)
    {
      M->p[i] = ((int *) array_data(array_indptr_)) [i];
    }
    
    for(unsigned int i = 0; i < M->nzmax; i++)
    {
      M->i[i] = ((int *) array_data(array_indices_)) [i];
    }
    
    memcpy(M->x, (double *) array_data(array_data_), M->nzmax * sizeof(double));
    
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
    freeSparse(M$argnum);
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


