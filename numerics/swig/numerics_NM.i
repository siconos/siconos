 // Matrices
%include "SparseMatrix.h"
%include "SparseBlockMatrix.h"
%include "NumericsMatrix.h"
%include "NumericsSparseMatrix.h"

// some extensions but numpy arrays should be used instead
%extend NumericsMatrix_
{
%fragment("NumericsMatrix");

  NumericsMatrix_(PyObject* o)
  {
    PyArrayObject* array_ = NULL;
    int array_ctrl_ = 0;
    PyArrayObject* array_i_ = NULL;
    int array_i_ctrl_ = 0;
    PyArrayObject* array_p_ = NULL;
    int array_p_ctrl_ = 0;
    int alloc_ctrl_ = 0;

    NumericsMatrix *tmpmat = NULL;
    NumericsMatrix *Mtmp = NM_convert_from_python(o, &tmpmat, &array_, &array_ctrl_, &array_i_, &array_i_ctrl_, &array_p_, &array_p_ctrl_, &alloc_ctrl_);

    if (!Mtmp) { return NULL; }

    NumericsMatrix *M = createNumericsMatrix(Mtmp->storageType, Mtmp->size0, Mtmp->size1);
    NM_copy(Mtmp, M);

    if (array_ctrl_ && array_) { Py_DECREF(array_); }
    if (array_i_ctrl_ && array_i_) { Py_DECREF(array_i_); }
    if (array_p_ctrl_ && array_p_) { Py_DECREF(array_p_); }

    if (tmpmat)
    {
      if (!NM_clean(tmpmat, alloc_ctrl_)) { return NULL; }
      free(tmpmat);
    }

    return M;
  }

  void set_matrix0(int i, int j, double v)
  {
    assert(self->matrix0);
    self->matrix0[i+j*self->size1] = v;
  }

  double get_matrix0(int i, int j)
  {
    assert(self->matrix0);
    return self->matrix0[i+j*self->size1];
  }


  PyObject * __setitem__(PyObject* index, double v)
  {
    int i, j;
    if (!self->matrix0)
    {
      PyErr_SetString(PyExc_RuntimeError, "The given matrix is not dense (matrix0 == NULL). For now only items on dense matrices can be set.");
      return NULL;
    }
    if (!PyArg_ParseTuple(index, "ii:NumericsMatrix___setitem__",&i,&j)) return NULL;
    NumericsMatrix_set_matrix0(self,i,j,v);
    return Py_BuildValue("");
  }

  PyObject * __getitem__(PyObject * index)
  {
    int i, j;
    if (!self->matrix0)
    {
      PyErr_SetString(PyExc_RuntimeError, "The given matrix is not dense (matrix0 == NULL). For now only items on dense matrices can be requested.");
      return NULL;
    }
    if (!PyArg_ParseTuple(index, "ii:NumericsMatrix___getitem__",&i,&j)) return NULL;
    return SWIG_From_double(NumericsMatrix_get_matrix0(self,i,j));
  }

  int __len__()
  {
    return self->size0 * self->size1;
  }

// useful? -- xhub
//  PyObject * __str__()
//  {
//    if (!self->matrix0)
//    {
//      PyErr_SetString(PyExc_RuntimeError, "The given matrix is not dense (matrix0 == NULL). Only dense matrix can be displayed");
//      return NULL;
//    }
//    std::stringstream result;
//    result << "[ ";
//    for (int i=0; i < self->size0; ++i)
//      {
//        if (i > 0) result << "  ";
//        result << "[";
//        for (int j=0; j < self->size1; ++j)
//          {
//            result << " " << NumericsMatrix_get_matrix0(self,i,j);
//            if (j < self->size1-1) result << ",";
//          }
//        result << " ]";
//        if (i < self->size0-1) result << "," << std::endl;
//      }
//    result << " ]" << std::endl;
//    %#if PY_MAJOR_VERSION < 3
//    return PyString_FromString(result.str().c_str());
//    %#else
//    return PyUnicode_FromString(result.str().c_str());
//    %#endif
//  }
//

  ~NumericsMatrix_()
  {
    freeNumericsMatrix($self);
    free($self);
  }

};


#define GET_ATTR(OBJ,ATTR)                                              \
  ATTR = PyObject_GetAttrString(OBJ,#ATTR);          \
  if(!ATTR)                                                             \
  {                                                                     \
    throw(std::invalid_argument("need a " #ATTR "attr")); \
  }

/* I'm not sure it's a good idea to duplicate thise here ... it is already defined in csparse.h */
typedef struct cs_sparse    /* matrix in compressed-column or triplet form */
{
  csi nzmax ;	    /* maximum number of entries */
  csi m ;	    /* number of rows */
  csi n ;	    /* number of columns */
  csi *p ;	    /* column pointers (size n+1) or col indices (size nzmax) */
  csi *i ;	    /* row indices, size nzmax */
  double *x ;	    /* numerical values, size nzmax */
  csi nz ;	    /* # of entries in triplet matrix, -1 for compressed-col */
} cs ;

%extend cs_sparse
{
%fragment("NumericsMatrix");

 // from a scipy.sparse matrix
 cs_sparse(PyObject *obj)
 {
   int array_data_ctrl_ = 0;
   int array_i_ctrl_ = 0;
   int array_p_ctrl_ = 0;
   int alloc_ctrl_ = 0;
   PyArrayObject *array_data_ = NULL, *array_i_ = NULL, *array_p_ = NULL;
   CSparseMatrix* M = NULL;

   CSparseMatrix Mtmp;
   CSparseMatrix* Mtmp_p = &Mtmp;

   int res = cs_convert_from_scipy_sparse(obj, &Mtmp_p, &array_data_, &array_data_ctrl_, &array_i_, &array_i_ctrl_, &array_p_, &array_p_ctrl_, &alloc_ctrl_);

   if (!res) { SWIG_fail; }
   else if (res < 0) { PyErr_SetString(PyExc_RuntimeError, "Error the matrix is not sparse!"); goto fail; }

   M = (CSparseMatrix *) malloc(sizeof(CSparseMatrix));
   if(!M) { PyErr_SetString(PyExc_RuntimeError, "Failed to allocate a cs_sparse"); goto fail; }

   // perform a deep copy since we do not have any mechanism to record the fact we use data from the python object
   NM_copy_sparse(Mtmp_p, M);

fail:
   if (array_data_ && array_data_ctrl_) { Py_DECREF(array_data_); }
   if (array_i_ && array_i_ctrl_) { Py_DECREF(array_i_); }
   if (array_p_ && array_p_ctrl_) { Py_DECREF(array_p_); }

   NM_clean_cs(Mtmp_p, alloc_ctrl_);

   return M;

 }

 ~cs_sparse()
 {
  cs_spfree($self);
 }
}

%extend SparseBlockStructuredMatrix_
{
 ~SparseBlockStructuredMatrix_()
 {
   freeSBM($self);
   free($self);
 }
}
