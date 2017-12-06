 // Matrices
%include "SparseMatrix.h"
%include "SparseBlockMatrix.h"
%include "NumericsMatrix.h"
%include "NumericsSparseMatrix.h"

%define %NM_convert_from_target(input, output, ACTION)
{
 SN_ARRAY_TYPE* array_ = NULL;
 int array_ctrl_ = 0;
 SN_ARRAY_TYPE* array_i_ = NULL;
 int array_i_ctrl_ = 0;
 SN_ARRAY_TYPE* array_p_ = NULL;
 int array_p_ctrl_ = 0;
 int alloc_ctrl_ = 0;
 NumericsMatrix *nummat = NULL;

#if defined(SWIGPYTHON)
   NumericsMatrix* Mtmp = NM_convert_from_python(input, &nummat, &array_, &array_ctrl_, &array_i_, &array_i_ctrl_, &array_p_, &array_p_ctrl_, &alloc_ctrl_);
#endif /* SWIGPYTHON */

#if defined(SWIGMATLAB)
  NumericsMatrix* Mtmp = NM_convert_from_matlab(input, &nummat, &array_, &array_ctrl_, &array_i_, &array_i_ctrl_, &array_p_, &array_p_ctrl_, &alloc_ctrl_);
#endif /* SWIGMATLAB */

   if (!Mtmp) { ACTION; }

  *(output) = NM_create(Mtmp->storageType, Mtmp->size0, Mtmp->size1);
  NM_copy(Mtmp, *(output));

  target_mem_mgmt(array_ctrl_,  array_);
  target_mem_mgmt(array_i_ctrl_,  array_i_);
  target_mem_mgmt(array_p_ctrl_,  array_p_);

  if (nummat)
  {
    if (!NM_clean(nummat, alloc_ctrl_)) { ACTION; }
    free(nummat);
  }

}
%enddef
// some extensions but numpy arrays should be used instead
%extend NumericsMatrix
{
%fragment("NumericsMatrix");

  NumericsMatrix(SN_OBJ_TYPE* o)
  {
    NumericsMatrix *M;
    %NM_convert_from_target(o, &M, return NULL)
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

#ifdef SWIGPYTHON
  SN_OBJ_TYPE * __setitem__(SN_OBJ_TYPE* index, double v)
  {
    int i, j;
    if (!self->matrix0)
    {
      PyErr_SetString(PyExc_RuntimeError, "The given matrix is not dense (matrix0 == NULL). For now only items on dense matrices can be set.");
      return NULL;
    }
    if (!PyArg_ParseTuple(index, "ii:NumericsMatrix__setitem__",&i,&j)) return NULL;
    NumericsMatrix_set_matrix0(self,i,j,v);
    return Py_BuildValue("");
  }

  SN_OBJ_TYPE * __getitem__(SN_OBJ_TYPE * index)
  {
    int i, j;
    if (!self->matrix0)
    {
      PyErr_SetString(PyExc_RuntimeError, "The given matrix is not dense (matrix0 == NULL). For now only items on dense matrices can be requested.");
      return NULL;
    }
    if (!PyArg_ParseTuple(index, "ii:NumericsMatrix__getitem__",&i,&j)) return NULL;
    return SWIG_From_double(NumericsMatrix_get_matrix0(self,i,j));
  }

  int __len__()
  {
    return self->size0 * self->size1;
  }

// useful? -- xhub
//  SN_OBJ_TYPE * __str__()
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
#endif /* SWIGPYTHON */

  ~NumericsMatrix()
  {
    NM_free($self);
    free($self);
  }

};

%extend SparseBlockStructuredMatrix
{
 ~SparseBlockStructuredMatrix()
 {
   SBM_free($self);
   free($self);
 }
}


