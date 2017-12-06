// XXX this is a hack --xhub
#undef SICONOS_INT64

%inline %{
#include "SiconosConfig.h"
%}

%define %SAFE_CAST_INT(pyvar, len, dest_array, alloc)
{
%#ifdef SICONOS_INT64
//        PyErr_Warn(PyExc_UserWarning, "Performance warning: the vector of indices or pointers is in int32, but siconos has 64-bits integers: we have to perform a conversion. Consider given sparse matrix in the right format");
        dest_array = (int64_t*) malloc(len * sizeof(int64_t));
        if(!dest_array) { SWIG_Error(SWIG_RuntimeError, "Allocation of i or p failed (triggered by conversion to int32)"); return 0; }
        for(unsigned i = 0; i < len; ++i)
        {
          dest_array[i] = (int32_t) pyvar[i];
        }
        alloc = true;
%#else
        dest_array = (int32_t *) pyvar;
%#endif
}
%enddef

%fragment("NumericsMatrix", "header")
{
  static int cs_convert_from_matlab(mxArray* obj, CSparseMatrix** m, mxArray** array_data_, int* array_data_ctrl_, mxArray** array_i_, int* array_i_ctrl_, mxArray** array_p_, int* array_p_ctrl_, int* alloc_ctrl)
  {

  assert(m);

  bool isspmat = mxIsSparse(obj);
  if (isspmat == false)
  {
    return -1;
  }
  else
  {

    // csc
    CSparseMatrix* M = (CSparseMatrix*) calloc(1, sizeof(CSparseMatrix));
    if(!M) { SWIG_Error(SWIG_RuntimeError, "Failed to allocate a cs_sparse"); return 0; }

    M->m = mxGetM(obj);
    M->n = mxGetN(obj);
    M->nzmax = mxGetNzmax(obj);
    *m = M;
    M->nz = -1;

    double* data_ = (double*)mxGetData(obj);
    mwIndex* row_ = mxGetIr(obj);
    mwIndex* col_ = mxGetJc(obj);

    M->x = data_;

    bool alloc_p = false;
    %SAFE_CAST_INT(col_, (M->n + 1), M->p, alloc_p);
    if (alloc_p) { *alloc_ctrl |= ALLOC_CTRL_P; };
    bool alloc_i = false;
    %SAFE_CAST_INT(row_, M->nzmax, M->i, alloc_i);
    if (alloc_i) { *alloc_ctrl |= ALLOC_CTRL_I; };

    return 1;
  }
  }

  static int NM_convert_from_matlab2(mxArray* obj, NumericsMatrix* m, mxArray** array_data_, int* array_data_ctrl_, mxArray** array_i_, int* array_i_ctrl_, mxArray** array_p_, int* array_p_ctrl_, int* alloc_ctrl)
  {
    CSparseMatrix* csm = NULL;
    int res = cs_convert_from_matlab(obj, &csm, array_data_, array_data_ctrl_, array_i_, array_i_ctrl_, array_p_, array_p_ctrl_, alloc_ctrl);
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
        SWIG_Error(SWIG_RuntimeError, "Unknown CSparseMatrix from cs_convert_from_matlab");
        return 0;
      }

      NM_update_size(m);
    }

    return res;
  }


  static NumericsMatrix* NM_convert_from_matlab(mxArray* obj, NumericsMatrix** tmpmat, mxArray** array_data_, int* array_ctrl, mxArray** array_i_, int* array_i_ctrl_, mxArray** array_p_, int* array_p_ctrl_, int* alloc_ctrl)
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
    if (!mxIsNumeric(obj) || array_numdims(obj) != 2)
    {
      SWIG_Error(SWIG_TypeError, "The given object does not have the right structure. We expect a 2 dimensional array");
      goto fail;
    }

    if (!mxIsSparse(obj))
    {

      out->storageType = NM_DENSE;
      out->size0 =  array_size(obj, 0);
      out->size1 =  array_size(obj, 1);
      out->matrix0 = (double *)array_data(obj);

      *array_data_ = obj;
    }
    else
    {
      int sp_conv = NM_convert_from_matlab2(obj, out, array_data_, array_ctrl, array_i_, array_i_ctrl_, array_p_, array_p_ctrl_, alloc_ctrl);
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
          SWIG_exception_fail(SWIG_TypeError, "Cannot build a NumericsMatrix from the given object");
        }
      }
    }
  }

  return out;

fail:
  if (*tmpmat) { free(*tmpmat); *tmpmat = NULL; }
  return NULL;
  }

static mxArray* cs_sparse_to_matlab(CSparseMatrix *M, bool copy)
{
    if(M->nz == -1)
    {
      mxArray* out_mat = mxCreateSparse(M->m, M->n, M->nzmax, mxREAL);
      if (!out_mat)
      {
        SWIG_Error(SWIG_RuntimeError, "cs_sparse_to_matlab :: Could not create sparse matrix, not enough memory");
        return NULL;
      }

      double* outx  = mxGetPr(out_mat);
      mwIndex* outi = mxGetIr(out_mat);
      mwIndex* outp = mxGetJc(out_mat);

      double* Mx = M->x;
      CS_INT* Mi = M->i;
      CS_INT* Mp = M->p;

      size_t nnz = Mp[M->n];
      memcpy(outx, Mx, nnz * sizeof(double));

      for (size_t i = 0; i < nnz; ++i)
      {
        outi[i] = (mwIndex)Mi[i];
      }

      for (size_t i = 0; i <= M->n; ++i)
      {
        outp[i] = (mwIndex)Mp[i];
      }

      return out_mat;
    }
    else
    {
      SWIG_Error(SWIG_RuntimeError, "cs_sparse_to_matlab :: Could not create sparse matrix");
      return NULL;
    }

}

  static mxArray* NM_to_matlab(NumericsMatrix* m)
  {
  if (m)
  {
    if (m->matrix0)
    {
      mxArray *obj;
      C_to_target_lang2(obj, m->size0,  m->size1, m->matrix0, return NULL);
      return obj;
    }
    else if(m->matrix1)
    {
      // matrix is sparse : return opaque pointer
      return SWIG_NewPointerObj(SWIG_as_voidptr(m->matrix1), $descriptor(SparseBlockStructuredMatrix *), 0);
    }
    else if(m->matrix2)
    {
      return cs_sparse_to_matlab(NM_csc(m), false);
    }
    else
    {
      SWIG_Error(SWIG_RuntimeError, "The given matrix is of unknown type. Please file a bug");
      return NULL;
    }
  }
  else
  {
     return NULL;
  }
  }


} // end fragment NumericsMatrix


