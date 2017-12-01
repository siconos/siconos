%inline %{
#define ALLOC_CTRL_I 0x1
#define ALLOC_CTRL_P 0x2

  static void NM_clean_cs(CSparseMatrix* m, int alloc_ctrl)
  {
    assert(m);
    if (alloc_ctrl & ALLOC_CTRL_P) { assert(m->p); free(m->p); }
    if (alloc_ctrl & ALLOC_CTRL_I) { assert(m->i); free(m->i); }
    // We do not own any data (we steal it from a numpy array)
    m->p = NULL;
    m->i = NULL;
    m->x = NULL;
  }

  static int NM_clean(NumericsMatrix* M, int alloc_ctrl)
  {
    switch (M->storageType)
    {
    case NM_DENSE:
    {
      // We do not own the data
      M->matrix0 = NULL;
      break;
    }
    case NM_SPARSE:
    {
      assert(M->matrix2);
      switch (M->matrix2->origin)
      {
      case NS_CSC:
      {
        NM_clean_cs(M->matrix2->csc, alloc_ctrl);
        free(M->matrix2->csc);
        M->matrix2->csc = NULL;
        break;
      }
      case NS_CSR:
      {
        NM_clean_cs(M->matrix2->csr, alloc_ctrl);
        free(M->matrix2->csr);
        M->matrix2->csr = NULL;
        break;
      }
      case NS_TRIPLET:
      {
        NM_clean_cs(M->matrix2->triplet, alloc_ctrl);
        free(M->matrix2->triplet);
        M->matrix2->triplet = NULL;
        break;
      }
      default:
      {
        SWIG_Error(SWIG_RuntimeError, "The origin of the sparse matrix is unknown!");
        return 0;
      }
      }
      if (M->matrix2->trans_csc) { free(M->matrix2->trans_csc); M->matrix2->trans_csc = NULL; }
      if (M->matrix2->csc) { free(M->matrix2->csc); M->matrix2->csc = NULL; }
      if (M->matrix2->csr) { free(M->matrix2->csr); M->matrix2->csr = NULL;}
      if (M->matrix2->triplet) { free(M->matrix2->triplet); M->matrix2->triplet = NULL;}
      NM_clearSparse(M);
      break;
    }
    case NM_SPARSE_BLOCK:
    {
      // We do not own the data
      M->matrix1 = NULL;
      break;
    }
    default:
    {
      SWIG_Error(SWIG_RuntimeError, "NM_clean: unknown matrix storageType!");
      return 0;
    }
    }
    return 1;

  }
%}



#ifdef SWIGPYTHON
%include Numerics_NM_python.i
#endif /* SWIGPYTHON */

#ifdef SWIGMATLAB
%include Numerics_NM_matlab.i
#endif /* SWIGMATLAB */

// Numpy array -> NumericsMatrix
%typemap(in, fragment="NumericsMatrix") (NumericsMatrix*) 
 // free in typemap(freearg)
 (SN_ARRAY_TYPE* array_ = NULL,
 int array_ctrl_ = 0,
 SN_ARRAY_TYPE* array_i_ = NULL,
 int array_i_ctrl_ = 0,
 SN_ARRAY_TYPE* array_p_ = NULL,
 int array_p_ctrl_ = 0,
 int alloc_ctrl_ = 0,
 NumericsMatrix *nummat = NULL)
{
#ifdef SWIGPYTHON
   $1 = NM_convert_from_python($input, &nummat, &array_, &array_ctrl_, &array_i_, &array_i_ctrl_, &array_p_, &array_p_ctrl_, &alloc_ctrl_);
#endif /* SWIGPYTHON */

#ifdef SWIGMATLAB
  $1 = NM_convert_from_matlab($input, &nummat, &array_, &array_ctrl_, &array_i_, &array_i_ctrl_, &array_p_, &array_p_ctrl_, &alloc_ctrl_);
#endif /* SWIGMATLAB */
   if (!$1) { SWIG_fail; }
}


%typemap(memberin) (NumericsMatrix*) {
 //  %typemap(memberin) (NumericsMatrix*)
 // perform a deep copy
 if (!$1) { $1 = NM_create($input->storageType, $input->size0, $input->size1); }
 NM_copy($input, $1);
}

%typemap(freearg) (NumericsMatrix*) {
  // %typemap(freearg) (NumericsMatrix*)
  target_mem_mgmt(array_ctrl_$argnum,  array_$argnum);
  target_mem_mgmt(array_i_ctrl_$argnum,  array_i_$argnum);
  target_mem_mgmt(array_p_ctrl_$argnum,  array_p_$argnum);

  if (nummat$argnum)
  {
    if (!NM_clean(nummat$argnum, alloc_ctrl_$argnum)) { return SN_SWIG_ERROR_CODE; }
    NM_free(nummat$argnum);
    free(nummat$argnum);
  }

}

%typemap(out, fragment="NumericsMatrix") (NumericsMatrix*) {
  if (strcmp("$symname", "new_NumericsMatrix"))
  {
#ifdef SWIGPYTHON
    $result = NM_to_python($1);
#elif defined(SWIGMATLAB)
    $result = NM_to_matlab($1);
#endif
    if (!$result) SWIG_fail;
  }
  else
  {
    $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), $descriptor(NumericsMatrix *), SWIG_POINTER_NEW |  0 );
  }
}

%typemap(freearg) (double *z)
{
 
}

// XXX unused ???
//%typemap(out) (double* q) {

//  if ($1)
//  {
//    SN_OBJ_TYPE* out
//    C_to_target_lang1(out, arg1->size, $1);
//    $result = out;
//  }
//  else
//   {
//     Py_INCREF(Py_None);
//     $result = Py_None;
//   }
// }

// SBM handling

%typemap(in) (const SparseBlockStructuredMatrix* const A, double *denseMat)
(int dims[2])
{
  int res1=0;
  void* temp = NULL;
  res1 = SWIG_ConvertPtr($input, &temp, $1_descriptor, 0 |  0);
  if (!SWIG_IsOK(res1))
  {
    %argument_fail(res1, "$type", $symname, $argnum);
  };
  SparseBlockStructuredMatrix* A = (SparseBlockStructuredMatrix *)temp;
  assert(A);
  dims[0] = A->blocksize0[A->blocknumber0-1];
  dims[1] = A->blocksize1[A->blocknumber1-1];
  $1 = A;

  $2 = (double *) malloc(dims[0] * dims[1] * sizeof(double));
}

// FIX: do not work
%newobject SBM_to_dense(const SparseBlockStructuredMatrix* const A, double *denseMat);
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
  
  SN_OBJ_TYPE* out;
  C_to_target_lang2(out, arg1->blocksize0[arg1->blocknumber0-1], arg1->blocksize1[arg1->blocknumber1-1], $1, SWIG_fail);
  $result = out;
}

// conversion string -> FILE
%typemap(in) (FILE *file=NULL)
{
  // %typemap(in) (FILE *file)
  int alloc = 1;
  char* cstr;
  int res = SWIG_AsCharPtrAndSize($input, &cstr, NULL, &alloc);
  if (!SWIG_IsOK(res)) {
    SWIG_Error(SWIG_ArgError(res), "in method unknown', argument " "1"" of type '" "char *""'");
  }
  $1 = fopen(cstr, "r");
  if (!$1)
  {
    SWIG_Error(SWIG_IOError, format_msg_concat("in method '" "$symname" "' cannot fopen file", cstr));
    if (alloc == SWIG_NEWOBJ) free(cstr);
    SWIG_fail;
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

  $result = SWIG_AppendOutput($result,
                                     SWIG_NewPointerObj(SWIG_as_voidptr($1), $1_descriptor, SWIG_POINTER_OWN));
}

#ifdef __cplusplus
// SWIG gives syntax error for CS_NAME(_sparse)
#ifdef SICONOS_INT64
#define CSparseMatrix cs_dl_sparse
#else
#define CSparseMatrix cs_di_sparse
#endif
#else
// SWIG gives syntax error for CS_NAME(_sparse)
#ifdef SICONOS_INT64
#define CSparseMatrix struct cs_dl_sparse
#else
#define CSparseMatrix struct cs_di_sparse
#endif
#endif


%typemap(in) (const SparseBlockStructuredMatrix* const A, CSparseMatrix *outSparseMat)
{
  void *swig_arp;
  int swig_res = SWIG_ConvertPtr($input,&swig_arp,$1_descriptor, 0 | 0);

  if (SWIG_IsOK(swig_res))
  {
    $1 = (SparseBlockStructuredMatrix*) swig_arp;
    $2 = (CSparseMatrix*) malloc(sizeof(CSparseMatrix));
    if(!$2) SWIG_fail;

    SBM_to_sparse_init_memory($1,$2);
  }
  else
    SWIG_fail;
}


%typemap(argout) (CSparseMatrix *outSparseMat)
{
#ifdef SWIGPYTHON
  SN_OBJ_TYPE* csrm = cs_sparse_to_csr_matrix($1, true);
#endif /* SWIGPYTHON */

#ifdef SWIGMATLAB
  SN_OBJ_TYPE* csrm = cs_sparse_to_matlab(NM_csr_to_csc($1), true);
#endif /* SWIGMATLAB */

  if (!csrm) { SWIG_fail; }
  $result = SWIG_AppendOutput($result, csrm);
  cs_spfree($1);
}

%typemap(out) (CSparseMatrix *)
{
#ifdef SWIGPYTHON
  if ($1->nz == -1)
  {
    $result = cs_sparse_to_csc_matrix($1, true);
  }
  else if ($1->nz == -2)
  {
    $result = cs_sparse_to_csr_matrix($1, true);
  }
  else if ($1->nz >= 0)
  {
    $result = cs_sparse_to_coo_matrix($1, true);
  }
  else
  {
    SWIG_exception_fail(SWIG_RuntimeError, "The given sparse matrix is of unknown type. Please file a bug");
  }
#endif /* SWIGPYTHON */

#ifdef SWIGMATLAB
  $result = cs_sparse_to_matlab($1, true);
#endif /* SWIGMATLAB */

 if (!$result) { SWIG_fail; }
}

%typemap(in) (CSparseMatrix*)
   (int array_data_ctrl_ = 0,
   int array_i_ctrl_ = 0,
   int array_p_ctrl_ = 0,
   SN_ARRAY_TYPE *array_data_ = NULL,
   SN_ARRAY_TYPE *array_i_ = NULL,
   SN_ARRAY_TYPE *array_p_ = NULL,
   int alloc_ctrl_ = 0,
   CSparseMatrix *M = NULL)
{
#ifdef SWIGPYTHON
  int res = cs_convert_from_scipy_sparse($input, &M, &array_data_, &array_data_ctrl_, &array_i_, &array_i_ctrl_, &array_p_, &array_p_ctrl_, &alloc_ctrl_);
#endif /* SWIGPYTHON */

#ifdef SWIGMATLAB
  int res = cs_convert_from_matlab($input, &M, &array_data_, &array_data_ctrl_, &array_i_, &array_i_ctrl_, &array_p_, &array_p_ctrl_, &alloc_ctrl_);
#endif /* SWIGMATLAB */

  if (!res) { SWIG_fail; }
  else if (res < 0) { SWIG_exception_fail(SWIG_RuntimeError, "Error the matrix is not sparse!"); }

  $1 = M;
}

%typemap(memberin) (CSparseMatrix*)
{
 // perform a deep copy
 if (!$1) { $1 = NM_csparse_alloc_for_copy($input); }
 NM_copy_sparse($input, $1);
}

%typemap(freearg) (CSparseMatrix*)
{

  target_mem_mgmt(array_data_ctrl_$argnum,  array_data_$argnum);
  target_mem_mgmt(array_i_ctrl_$argnum,  array_i_$argnum);
  target_mem_mgmt(array_p_ctrl_$argnum,  array_p_$argnum);

  if(M$argnum) { NM_clean_cs(M$argnum, alloc_ctrl_$argnum); cs_spfree(M$argnum); }
}

/* %inline */
/* %{ */
/*   static void getSBM(SparseBlockStructuredMatrix* M, SparseBlockStructuredMatrix* outSBM) */
/*   { */
/*     outSBM=M; */
/*   } */
/* %} */

#undef CSparseMatrix

// issue with typemap out and is useless for now
// convert matrix to scipy.sparse.csc and do the job there
%ignore SBM_row_to_dense;
