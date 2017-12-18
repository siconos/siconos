

#ifdef SWIGPYTHON
%inline %{

#define CHECK_ARRAY(X) \
!require_native(X)

#define CHECK_ARRAY_VECTOR(X) \
CHECK_ARRAY(X) || !require_contiguous(X) || !require_dimensions(X, 1)

#define CHECK_ARRAY_MATRIX(X) \
CHECK_ARRAY(X) || !require_fortran(X) || !require_dimensions(X, 2)

#define CHECK_ARRAY_SIZE(req_size, array, indx) (req_size == array_size(array, indx))

#define obj_to_sn_array(obj, alloc) obj_to_array_fortran_allow_conversion(obj, NPY_DOUBLE, alloc);
#define obj_to_sn_vector(obj, alloc) obj_to_array_contiguous_allow_conversion(obj, NPY_DOUBLE, alloc);
%}
#endif /* SWIGPYTHON */

#ifdef SWIGMATLAB
%inline %{

#include "NM_conversions.h"
static inline long int array_size(mxArray* m, int indx)
{
  const mwSize *dim_array;
  dim_array = mxGetDimensions(m);
  return (long int)dim_array[indx];
}

#define array_numdims(X) (int)mxGetNumberOfDimensions(X)

#define CHECK_ARRAY(X) false

#define CHECK_ARRAY_VECTOR(X) !(array_numdims(X) == 2 && ((array_size(X, 0) == 1 && (array_size(X, 1) > 0)) || (array_size(X, 0) > 0 && (array_size(X, 1) == 1))))

#define CHECK_ARRAY_MATRIX(X) !(array_numdims(X) == 2 && array_size(X, 0) > 0 && array_size(X, 1) > 0)

#define CHECK_ARRAY_SIZE(req_size, array, indx) (req_size == array_size(array, indx))

#define array_data(X) mxGetData(X)

// XXX maybe we need to copy stuff here? -- xhub
#define obj_to_sn_array(obj, dummy) (mxIsDouble(obj)) ? obj : NULL
#define obj_to_sn_vector(obj, dummy) (mxIsDouble(obj) && !mxIsSparse(obj)) ? obj : NULL
#define obj_to_sn_vector_int(obj, dummy) sizeof(int) == 8 ? (mxIsInt32(obj) ? obj : NULL) : (mxIsInt64(obj) ? obj : NULL)
%}
#endif /* SWIGMATLAB */

// vectors of size problem_size from given *Problem as first input
// no conversion => inout array XXX FIX issue here
%typemap(in) (double *z) (SN_ARRAY_TYPE* array = NULL, int is_new_object = 0) {

  array = obj_to_sn_vector($input, &is_new_object);

  if (!array)
  {
   SWIG_exception_fail(SWIG_TypeError, "Could not get a SN_ARRAY_TYPE from the python object");
  }

  if (CHECK_ARRAY_VECTOR(array))
  {
   SWIG_exception_fail(SWIG_TypeError, "The given object does not have the right structure. We expect a vector (or list, tuple, ...)");
  }

  $1 = (double *) array_data(array);

 }

%typemap(freearg) (double *z)
{
  target_mem_mgmt(is_new_object$argnum, array$argnum);
}

// list of matrices problemSizex3
%typemap(in) (double *blocklist3x3) (SN_ARRAY_TYPE* array=NULL, int is_new_object=0) {

  array = obj_to_sn_array($input, &is_new_object);

  if (!array)
  {
   SWIG_exception_fail(SWIG_TypeError, "Could not get a SN_ARRAY_TYPE from the python object");
  }

  if (CHECK_ARRAY(array))
  {
   SWIG_exception_fail(SWIG_TypeError, "The given object does not have the right structure. We expect a vector (or list, tuple, ...)");
  }

  SN_ARRAY_INT_TYPE array_len[2] = {0,0};

  if (! *p_problem_size1)
  {
    if (array_numdims(array) > 1)
    {
      *p_problem_size1 = fmax(array_size(array,0), array_size(array,1)) / 3;
    }
    else if (array_numdims(array) == 1)
    {
      *p_problem_size1 = array_size(array,0) / 3;
    }
    else
      SWIG_fail;

    if (*p_problem_size1 % 3 != 0) SWIG_fail;

    if (*p_problem_size1 / 3 == 0) SWIG_fail;

    /* number_of_contacts1 = *p_problem_size1 / 3; */

  }

  assert (*p_problem_size1);


  if (array_numdims(array) > 1)
  {
    array_len[0] = *p_problem_size1 * 3;
    array_len[1] = 1;
  }
  else
  {
    array_len[0] = *p_problem_size1 * 3;
  }

  if (CHECK_ARRAY_VECTOR(array) || !CHECK_ARRAY_SIZE(array_len[0], array, 0) || !(array_numdims(array) > 1 && CHECK_ARRAY_SIZE(array_len[1], array, 1))) SWIG_fail;

  $1 = (double *) array_data(array);

 }

%typemap(freearg) (double *blocklist3x3)
{
  target_mem_mgmt(is_new_object$argnum, array$argnum)
}


// matrices problemSizexproblemSize
%typemap(in) (double *blockarray3x3) (SN_ARRAY_TYPE* array=NULL, int is_new_object=0) {

  array = obj_to_sn_array($input, &is_new_object);
  if (!array)
  {
   SWIG_exception_fail(SWIG_TypeError, "Could not get a SN_ARRAY_TYPE from the python object");
  }

  if (CHECK_ARRAY(array))
  {
   SWIG_exception_fail(SWIG_TypeError, "The given object does not have the right structure. We expect a vector (or list, tuple, ...)");
  }


  SN_ARRAY_INT_TYPE array_len[2] = {0,0};

  if (! *p_problem_size1)
  {
    if (array_numdims(array) > 1)
    {
      *p_problem_size1 = array_size(array,0); // assume square matrix
    }
    else if (array_numdims(array) == 1)
    {
      *p_problem_size1 = isqrt(array_size(array,0));
    }
    else
      SWIG_fail;


    if (*p_problem_size1 % 3 != 0) SWIG_fail;

    if (*p_problem_size1 / 3 == 0) SWIG_fail;

    /* number_of_contacts1 = *p_problem_size1 / 3; */

  }

  assert (*p_problem_size1);

  if (array_numdims(array) > 1)
  {
    array_len[0] = *p_problem_size1;
    array_len[1] = *p_problem_size1;
  }
  else
  {
    array_len[0] = *p_problem_size1 * *p_problem_size1;
  }

  if (CHECK_ARRAY_VECTOR(array)  || (array_numdims == 2 && CHECK_ARRAY_MATRIX(array)) || !CHECK_ARRAY_SIZE(array_len[0], array, 0) || !(array_numdims(array) > 1 && CHECK_ARRAY_SIZE(array_len[1], array, 1))) SWIG_fail;

  $1 = (double *) array_data(array);

 }

%typemap(freearg) (double *blockarray3x3)
{
  if (is_new_object$argnum && array$argnum)
    { Py_DECREF(array$argnum); }
}

// vectors of size problem_size from problemSize as first input
%typemap(in) (double *blocklist3) (SN_ARRAY_TYPE* array=NULL, int is_new_object=0) {

  array = obj_to_sn_vector($input, &is_new_object);
  if (!array)
  {
   SWIG_exception_fail(SWIG_TypeError, "Could not get a SN_ARRAY_TYPE from the python object");
  }

  if (CHECK_ARRAY(array))
  {
   SWIG_exception_fail(SWIG_TypeError, "The given object does not have the right structure. We expect a vector (or list, tuple, ...)");
  }


  SN_ARRAY_INT_TYPE array_len[2] = {0,0};

  if (! *p_problem_size1)
  {
    if (array_numdims(array) == 1)
    {
      *p_problem_size1 = array_size(array,0);
    }
    else if (array_numdims(array) > 1)
    {
      *p_problem_size1 = fmax(array_size(array,0), array_size(array,1));
    }

    if (*p_problem_size1 % 3 != 0) SWIG_fail;

    if (*p_problem_size1 / 3 == 0) SWIG_fail;


    /* number_of_contacts1 = *p_problem_size1 / 3; */

  }

  assert (*p_problem_size1);

  if (array_numdims(array) == 1)
  {
    array_len[0] = *p_problem_size1;
  }
  else
  {
    array_len[0] = *p_problem_size1;
    array_len[1] = 1;
  }

  if (CHECK_ARRAY_VECTOR(array) || !CHECK_ARRAY_SIZE(array_len[0], array, 0) || !(array_numdims(array) > 1 && CHECK_ARRAY_SIZE(array_len[1], array, 1))) SWIG_fail;

  $1 = (double *) array_data(array);

 }

%typemap(freearg) (double *blocklist3)
{
  target_mem_mgmt(is_new_object$argnum, array$argnum);
}

// vectors of size problem_size

// 1 : numinputs=0 mandatory to avoid arg
%typemap(in, numinputs=0) (double *output_blocklist3) (SN_OBJ_TYPE* array=NULL)
{
    // %typemap(in, numinputs=0)
    // we cannot get problem_size here as numinputs=0 => before
    // numinputs=1, how can we change this ??
}

// 2 : check must be done after in
%typemap(check) (double *output_blocklist3)
{
  if (*p_problem_size1)
  {
    C_to_target_lang2_alloc($1, array$argnum, *p_problem_size1, 1, SWIG_fail)
  }

}

// 3 : return arg
%typemap(argout) (double *output_blocklist3)
{
  if (*p_problem_size1)
  {
     $result = SWIG_AppendOutput($result,(SN_OBJ_TYPE *)array$argnum);
  }

}

// 3x3 matrices

// 1 : numinputs=0 mandatory to avoid arg
%typemap(in, numinputs=0) (double *output_blocklist3x3) (SN_OBJ_TYPE* array=NULL)
{
    // %typemap(in, numinputs=0)
    // we cannot get problem_size here as numinputs=0 => before
    // numinputs=1, how can we change this ??
}

%typemap(in, numinputs=0) (double *output_blockarray3x3) (SN_OBJ_TYPE* array=NULL)
{
    // %typemap(in, numinputs=0)
    // we cannot get problem_size here as numinputs=0 => before
    // numinputs=1, how can we change this ??
}

// 2 : check must be done after in
%typemap(check) (double *output_blocklist3x3)
{
  if (*p_problem_size1)
  {
    C_to_target_lang2_alloc($1, array$argnum, (*p_problem_size1) * 3, 1, SWIG_fail)
  }

}

%typemap(check) (double *output_blockarray3x3)
{
  if (*p_problem_size1)
  {
    C_to_target_lang2_alloc($1, array$argnum, (*p_problem_size1), *p_problem_size1, SWIG_fail)
  }

}

// 3 : return arg
%typemap(argout) (double *output_blocklist3x3)
{
  if (*p_problem_size1)
  {
    $result = SWIG_AppendOutput($result,(SN_OBJ_TYPE *)array$argnum);
  }

}

%typemap(argout) (double *output_blockarray3x3)
{
  if (*p_problem_size1)
  {
    $result = SWIG_AppendOutput($result,(SN_OBJ_TYPE *)array$argnum);
  }
}

%typemap(in, numinputs=0) (double *output3) 
{
}

%typemap(argout) (double *output3)
{
  $result = SWIG_AppendOutput($result, SWIG_From_double($1[0]));
  $result = SWIG_AppendOutput($result, SWIG_From_double($1[1]));
  $result = SWIG_AppendOutput($result, SWIG_From_double($1[2]));
}

// other names that must be transformed this way
%apply (double *z) { (double *r) };

%apply (double *z) { (double *u) };

%apply (double *z) { (double *w) };

%apply (double *z) { (double *x) };

%apply (double *z) { (double *F) };

%apply (double *z) { (double *Fmcp) };

%apply (double *z) { (double *zlem) };

%apply (double *z) { (double *wlem) };

%apply (double *z) { (double *reaction) };

%apply (double *z) { (double *velocity) };

%apply (double *z) { (double *globalVelocity) };

%apply (double *z) { (double *mu) };

%apply (double *z) { (double *q) };

%apply (double *z) { (double *b) };

%apply (double *blocklist3) { (double *vect3D) };

%apply (double *blocklist3) { (double *reaction3D) };

%apply (double *blocklist3) { (double *velocity3D) };

%apply (double *blocklist3) { (double *rho3D) };

%apply (double *blocklist3) { (double *blocklist3_1) };

%apply (double *blocklist3) { (double *blocklist3_2) };

%apply (double *blocklist3) { (double *blocklist3_3) };

%apply (double *blocklist3x3) { (double *blocklist3x3_1) };

%apply (double *blocklist3x3) { (double *blocklist3x3_2) };

%apply (double *blocklist3x3) { (double *blocklist3x3_3) };

%apply (double *output_blocklist3) { (double *output_blocklist3_1) };

%apply (double *output_blocklist3) { (double *output_blocklist3_2) };

%apply (double *output_blocklist3) { (double *output_blocklist3_3) };

%apply (double *output_blocklist3) { (double *output_blocklist3_4) };

%apply (double *output_blocklist3x3) { (double *output_blocklist3x3_1) };

%apply (double *output_blocklist3x3) { (double *output_blocklist3x3_2) };
