// -*- C++ -*-
// Siconos-Front-End , Copyright INRIA 2005-2012.
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.	
// Siconos is a free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// Siconos is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Siconos; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
//
// Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr 
//

// Siconos.i - SWIG interface for Siconos
%module(directors="1") Numerics

%include start.i

%{
#include "SiconosNumerics.h"
#include "NumericsConfig.h"
#include "SolverOptions.h"
#include "SparseMatrix.h"
#include "SparseBlockMatrix.h"
#include "FrictionContact3D_Solvers.h"
#include "Friction_cst.h"
#include "frictionContact_test_function.h"
#include "FrictionContact3D_AlartCurnier.h"
#include "FrictionContact3D_localAlartCurnier.h"
#include "FrictionContact3D_compute_error.h"
#include "fclib_interface.h"
#include "Numerics_functions.h"
#include "box.h"

#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/cat.hpp>

#include <limits>
%}

#ifdef WITH_IO
%{
#include <SiconosFullNumerics.hpp>
%}
#endif
%include "picklable.i"

%include "std_string.i"

 // needed macros
 %include "NumericsConfig.h"



 // more convenient
 %rename (LCP) LinearComplementarityProblem;
 %rename (MLCP) MixedLinearComplementarityProblem;
 %rename (MCP) MixedComplementarityProblem;
 %rename (VI) VariationalInequality;

 %ignore lcp_compute_error_only;

 %apply (int DIM1 , int DIM2 , double* IN_FARRAY2)
 {(int nrows, int ncols, double* data          )};

 %apply (int DIM1 , double* INPLACE_ARRAY1)
 {(int sizex, double *x          )};

 %apply (int DIM1 , double* INPLACE_ARRAY1)
 {(int sizey, double *y          )};

 %apply (int DIM1 , double* INPLACE_ARRAY1)
 {(int sizeq, double *q         )};

 %apply (int DIM1 , double* INPLACE_ARRAY1)
 {(int sizez, double *z        )};

 %apply (int DIM1 , double* INPLACE_ARRAY1)
 {(int sizew, double *w        )};

 %apply (double IN_ARRAY1[ANY])
 {(double reaction[3])}

 %apply (double IN_ARRAY1[ANY])
 {(double velocity[3])}

 %apply (double IN_ARRAY1[ANY])
 {(double rho[3])}

 %apply (double ARGOUT_ARRAY1[ANY])
 {(double ACresult[3])}

 %apply (double ARGOUT_ARRAY1[ANY])
 {(double result[3])}

 %apply (double ARGOUT_ARRAY1[ANY])
 {(double A[9])}

 %apply (double ARGOUT_ARRAY1[ANY])
 {(double B[9])}


 // Handle standard exceptions
 %include "exception.i"
 %exception
 {
   try
   {
     $action
   }
   catch (const std::invalid_argument& e)
   {
     SWIG_exception(SWIG_ValueError, e.what());
   }
   catch (const std::out_of_range& e)
   {
     SWIG_exception(SWIG_IndexError, e.what());
   }
 }

%include solverOptions.i



 %typemap(in) (LinearComplementarityProblem*) (npy_intp problem_size) {
   void *lcp;
   int res = SWIG_ConvertPtr($input, &lcp,SWIGTYPE_p_LinearComplementarityProblem, 0 |  0 );
   if (!SWIG_IsOK(res)) SWIG_fail;

   problem_size=((LinearComplementarityProblem *) lcp)->size;

   $1 = (LinearComplementarityProblem *) lcp;
 }

 %typemap(in) (MixedLinearComplementarityProblem*) (npy_intp mlcproblem_size) {
   void *mlcp;
   int res = SWIG_ConvertPtr($input, &mlcp,SWIGTYPE_p_MixedLinearComplementarityProblem, 0 |  0 );
   if (!SWIG_IsOK(res)) SWIG_fail;

   mlcproblem_size=((MixedLinearComplementarityProblem *) mlcp)->n +((MixedLinearComplementarityProblem *) mlcp)->m;

   $1 = (MixedLinearComplementarityProblem *) mlcp;
 }

 // problemSize given as first arg => set by first numpy array length
 // in remaining args
 %typemap(in, numinputs=0) 
   (unsigned int problemSize) 
   (unsigned int *p_problem_size, npy_intp number_of_contacts)
 {
   // the first array length sets problemSize
   p_problem_size = &$1;
   *p_problem_size = 0;
   number_of_contacts = 0;  
 }

%typemap(in) 
  (FrictionContactProblem*) 
  (npy_intp problem_size, npy_intp problem_dimension, npy_intp number_of_contacts) 
{
  void *fcp;
  int res = SWIG_ConvertPtr($input, &fcp,SWIGTYPE_p_FrictionContactProblem, 0 |  0 );
  if (!SWIG_IsOK(res)) SWIG_fail;

  problem_dimension=((FrictionContactProblem *) fcp)->dimension;
  number_of_contacts=((FrictionContactProblem *) fcp)->numberOfContacts;
  problem_size=((FrictionContactProblem *) fcp)->numberOfContacts * problem_dimension;

  $1 = (FrictionContactProblem*) fcp;
}

%typemap(in) 
  (GlobalFrictionContactProblem*) 
  (npy_intp problem_size, npy_intp problem_dimension, npy_intp number_of_contacts) 
{
  void *fcp;
  int res = SWIG_ConvertPtr($input, &fcp,SWIGTYPE_p_FrictionContactProblem, 0 |  0 );
  if (!SWIG_IsOK(res)) SWIG_fail;

  problem_dimension=((GlobalFrictionContactProblem *) fcp)->dimension;
  number_of_contacts=((GlobalFrictionContactProblem *) fcp)->numberOfContacts;
  problem_size=((GlobalFrictionContactProblem *) fcp)->numberOfContacts * problem_dimension;

  $1 = (GlobalFrictionContactProblem*) fcp;
}

// vectors of size problem_size from given *Problem as first input
// no conversion => inout array XXX FIX issue here
%typemap(in) (double *z) (PyArrayObject* array=NULL, int is_new_object=0) {
  //typemap(in) (double *z) (PyArrayObject* array=NULL, int is_new_object=0)

  array = obj_to_array_allow_conversion($input, NPY_DOUBLE, &is_new_object);

  if (!array 
      || !require_native(array) ) 
    SWIG_fail;
  
  npy_intp array_len[2] = {0,0};

  array_len[0] = array_size(array,0);

  if (array_numdims(array) > 1)
    array_len[1] = array_size(array,1);

  if(!require_size(array, array_len, array_numdims(array))) 
    SWIG_fail;
  
  $1 = (double *) array_data(array);

 }

%typemap(freearg) (double *z)
{
  if (is_new_object$argnum && array$argnum)
    { Py_DECREF(array$argnum); }
}

// list of matrices problemSizex3
%typemap(in) (double *blocklist3x3) (PyArrayObject* array=NULL, int is_new_object=0) {

  array = obj_to_array_contiguous_allow_conversion($input, NPY_DOUBLE,&is_new_object);

  npy_intp array_len[2] = {0,0};

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
    
    number_of_contacts1 = *p_problem_size1 / 3;
    
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
    
  if (!array 
      || !require_native(array) || !require_contiguous(array)
      || !require_size(array, array_len, array_numdims(array))) SWIG_fail;
  
  $1 = (double *) array_data(array);

 }

%typemap(freearg) (double *blocklist3x3)
{
  if (is_new_object$argnum && array$argnum)
    { Py_DECREF(array$argnum); }
}


// matrices problemSizexproblemSize
%typemap(in) (double *blockarray3x3) (PyArrayObject* array=NULL, int is_new_object=0) {

  array = obj_to_array_fortran_allow_conversion($input, NPY_DOUBLE,&is_new_object);

  npy_intp array_len[2] = {0,0};

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
    
    number_of_contacts1 = *p_problem_size1 / 3;
    
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

  if (!array 
      || !require_native(array) || !require_fortran(array)
      || !require_size(array, array_len, array_numdims(array))) SWIG_fail;
  
  $1 = (double *) array_data(array);

 }

%typemap(freearg) (double *blockarray3x3)
{
  if (is_new_object$argnum && array$argnum)
    { Py_DECREF(array$argnum); }
}

// vectors of size problem_size from problemSize as first input
%typemap(in) (double *blocklist3) (PyArrayObject* array=NULL, int is_new_object=0) {

  array = obj_to_array_contiguous_allow_conversion($input, NPY_DOUBLE,&is_new_object);

  npy_intp array_len[2] = {0,0};

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
    
    
    number_of_contacts1 = *p_problem_size1 / 3;
    
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

  if (!array 
      || !require_native(array) || !require_contiguous(array)
      || !require_size(array, array_len, array_numdims(array))) SWIG_fail;
  
  $1 = (double *) array_data(array);

 }

%typemap(freearg) (double *blocklist3)
{
  if (is_new_object$argnum && array$argnum)
    { Py_DECREF(array$argnum); }
}

// vectors of size problem_size

// 1 : numinputs=0 mandatory to avoid arg
%typemap(in, numinputs=0) (double *output_blocklist3) (PyObject* array=NULL)
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
    
    npy_intp dims[2] = { *p_problem_size1, 1};
    
    array$argnum = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (!array$argnum) SWIG_fail;
    $1 = ($1_ltype) array_data(array$argnum);
  }
  
}

// 3 : return arg
%typemap(argout) (double *output_blocklist3)
{
  if (*p_problem_size1)
  {
     $result = SWIG_Python_AppendOutput($result,(PyObject *)array$argnum);
  }
  
}

// info param
%typemap(in, numinputs=0) (int *info) (int temp_info = -1)
{
  // a default initialization : solver may stop if *info = 0 (checkTrivialCase)
  // checkTrivialCase => better if directly in solvers, not in driver.
  $1 = &temp_info;
}

// 3x3 matrices

// 1 : numinputs=0 mandatory to avoid arg
%typemap(in, numinputs=0) (double *output_blocklist3x3) (PyObject* array=NULL)
{
    // %typemap(in, numinputs=0)
    // we cannot get problem_size here as numinputs=0 => before
    // numinputs=1, how can we change this ??
}

%typemap(in, numinputs=0) (double *output_blockarray3x3) (PyObject* array=NULL)
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
    
    npy_intp dims[2] = { *p_problem_size1 * 3, 1 };
    
    array$argnum = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    // block list : require_fortran useless?
    if (!array$argnum) SWIG_fail;
    PyArrayObject *array = (PyArrayObject*) array$argnum;
    if (!array || !require_fortran(array)) SWIG_fail;
    $1 = ($1_ltype) array_data(array);
  }
  
}

%typemap(check) (double *output_blockarray3x3) 
{
  if (*p_problem_size1)
  {
    
    npy_intp dims[2] = { *p_problem_size1, *p_problem_size1};
    
    array$argnum = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    PyArrayObject *array = (PyArrayObject*) array$argnum;
    if (!array || !require_fortran(array)) SWIG_fail;
    $1 = ($1_ltype) array_data(array);
  }
  
}

// 3 : return arg
%typemap(argout) (double *output_blocklist3x3)
{
  if (*p_problem_size1)
  {
    $result = SWIG_Python_AppendOutput($result,(PyObject *)array$argnum);
  }
  
}

%typemap(argout) (double *output_blockarray3x3)
{
  if (*p_problem_size1)
  {
    $result = SWIG_Python_AppendOutput($result,(PyObject *)array$argnum);
  }
}


// vectors of size numberOfContacts
%typemap(in) (double *mu) (PyArrayObject* array=NULL, int is_new_object=0) {

  npy_intp array_len[2] = {0,0};

  if (number_of_contacts1)
  {
    array_len[0] = number_of_contacts1;
    array_len[1] = 1;

    array = obj_to_array_contiguous_allow_conversion($input, NPY_DOUBLE,&is_new_object);
    
    if (!array
        || !require_native(array) || !require_contiguous(array) || !require_fortran(array)
        || !require_size(array, array_len, array_numdims(array))) SWIG_fail;
    
    $1 = (double *) array_data(array);
    
  }
 }

%typemap(freearg) (double *mu)
{
  if (is_new_object$argnum && array$argnum)
    { Py_DECREF(array$argnum); }
}


// other names that must be transformed this way
%apply (double *z) { (double *w) }; 

%apply (double *z) { (double *x) }; 

%apply (double *z) { (double *F) }; 

%apply (double *z) { (double *Fmcp) }; 

%apply (double *z) { (double *zlem) };

%apply (double *z) { (double *wlem) };

%apply (double *z) { (double *reaction) };

%apply (double *z) { (double *velocity) };

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

%inline
%{
 unsigned int isqrt(unsigned int n)
  {
    unsigned int c = 0x8000;
    unsigned int g = 0x8000;

    for(;;)
    {
      if(g*g > n)
      g ^= c;
      c >>= 1;
      if(c == 0)
        return g;
      g |= c;
    }
  }
%}

// signatures
%feature("autodoc", 1);

// generated docstrings from doxygen xml output
%include Numerics-docstrings.i
 
%fragment("NumPy_Fragments");

// LCP
%include "LinearComplementarityProblem.h"
%include "LCP_Solvers.h"
%include "lcp_cst.h"
%include "SolverOptions.h"
%include "NumericsOptions.h"
%include "frictionContact_test_function.h"

//Relay
%include "relay_cst.h"


// redefine typemap on q for MLCP
%typemap(out) (double* q) {
  npy_intp dims[2];

  dims[0] = arg1->m + arg1->n;
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



// MLCP
%include "MixedLinearComplementarityProblem.h"
%include "MLCP_Solvers.h"
%include "mlcp_cst.h"

%inline %{


#include <stdio.h>
  FrictionContactProblem* frictionContactProblemFromFile
    (const char * filename)
  {
    FILE * finput = fopen(filename, "r");
    if (finput) 
    {
      FrictionContactProblem* problem = 
        (FrictionContactProblem *) malloc(sizeof(FrictionContactProblem));
      if (frictionContact_newFromFile(problem,finput))
      {
      fprintf(stderr, "frictionContactProblemFromFile: cannot load %s\n",filename);
      free(problem);
      return 0;
      }
      else
        return problem;
    }
    else
    {
      fprintf(stderr, "frictionContactProblemFromFile: cannot open %s\n",filename);
      return 0;
    }
    fclose(finput);
  }
  
  MixedLinearComplementarityProblem* mixedLinearComplementarityProblemFromFile
    (const char * filename)
  {
    FILE * finput = fopen(filename, "r");
    if (finput) 
    {
      MixedLinearComplementarityProblem* problem = 
        (MixedLinearComplementarityProblem *) malloc(sizeof(MixedLinearComplementarityProblem));
      if (mixedLinearComplementarity_newFromFile(problem,finput))
      {
      fprintf(stderr, "mixedLinearComplementarityProblemFromFile: cannot load %s\n",filename);
      free(problem);
      return 0;
      }
      else
        return problem;
    }
    else
    {
      fprintf(stderr, "mixedLinearComplementarityProblemFromFile: cannot open %s\n",filename);
      return 0;
    }
    fclose(finput);
  }

static PyObject *my_callback_NablaFmcp = NULL;

static PyObject * set_my_callback_NablaFmcp(PyObject *o)
{
  PyObject *result = NULL;
  PyObject *temp;
  if (!PyCallable_Check(o)) {
    PyErr_SetString(PyExc_TypeError, "parameter must be callable");
    return NULL;
  }
  Py_XINCREF(o);         /* Add a reference to new callback */
  Py_XDECREF(my_callback_NablaFmcp);  /* Dispose of previous callback */
  my_callback_NablaFmcp = o;       /* Remember new callback */
  
  /* Boilerplate to return "None" */
  Py_INCREF(Py_None);
  result = Py_None;

  return result;
}

static void  my_call_to_callback_NablaFmcp (int size, double *z, double *nablaF)
{  
//  printf("I am in my_call_to_callback_NablaFmcp (int size, double *z, double *NablaF)\n");

  npy_intp this_matrix_dim[1];
  this_matrix_dim[0]=size;
  
  PyObject* pyarray = FPyArray_SimpleNewFromData(1,this_matrix_dim, NPY_DOUBLE, z);   
  PyObject* tuple = PyTuple_New(1);
  PyTuple_SetItem(tuple, 0, pyarray);  
  PyObject* result; 

  if (PyCallable_Check(my_callback_NablaFmcp))
  {
    result = PyObject_CallObject(my_callback_NablaFmcp, tuple);
  }
  else
  {
    PyErr_SetString(PyExc_TypeError, "parameter must be callable");
  }

   // PyTuple_SetItem steals a reference to the object
   // I'm leaving this commented so that people don't do the mistake twice -- xhub
   //Py_DECREF(pyarray);
   Py_DECREF(tuple);

//Comment :  it will be safier to use  obj_to_array_fortran_allow_conversion

  if (is_array(result))
  {
    if (array_size(result,0) != size || array_size(result,1) != size )
    {
      char message[240];
      sprintf(message, "Wrong size for  the return value of callback function. Expected size is %i x %i", size,size);
      PyErr_SetString(PyExc_RuntimeError,message);
    }
    else if (array_numdims(result) != 2)
    {
      char message[240];
      sprintf(message, "Wrong dimension for  the return value of callback function. Expected dimension is 2");
      PyErr_SetString(PyExc_RuntimeError,message);
    }
    else
    { 

      int is_new_object0=0;      
      PyArrayObject* result2 = make_fortran((PyArrayObject *)result, &is_new_object0);
      // if (is_new_object0)
      // {
      //   Py_DECREF(result);
      //   printf ("the object is new !!\n");
      // }
      memcpy(nablaF, (double *)array_data(result2), size*size * sizeof(double));
      
    }
   Py_DECREF(result);
  }
  else
  {      
    const char * desired_type = typecode_string(NPY_DOUBLE);
    const char * actual_type  = pytype_string(result);
    PyErr_Format(PyExc_TypeError,
                 "Array of type '%s' required as return value fo callback function. A '%s' was returned",
                   desired_type, actual_type);
    if (result != NULL) Py_DECREF(result); // things can go really south ...
  }

  return;

}

static PyObject *my_callback_Fmcp = NULL;

static PyObject * set_my_callback_Fmcp(PyObject *o1)
{
  PyObject *result = NULL;
  PyObject *temp;
  if (!PyCallable_Check(o1)) {
    PyErr_SetString(PyExc_TypeError, "parameter must be callable");
    return NULL;
  }
  Py_XINCREF(o1);         /* Add a reference to new callback */
  Py_XDECREF(my_callback_Fmcp);  /* Dispose of previous callback */
  my_callback_Fmcp = o1;       /* Remember new callback */
  
  /* Boilerplate to return "None" */
  Py_INCREF(Py_None);
  result = Py_None;

  return result;
}

static void  my_call_to_callback_Fmcp (int size, double *z, double *F)
{

//  printf("I am in my_call_to_callback_Fmcp (int size, double *z, double *F)\n");

  npy_intp this_matrix_dim[1];
  this_matrix_dim[0]=size;
  
  PyObject* pyarray = FPyArray_SimpleNewFromData(1,this_matrix_dim, NPY_DOUBLE, z);   
  PyObject* tuple = PyTuple_New(1);
  PyTuple_SetItem(tuple, 0, pyarray);  
  PyObject* result; 

  if (PyCallable_Check(my_callback_Fmcp))
  {
    result = PyObject_CallObject(my_callback_Fmcp, tuple);
  }
  else
  {
    PyErr_SetString(PyExc_TypeError, "parameter must be callable");
  }

  // PyTuple_SetItem steals a reference to the object
  // I'm leaving this commented so that people don't do the mistake twice -- xhub
  //Py_DECREF(pyarray);
  Py_DECREF(tuple);

  if (is_array(result))
  {
    if (array_size(result,0) != size)
    {
      char message[240];
      sprintf(message, "Wrong size for  the return value of callback function. Expected size is %i", size);
      PyErr_SetString(PyExc_RuntimeError,message);
    }
    else if (array_numdims(result) != 1)
    {
      char message[240];
      sprintf(message, "Wrong dimension for  the return value of callback function. Expected dimension is 1");
      PyErr_SetString(PyExc_RuntimeError,message);
    }
    else
    { 
      int is_new_object0=0;      
      PyArrayObject* result2 = make_fortran((PyArrayObject *)result, &is_new_object0);
      memcpy(F, (double *)array_data(result2), size * sizeof(double));
    }
    Py_DECREF(result);
  }
  else
  {
    const char * desired_type = typecode_string(NPY_DOUBLE);
    const char * actual_type  = pytype_string(result);
    PyErr_Format(PyExc_TypeError,
                 "Array of type '%s' required as return value fo callback function. A '%s' was returned",
                   desired_type, actual_type);
    if (result != NULL) Py_DECREF(result); // things can go really south ...
  }
  
  return;

}




 %}



// MCP 
%include "MixedComplementarityProblem.h"
%include "MCP_Solvers.h"
%include "MCP_cst.h"

%extend NumericsOptions
{
  NumericsOptions()
  {
    NumericsOptions *numerics_options;
    numerics_options = (NumericsOptions *) malloc(sizeof(NumericsOptions));
    return numerics_options;
  }

  ~NumericsOptions()
  {
    delete($self);
  }
}

%extend SolverOptions
{

  SolverOptions(FRICTION_SOLVER id)
  {
    SolverOptions *SO;
    SO = (SolverOptions *) malloc(sizeof(SolverOptions));
    switch(id)
    {
    case SICONOS_FRICTION_2D_NSGS:
   
    case SICONOS_FRICTION_2D_PGS:
    case SICONOS_FRICTION_2D_CPG:
    case SICONOS_FRICTION_2D_LATIN:
    {
      frictionContact2D_setDefaultSolverOptions(SO, id);
      break;
    }      
    // 3D
    default:
    {
      frictionContact3D_setDefaultSolverOptions(SO, id);
    }
    }
    
    return SO;
  }

  SolverOptions()
    {
      SolverOptions *SO;
      SO = (SolverOptions *) malloc(sizeof(SolverOptions));
      return SO;
    }

  SolverOptions(LinearComplementarityProblem* lcp, LCP_SOLVER id)
  {
    SolverOptions *SO;
    SO = (SolverOptions *) malloc(sizeof(SolverOptions));
    linearComplementarity_setDefaultSolverOptions(lcp, SO, id);
    return SO;
  }
  
  SolverOptions(MixedLinearComplementarityProblem* mlcp, MLCP_SOLVER id)
  {
    SolverOptions *SO;
    SO = (SolverOptions *) malloc(sizeof(SolverOptions));
    SO->solverId=id;
    mixedLinearComplementarity_setDefaultSolverOptions(mlcp, SO);
    return SO;
  }

  SolverOptions(MixedComplementarityProblem* mlcp, MCP_SOLVER id)
  {
    SolverOptions *SO;
    SO = (SolverOptions *) malloc(sizeof(SolverOptions));
    SO->solverId=id;
    mixedComplementarity_setDefaultSolverOptions(mlcp, SO);
    return SO;
  }

  SolverOptions(MixedComplementarityProblem2* mcp, MCP_SOLVER id)
  {
    SolverOptions *SO;
    SO = (SolverOptions *) malloc(sizeof(SolverOptions));
    SO->solverId=id;
    mixedComplementarity_setDefaultSolverOptions((MixedComplementarityProblem*)mcp, SO);
    return SO;
  }

  SolverOptions(VariationalInequality* vi, VI_SOLVER id)
  {
    SolverOptions *SO;
    SO = (SolverOptions *) malloc(sizeof(SolverOptions));
    SO->solverId=id;
    variationalInequality_setDefaultSolverOptions(SO, id);
    return SO;
  }

  // SolverOptions(FRICTION_SOLVER id)
  // {
  //   return BOOST_PP_CAT(FE_SWIG_INTERNAL_MEMBER,SolverOptions_makeSolverOptions)(NULL, NULL, id);
  // }

  // SolverOptions(FrictionContactProblem* fcp, FRICTION_SOLVER id)
  // {
  //   return BOOST_PP_CAT(FE_SWIG_INTERNAL_MEMBER,SolverOptions_makeSolverOptions)(NULL, fcp, id);
  // }

  ~SolverOptions() 
    { 
      deleteSolverOptions($self);
    }

};

%extend LinearComplementarityProblem
{
  LinearComplementarityProblem(PyObject *o1, PyObject *o2)
    {

      int is_new_object1=0;
      int is_new_object2=0;
      PyArrayObject* array = obj_to_array_fortran_allow_conversion(o1, NPY_DOUBLE,&is_new_object1);
      PyArrayObject* vector = obj_to_array_contiguous_allow_conversion(o2, NPY_DOUBLE, &is_new_object2); 
      LinearComplementarityProblem *LC;
      // return pointer : free by std swig destructor
      LC = (LinearComplementarityProblem *) malloc(sizeof(LinearComplementarityProblem));
      NumericsMatrix *M = (NumericsMatrix *) malloc(sizeof(NumericsMatrix));

      M->storageType = 0;
      M->size0 = array_size(array,0);
      M->size1 = array_size(array,1);
      M->matrix0 = (double *) malloc(M->size0*M->size1*sizeof(double));
      memcpy(M->matrix0,array_data(array),M->size0*M->size1*sizeof(double));
      LC->size = M->size0;
      LC->M = M;
      LC->q = (double *) malloc(M->size0*sizeof(double));
      memcpy(LC->q,array_data(vector),M->size0*sizeof(double));

      // python mem management
      if(is_new_object1 && array)
      {
        Py_DECREF(array);
      }

      if(is_new_object2 && vector)
      {
        Py_DECREF(vector);
       }

      return LC;

    }
  

  ~LinearComplementarityProblem()
  {
    freeLinearComplementarityProblem($self);
  }

};


%exception MixedLinearComplementarityProblem {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}

%extend MixedLinearComplementarityProblem
{
  MixedLinearComplementarityProblem()
   {
     MixedLinearComplementarityProblem* MLCP;
     MLCP =  (MixedLinearComplementarityProblem *) malloc(sizeof(MixedLinearComplementarityProblem));
     return MLCP;
   }

  MixedLinearComplementarityProblem(PyObject *dim, PyObject *o1, PyObject *o2)
    {

      int is_new_object1=0;
      int is_new_object2=0;

      PyArrayObject* array = obj_to_array_fortran_allow_conversion(o1, NPY_DOUBLE,&is_new_object1);
      PyArrayObject* vector = obj_to_array_contiguous_allow_conversion(o2, NPY_DOUBLE, &is_new_object2); 
      
      
      MixedLinearComplementarityProblem *MLCP;
      // return pointer : free by std swig destructor
      MLCP = (MixedLinearComplementarityProblem *) malloc(sizeof(MixedLinearComplementarityProblem));
      NumericsMatrix *M = (NumericsMatrix *) malloc(sizeof(NumericsMatrix));

      M->storageType = 0;
      M->size0 = array_size(array,0);
      M->size1 = array_size(array,1);
      if ( array_size(array,0) !=  array_size(array,1))
      {
        PyErr_Format(PyExc_ValueError,
                     "A non square matrix (%ld,%ld) has been given",
                     array_size(array,0), array_size(array,1));
      }
      
      M->matrix0 = (double *) malloc(M->size0*M->size1*sizeof(double));
      memcpy(M->matrix0,array_data(array),M->size0*M->size1*sizeof(double));

      MLCP->n = (int) PyInt_AsLong(dim);
      MLCP->m = M->size0 - MLCP->n;
      MLCP->blocksRows = (int *) malloc(3*sizeof(int));
      MLCP->blocksIsComp = (int *) malloc(2*sizeof(int));
       
      
      MLCP->blocksRows[0]=0;
      MLCP->blocksRows[1]=MLCP->n;
      MLCP->blocksRows[2]=MLCP->n+MLCP->m;
      MLCP->blocksIsComp[0]=0;
      MLCP->blocksIsComp[1]=1;
      
 
      MLCP->isStorageType1 = 1;     
      MLCP->isStorageType2 = 0;
      MLCP->M = M;
      MLCP->A = NULL;
      MLCP->B = NULL;
      MLCP->C = NULL;
      MLCP->D = NULL;
      MLCP->a = NULL;
      MLCP->b = NULL;

      if ( array_size(array,0) !=  array_size(vector,0))
      {
        //printf("size of q = %i\n",  array_size(vector,0));
        //printf("size of M = %i\n",  array_size(array,0));
        
        PyErr_Format(PyExc_ValueError,
                     "Matrix and vector of incompatible lengths (%ld != %ld) ",
                     array_size(array,0), array_size(vector,0) );
      }
      MLCP->q = (double *) malloc(M->size0*sizeof(double));
      memcpy(MLCP->q,array_data(vector),M->size0*sizeof(double));

      // python mem management
      if(is_new_object1 && array)
      {
        Py_DECREF(array);
      }

      if(is_new_object2 && vector)
      {
        Py_DECREF(vector);
      }

      return MLCP;

    }
  
 

  ~MixedLinearComplementarityProblem()
  {
    freeMixedLinearComplementarityProblem($self);
  }
  
  // MixedLinearComplementarityProblem * newFromFilename(PyObject * o1)
  // {
  //   int result;
  //   MixedLinearComplementarityProblem *MLCP;
  //   // return pointer : free by std swig destructor
  //   MLCP = (MixedLinearComplementarityProblem *) malloc(sizeof(MixedLinearComplementarityProblem));
    
  //   char *arg1 = (char *) 0 ;
  //   int res1 ;
  //   char *buf1 = 0 ;
  //   int alloc1 = 0 ;
    
  //   res1 = SWIG_AsCharPtrAndSize(o1, &buf1, NULL, &alloc1);
  //   // if (!SWIG_IsOK(res1)) {
  //   //   SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "MixedLinearComplementarity_newFromFilename" "', argument " "1"" of type '" "char *""'");
  //   // }
  //   arg1 = reinterpret_cast< char * >(buf1);
  //   {
  //     try
  //     {
  //       result = (int)mixedLinearComplementarity_newFromFilename(MLCP,arg1);
  //     }
  //     catch (const std::invalid_argument& e)
  //     {
  //       // SWIG_exception(SWIG_ValueError, e.what());
  //     }
  //   }
    
  //   return MLCP;
    
  // }

};

// %callback("%s_cb");
// void tutu(double z);
// %nocallback;




%extend MixedComplementarityProblem
{

 

  MixedComplementarityProblem()
   {
     MixedComplementarityProblem* MCP;
     MCP =  (MixedComplementarityProblem *) malloc(sizeof(MixedComplementarityProblem));
     MCP->Fmcp=NULL;
     MCP->nablaFmcp=NULL;
     MCP->computeFmcp=NULL;
     MCP->computeNablaFmcp=NULL;
     return MCP;
   }

  void set_computeFmcp(PyObject *o)
  {
    set_my_callback_Fmcp(o);
    $self->computeFmcp = (my_call_to_callback_Fmcp);
  }
  
  void set_computeNablaFmcp(PyObject *o)
  {

    set_my_callback_NablaFmcp(o);
    $self->computeNablaFmcp = (my_call_to_callback_NablaFmcp);
  }
  
  void test_call_to_callback()
  {
    printf("I am in test_call_to_callback()\n");
    
    int size =   $self->sizeEqualities +  $self->sizeInequalities;

    double * z = (double *)malloc(size*sizeof(double));
    double * F = (double *)malloc(size*sizeof(double));
    double * nablaF = (double *)malloc(size*size*sizeof(double));
    
    for (int i=0; i < size; i++) z[i]=i;
    printf("Input \n");
    for (int i=0; i < size; i++) printf("z[%i] = %lf\t", i, z[i]);
    printf("\n");
    $self->computeFmcp(size,z,F);
    if  (!PyErr_Occurred())
    {
      $self->computeNablaFmcp(size,z,nablaF);
    }
    printf("Output \n");
    for (int i=0; i < size; i++) printf("F[%i] =  %lf\t", i, F[i]); 
    printf("\n");
    for (int i=0; i < size*size; i++) printf("nablaF[%i] =  %lf\t", i, nablaF[i]);
    
    printf("\n");
    free(z);
    free(F);
    free(nablaF);
    
    

    printf("I leave test_call_to_callback()\n");
  }
 


  MixedComplementarityProblem(PyObject *sizeEq, PyObject *sizeIneq, PyObject *o1, PyObject *o2)
  {
     MixedComplementarityProblem* MCP;
     MCP =  (MixedComplementarityProblem *) malloc(sizeof(MixedComplementarityProblem));

     MCP->sizeEqualities = (int) PyInt_AsLong(sizeEq);
     MCP->sizeInequalities = (int) PyInt_AsLong(sizeIneq);
     int size =  MCP->sizeEqualities +  MCP->sizeInequalities;
     
     if (size<1)
     {
       PyErr_SetString(PyExc_RuntimeError, "sizeEqualities + sizeInequalities has to be positive");
       MCP->Fmcp = NULL;
       MCP->nablaFmcp = NULL;
       freeMixedComplementarityProblem(MCP);
       return NULL;
     }
     else
     {
       MCP->Fmcp = (double *) malloc(size*sizeof(double));
       MCP->nablaFmcp = (double *) malloc(size*size*sizeof(double));
     }
     
     if (PyCallable_Check(o1)) 
     {
       set_my_callback_Fmcp(o1);
       MCP->computeFmcp = (my_call_to_callback_Fmcp);
     }
     else
     {
       PyErr_SetString(PyExc_TypeError, "argument 3 must be callable");
       free(MCP->Fmcp);
       free(MCP->nablaFmcp);
       MCP->Fmcp = NULL;
       MCP->nablaFmcp = NULL;
       freeMixedComplementarityProblem(MCP);
       return NULL;
     }

     
     if (PyCallable_Check(o2))
     {
       set_my_callback_NablaFmcp(o2);
       MCP->computeNablaFmcp = (my_call_to_callback_NablaFmcp);
     }
     else
     {
       PyErr_SetString(PyExc_TypeError, "argument 4 must be callable");
       free(MCP->Fmcp);
       free(MCP->nablaFmcp);
       MCP->Fmcp = NULL;
       MCP->nablaFmcp = NULL;
       freeMixedComplementarityProblem(MCP);
       return NULL;
     }

     return MCP;
   }

  ~MixedComplementarityProblem()
  {
    free($self->Fmcp);
    free($self->nablaFmcp);
    $self->Fmcp = NULL;
    $self->nablaFmcp = NULL;
    freeMixedComplementarityProblem($self);
  }
};

%{

  static void call_py_compute_Fmcp(void *env, int n1, int n2, double* z, double* Fmcp)
  {
    npy_intp dim[1];
    dim[0] = n1 + n2;

    PyObject* py_z = FPyArray_SimpleNewFromData(1, dim, NPY_DOUBLE, z);
    PyObject* py_Fmcp = FPyArray_SimpleNewFromData(1, dim, NPY_DOUBLE, Fmcp);

    PyObject* py_n1 = PyInt_FromLong(n1);
    PyObject* py_n2 = PyInt_FromLong(n2);

    PY_CALL_METHOD_OR_FUNCTION(env, "compute_Fmcp", env_compute_function, py_n1, py_n2, py_z, py_Fmcp)

    Py_DECREF(py_n1);
    Py_DECREF(py_n2);

  };

  static void call_py_compute_nabla_Fmcp(void *env, int n1, int n2, double* z, double* nabla_Fmcp)
  {
    npy_intp dim[1];
    dim[0] = n1 + n2;

    npy_intp dim2[2];
    dim2[0] = n1 + n2;
    dim2[1] = n1 + n2;

    PyObject* py_z = FPyArray_SimpleNewFromData(1, dim, NPY_DOUBLE, z);
    PyObject* py_nabla_Fmcp = FPyArray_SimpleNewFromData(2, dim2, NPY_DOUBLE, nabla_Fmcp);

    PyObject* py_n1 = PyInt_FromLong(n1);
    PyObject* py_n2 = PyInt_FromLong(n2);

    PY_CALL_METHOD_OR_FUNCTION(env, "compute_nabla_Fmcp", env_compute_jacobian, py_n1, py_n2, py_z, py_nabla_Fmcp)
    Py_DECREF(py_n1);
    Py_DECREF(py_n2);

  };


%}

%{

  static void call_py_compute_F(void *env, int n, double* z, double* F)
  {
    npy_intp dim[1];
    dim[0] = n;


    PyObject* py_z = FPyArray_SimpleNewFromData(1, dim, NPY_DOUBLE, z);
    PyObject* py_F = FPyArray_SimpleNewFromData(1, dim, NPY_DOUBLE, F);

    PyObject* py_n = PyInt_FromLong(n);

    PY_CALL_METHOD_OR_FUNCTION(env, "compute_Fmcp", env_compute_function, py_n, py_z, py_F)

    Py_DECREF(py_n);
  };

  static void call_py_compute_nabla_F(void *env, int n, double* z, double* nabla_F)
  {
    npy_intp dim[1];
    dim[0] = n;

    npy_intp dim2[2];
    dim2[0] = n;
    dim2[1] = n;


    PyObject* py_z = FPyArray_SimpleNewFromData(1, dim, NPY_DOUBLE, z);
    PyObject* py_nabla_F = FPyArray_SimpleNewFromData(2, dim2, NPY_DOUBLE, nabla_F);

    PyObject* py_n = PyInt_FromLong(n);

    PY_CALL_METHOD_OR_FUNCTION(env, "compute_nabla_Fmcp", env_compute_jacobian, py_n, py_z, py_nabla_F)

    Py_DECREF(py_n);
  };


%}



%extend MixedComplementarityProblem2_
{
  MixedComplementarityProblem2()
   {
     MixedComplementarityProblem2* MCP;
     MCP = (MixedComplementarityProblem2 *) malloc(sizeof(MixedComplementarityProblem2));
     MCP->nabla_Fmcp = NULL;
     MCP->compute_Fmcp = &call_py_compute_Fmcp;
     MCP->compute_nabla_Fmcp = &call_py_compute_nabla_Fmcp;
     MCP->env = NULL;

     return MCP;
   }

  MixedComplementarityProblem2(PyObject* n1, PyObject* n2, PyObject* py_compute_class)
  {
     MixedComplementarityProblem2* MCP;
     MCP =  (MixedComplementarityProblem2 *) malloc(sizeof(MixedComplementarityProblem2));

     MCP->compute_Fmcp = &call_py_compute_Fmcp;
     MCP->compute_nabla_Fmcp = &call_py_compute_nabla_Fmcp;
     MCP->n1 = (int) PyInt_AsLong(n1);
     MCP->n2 = (int) PyInt_AsLong(n2);
     int size =  MCP->n1 +  MCP->n2;

     if (size<1)
     {
       PyErr_SetString(PyExc_RuntimeError, "sizeEqualities + sizeInequalities has to be positive");
       free(MCP);
       return NULL;
     }
     else
     {
       MCP->nabla_Fmcp = (double *) calloc(size*size, sizeof(double));
     }

     PyObject* method_compute_Fmcp = PyObject_GetAttrString(py_compute_class, "compute_Fmcp");
     PyObject* method_compute_nabla_Fmcp = PyObject_GetAttrString(py_compute_class, "compute_nabla_Fmcp");
     if (PyCallable_Check(method_compute_Fmcp) && PyCallable_Check(method_compute_nabla_Fmcp))
     {
       MCP->env = (void*) malloc(sizeof(class_env_python));
       class_env_python* mcp_env_python = (class_env_python*) MCP->env;
       mcp_env_python->id = ENV_IS_PYTHON_CLASS;
       mcp_env_python->class_object = py_compute_class;
       Py_DECREF(method_compute_Fmcp);
       Py_DECREF(method_compute_nabla_Fmcp);
     }
     else
     {
       Py_XDECREF(method_compute_Fmcp);
       Py_XDECREF(method_compute_nabla_Fmcp);
       PyErr_SetString(PyExc_TypeError, "argument 2 must be have a method compute_Fmcp and a method compute_nabla_Fmcp");
       free(MCP->nabla_Fmcp);
       free(MCP);
       PyErr_PrintEx(0);
       exit(1);
     }

     return MCP;
   }

  MixedComplementarityProblem2(PyObject* n1, PyObject* n2, PyObject* py_compute_Fmcp, PyObject* py_compute_nabla_Fmcp)
  {
     MixedComplementarityProblem2* MCP;
     MCP =  (MixedComplementarityProblem2 *) malloc(sizeof(MixedComplementarityProblem2));

     MCP->compute_Fmcp = &call_py_compute_Fmcp;
     MCP->compute_nabla_Fmcp = &call_py_compute_nabla_Fmcp;
     MCP->n1 = (int) PyInt_AsLong(n1);
     MCP->n2 = (int) PyInt_AsLong(n2);
     int size =  MCP->n1 +  MCP->n2;

     if (size<1)
     {
       PyErr_SetString(PyExc_RuntimeError, "sizeEqualities + sizeInequalities has to be positive");
       free(MCP);
       return NULL;
     }
     else
     {
       MCP->nabla_Fmcp = (double *) malloc(size*size*sizeof(double));
     }

     MCP->env = (void*) malloc(sizeof(functions_env_python));
     functions_env_python* mcp_env_python = (functions_env_python*) MCP->env;
     mcp_env_python->id = ENV_IS_PYTHON_FUNCTIONS;

     if (PyCallable_Check(py_compute_Fmcp)) 
     {
       mcp_env_python->env_compute_function = py_compute_Fmcp;
     }
     else
     {
       PyErr_SetString(PyExc_TypeError, "argument 3 must be callable");
       free(MCP->nabla_Fmcp);
       free(MCP->env);
       free(MCP);
       PyErr_PrintEx(0);
       exit(1);
     }


     if (PyCallable_Check(py_compute_nabla_Fmcp))
     {
       mcp_env_python->env_compute_jacobian = py_compute_nabla_Fmcp;
     }
     else
     {
       PyErr_SetString(PyExc_TypeError, "argument 4 must be callable");
       free(MCP->nabla_Fmcp);
       free(MCP->env);
       free(MCP);
       PyErr_PrintEx(0);
       exit(1);
     }

     return MCP;
   }
  ~MixedComplementarityProblem2()
  {
    free($self->nabla_Fmcp);
    free($self->env);
    free($self);
  }
};

%include "VariationalInequality.h"
%include "VariationalInequality_Solvers.h"
%include "VI_cst.h"

%extend VariationalInequality
{
  VariationalInequality()
   {
     VariationalInequality* vi;
     vi = (VariationalInequality *) malloc(sizeof(VariationalInequality));
     vi->nabla_F = NULL;
     vi->F = &call_py_compute_F;
     vi->compute_nabla_F = &call_py_compute_nabla_F;
     vi->env = NULL;

     return vi;
   }

  VariationalInequality(PyObject* n, PyObject* py_compute)
  {
     VariationalInequality* vi;
     vi =  (VariationalInequality *) malloc(sizeof(VariationalInequality));

     vi->size = (int) PyInt_AsLong(n);
     vi->set = NULL;
     vi->env = NULL;
     vi->compute_nabla_F = NULL;
     vi->ProjectionOnX = NULL;
     vi->nabla_F = NULL;
     vi->normVI = std::numeric_limits<double>::quiet_NaN();
     vi->istheNormVIset = 0;


     if (vi->size < 1)
     {
       PyErr_SetString(PyExc_RuntimeError, "the size of the VI has to be positive");
       free(vi);
       return NULL;
     }


     PyObject* method_compute_F = PyObject_GetAttrString(py_compute, "compute_F");
     PyObject* method_compute_nabla_F = PyObject_GetAttrString(py_compute, "compute_nabla_F");

     if (PyCallable_Check(method_compute_F) && PyCallable_Check(method_compute_nabla_F))
     {
       vi->env = (void*) malloc(sizeof(class_env_python));
       class_env_python* vi_env_python = (class_env_python*) vi->env;
       vi_env_python->id = ENV_IS_PYTHON_CLASS;
       vi_env_python->class_object = py_compute;
       Py_DECREF(method_compute_F);
       Py_DECREF(method_compute_nabla_F);
     }
     else
     {
       if (PyCallable_Check(py_compute))
       {
         vi->F = &call_py_compute_F;
         vi->env = (void*) malloc(sizeof(functions_env_python));
         functions_env_python* vi_env_python = (functions_env_python*) vi->env;
         vi_env_python->id = -1;
         vi_env_python->env_compute_function = py_compute;
       }
       else
       {
         Py_XDECREF(method_compute_F);
         Py_XDECREF(method_compute_nabla_F);
         PyErr_SetString(PyExc_TypeError, "argument 2 must either be an object with a method compute_F and a method compute_nabla_F or a callable function");
         free(vi->nabla_F);
         free(vi);
         PyErr_PrintEx(0);
         exit(1);
       }
     }

     return vi;
   }


   void set_compute_nabla_F(PyObject* py_compute_nabla_F)
   {
     if (PyCallable_Check(py_compute_nabla_F))
     {
       $self->compute_nabla_F = &call_py_compute_nabla_F;
       functions_env_python* vi_env_python = (functions_env_python*) $self->env;
       vi_env_python->id = ENV_IS_PYTHON_FUNCTIONS;
       vi_env_python->env_compute_jacobian = py_compute_nabla_F;
       $self->nabla_F = (double*) calloc($self->size*$self->size, sizeof(double));
     }
     else
     {
       PyErr_SetString(PyExc_TypeError, "argument 1 must be callable");
       PyErr_PrintEx(0);
     }
   }

   void set_box_constraints(PyObject* box_lower_bound, PyObject* box_upper_bound)
   {
     if ((PyObject_Length(box_lower_bound) == $self->size) && (PyObject_Length(box_upper_bound) == $self->size))
     {
       $self->set = (void *)malloc(sizeof(box_constraints));
       box_constraints* box_c = (box_constraints*)$self->set;
       box_c->id = VI_SET_IS_BOX;
       box_c->lb = (double*)malloc($self->size*sizeof(double));
       box_c->ub = (double*)malloc($self->size*sizeof(double));

       if (!convert_darray(box_lower_bound, box_c->lb) || !convert_darray(box_upper_bound, box_c->ub))
       {
         PyErr_PrintEx(0);
         exit(1);
       }
     }
     else
     {
       PyErr_SetString(PyExc_TypeError, "The arguments do not have the right length");
       PyErr_PrintEx(0);
       exit(1);
     }
   }



  ~VariationalInequality()
  {
    if ($self->nabla_F)
    {
      free($self->nabla_F);
    }
    if ($self->set)
    {
      //black magic
      switch (((generic_set*)$self->set)->id)
      {
        case VI_SET_IS_BOX:
        {
          box_constraints* box_c = (box_constraints*)$self->set;
          free_box(box_c);
          break;
        }
        default:
          PyErr_SetString(PyExc_TypeError, "unknown set type");
      }
      free($self->set);
      $self->set = NULL;
    }
    free($self);
  }
};


%typemap(out) (double* q) {
  npy_intp dims[2];

  dims[0] = arg1->dimension * arg1->numberOfContacts;
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


%typemap(out) (double* mu) {
  npy_intp dims[2];

  dims[0] = arg1->numberOfContacts;
  dims[1] = 1;
  if (arg1->q)
  {
    PyObject *obj = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, arg1->mu);
    PyArrayObject *array = (PyArrayObject*) obj;
    if (!array || !require_fortran(array)) SWIG_fail;
    $result = obj;
  }
  else
    SWIG_fail;
 }





// FrictionContact
%ignore LocalNonsmoothNewtonSolver; //signature problem (should be SolverOption
                          //instead of *iparam, *dparam).
%ignore DampedLocalNonsmoothNewtonSolver; // signature problem idem.

%ignore frictionContactProblem_new; // signature issue with mu param


%include FrictionContactProblem.i

%include "FrictionContact3D_Solvers.h"
%include "FrictionContact2D_Solvers.h"
%include "Friction_cst.h"
%include "FrictionContact3D_AlartCurnier.h"
%include "FrictionContact3D_localAlartCurnier.h"
%include "FrictionContact3D_compute_error.h"
%include "fclib_interface.h"

%extend FrictionContactProblem
{
  FrictionContactProblem()
  {
    FrictionContactProblem * FCP = (FrictionContactProblem *) malloc(sizeof(FrictionContactProblem));
    
    return FCP;
  }


  /* copy constructor */
  FrictionContactProblem(PyObject *o)
  {
    void * vp;
    FrictionContactProblem* fcp;
    FrictionContactProblem* FCP;
    int res = SWIG_ConvertPtr(o, &vp,SWIGTYPE_p_FrictionContactProblem, 0 |  0 );
    if (!SWIG_IsOK(res)) return 0;
    fcp = (FrictionContactProblem *) vp;
    FCP = (FrictionContactProblem *) malloc(sizeof(FrictionContactProblem));
    FCP->dimension = fcp->dimension;
    FCP->M = fcp->M;
    FCP->numberOfContacts = fcp->numberOfContacts;
    FCP->q =  fcp->q;
    FCP->mu = fcp->mu;

    Py_INCREF(o);

    return FCP;
  }

  /* */
  FrictionContactProblem(PyObject *dim, PyObject *numberOfContacts, PyObject *M, PyObject *q, PyObject *mu)
  {
    FrictionContactProblem * FC = (FrictionContactProblem *) malloc(sizeof(FrictionContactProblem));
    FC->dimension = PyInt_AsLong(dim);
    FC->numberOfContacts = PyInt_AsLong(numberOfContacts);
    
    {
      void * _M;
      int res = SWIG_ConvertPtr(M, &_M,SWIGTYPE_p_NumericsMatrix, 0 |  0 );
      if (!SWIG_IsOK(res)) return 0;
      FC->M = (NumericsMatrix *) _M;
    }
    {
      void * _q;
      int res = SWIG_ConvertPtr(M, &_q,SWIGTYPE_p_double, 0 |  0 );
      if (!SWIG_IsOK(res)) return 0;
      FC->q = (double *) _q;
    }
    {
      void * _mu;
      int res = SWIG_ConvertPtr(M, &_mu,SWIGTYPE_p_double, 0 |  0 );
      if (!SWIG_IsOK(res)) return 0;
      FC->mu = (double *) _mu;
    }

    return FC;
    
  }

  NumericsMatrix* rawM()
  {
    return $self->M;
  }

  FrictionContactProblem(PyObject *dim, PyObject *o1, PyObject *o2, PyObject *o3)
    {

      int is_new_object1=0;
      int is_new_object2=0;
      int is_new_object3=0; 

      PyArrayObject* array = obj_to_array_fortran_allow_conversion(o1, NPY_DOUBLE,&is_new_object1);
      PyArrayObject* vector = obj_to_array_contiguous_allow_conversion(o2, NPY_DOUBLE, &is_new_object2);
      PyArrayObject* mu_vector = obj_to_array_contiguous_allow_conversion(o3, NPY_DOUBLE, &is_new_object3); 
      FrictionContactProblem *FC;
      // return pointer : free by std swig destructor
      FC = (FrictionContactProblem *) malloc(sizeof(FrictionContactProblem));
      NumericsMatrix *M = (NumericsMatrix *) malloc(sizeof(NumericsMatrix));

      M->storageType = 0;
      M->size0 = array_size(array,0);
      M->size1 = array_size(array,1);
      M->matrix0 = (double *) malloc(M->size0*M->size1*sizeof(double));
      memcpy(M->matrix0,array_data(array),M->size0*M->size1*sizeof(double));
      FC->dimension = (int) PyInt_AsLong(dim);
      FC->numberOfContacts = M->size0 / FC->dimension;
      FC->M = M;
      FC->q = (double *) malloc(M->size0*sizeof(double));
      memcpy(FC->q,array_data(vector),M->size0*sizeof(double));
      FC->mu = (double *) malloc(FC->numberOfContacts*sizeof(double));
      memcpy(FC->mu,array_data(mu_vector),FC->numberOfContacts*sizeof(double));
   

      // python mem management
      if(is_new_object1 && array)
      {
        Py_DECREF(array);
      }
  
      if(is_new_object2 && vector)
      {
        Py_DECREF(vector);
      }

      if(is_new_object3 && mu_vector)
      {
        Py_DECREF(mu_vector);
      }

      return FC;
    }

  ~FrictionContactProblem()
  {
    freeFrictionContactProblem($self);
  }

};


%extend GlobalFrictionContactProblem
{

  GlobalFrictionContactProblem(PyObject *dim, PyObject *o1, PyObject *o2, PyObject *o3)
    {

      int is_new_object1=0;
      int is_new_object2=0;
      int is_new_object3=0; 

      PyArrayObject* array = obj_to_array_fortran_allow_conversion(o1, NPY_DOUBLE,&is_new_object1);
      PyArrayObject* vector = obj_to_array_contiguous_allow_conversion(o2, NPY_DOUBLE, &is_new_object2);
      PyArrayObject* mu_vector = obj_to_array_contiguous_allow_conversion(o3, NPY_DOUBLE, &is_new_object3); 
      GlobalFrictionContactProblem *FC;
      // return pointer : free by std swig destructor
      FC = (GlobalFrictionContactProblem *) malloc(sizeof(GlobalFrictionContactProblem));
      NumericsMatrix *M = (NumericsMatrix *) malloc(sizeof(NumericsMatrix));

      M->storageType = 0;
      M->size0 = array_size(array,0);
      M->size1 = array_size(array,1);
      M->matrix0 = (double *) malloc(M->size0*M->size1*sizeof(double));
      memcpy(M->matrix0,array_data(array),M->size0*M->size1*sizeof(double));
      FC->dimension = (int) PyInt_AsLong(dim);
      FC->numberOfContacts = M->size0 / FC->dimension;
      FC->M = M;
      FC->q = (double *) malloc(M->size0*sizeof(double));
      memcpy(FC->q,array_data(vector),M->size0*sizeof(double));
      FC->mu = (double *) malloc(FC->numberOfContacts*sizeof(double));
      memcpy(FC->mu,array_data(mu_vector),FC->numberOfContacts*sizeof(double));
   

      // python mem management
      if(is_new_object1 && array)
      {
        Py_DECREF(array);
      }
  
      if(is_new_object2 && vector)
      {
        Py_DECREF(vector);
      }

      if(is_new_object3 && mu_vector)
      {
        Py_DECREF(mu_vector);
      }

      return FC;
    }

  ~GlobalFrictionContactProblem()
  {
    freeGlobalFrictionContactProblem($self);
  }

};

//GenericMechanical
//%include GMPReduced.h
//%include GenericMechanicalProblem.h
//%include GenericMechanical_Solvers.h
%{
#include <GenericMechanical_cst.h>
%}
%include GenericMechanical_cst.h


 // Matrices
%include "SparseMatrix.h"
%include "SparseBlockMatrix.h"
%include "NumericsMatrix.h"

// some extensions but numpy arrays should be used instead 
%extend NumericsMatrix
{

  NumericsMatrix(int nrows, int ncols, double *data)
    {
      NumericsMatrix *M;
      
      // return pointer : free by std swig destructor
      M = (NumericsMatrix *) malloc(sizeof(NumericsMatrix));
      M->storageType = 0;
      M->size0 = nrows;
      M->size1 = ncols;
      M->matrix0 = data;
      return M;
    }

  void set_matrix0(int i, int j, double v)
  {
    self->matrix0[i+j*self->size1] = v;
  }

  double get_matrix0(int i, int j)
  {
    return self->matrix0[i+j*self->size1];
  }


  PyObject * __setitem__(PyObject* index, double v)
  {
    int i, j;
    if (!PyArg_ParseTuple(index, "ii:NumericsMatrix___setitem__",&i,&j)) return NULL;
    NumericsMatrix_set_matrix0(self,i,j,v);
    return Py_BuildValue("");
  }

  PyObject * __getitem__(PyObject * index)
  {
    int i, j;
    if (!PyArg_ParseTuple(index, "ii:NumericsMatrix___getitem__",&i,&j)) return NULL;
    return SWIG_From_double(NumericsMatrix_get_matrix0(self,i,j));
  }

  int __len__()
  {
    return self->size0 * self->size1;
  }

  PyObject * __str__()
  {
    std::stringstream result;
    result << "[ ";
    for (int i=0; i < self->size0; ++i)
      {
        if (i > 0) result << "  ";
        result << "[";
        for (int j=0; j < self->size1; ++j)
          {
            result << " " << NumericsMatrix_get_matrix0(self,i,j);
            if (j < self->size1-1) result << ",";
          }
        result << " ]";
        if (i < self->size0-1) result << "," << std::endl;
      }
    result << " ]" << std::endl;
    %#if PY_MAJOR_VERSION < 3
    return PyString_FromString(result.str().c_str());
    %#else
    return PyUnicode_FromString(result.str().c_str());
    %#endif
  }
  
}; 


#define GET_ATTR(OBJ,ATTR)                                              \
  ATTR = PyObject_GetAttrString(OBJ,BOOST_PP_STRINGIZE(ATTR));          \
  if(!ATTR)                                                             \
  {                                                                     \
    throw(std::invalid_argument(BOOST_PP_STRINGIZE(need a ATTR attr))); \
  }

typedef struct cs_sparse    /* matrix in compressed-column or triplet form */
{
  int nzmax ;	    /* maximum number of entries */
  int m ;	    /* number of rows */
  int n ;	    /* number of columns */
  int *p ;	    /* column pointers (size n+1) or col indices (size nzmax) */
  int *i ;	    /* row indices, size nzmax */
  double *x ;	    /* numerical values, size nzmax */
  int nz ;	    /* # of entries in triplet matrix, -1 for compressed-col */
} cs ;

%extend cs_sparse
{
  %fragment("NumPy_Fragments");

  // from a scipy.sparse csr
  cs_sparse(PyObject *obj)
  {

    PyObject *shape,*nnz,*data,*indices,*indptr;
    int is_new_object1=0;
    int is_new_object2=0;
    int is_new_object3=0;
    PyArrayObject *array_data, *array_indices, *array_indptr;
    cs_sparse* M;

    try
    {

      M = (cs_sparse *) malloc(sizeof(cs_sparse));      
      
      GET_ATTR(obj,shape);
      GET_ATTR(obj,nnz);
      GET_ATTR(obj,data);
      GET_ATTR(obj,indices);
      GET_ATTR(obj,indptr);
      
      int dim0, dim1, nzmax;

      GET_INTS(shape,0,dim0);
      GET_INTS(shape,1,dim1);
//      GET_INT(nnz,nzmax); fail: type is numpy.int32!
      nzmax = PyInt_AsLong(nnz);

      array_data = obj_to_array_allow_conversion(data, NPY_DOUBLE, &is_new_object1);
      array_indices = obj_to_array_allow_conversion(indices, NPY_INT32, &is_new_object2);
      array_indptr = obj_to_array_allow_conversion(indptr, NPY_INT32, &is_new_object3);
      
      M->m = dim0;
      M->n = dim1;

      M->nzmax = nzmax;

      M->nz = -2; // csr only for the moment
      
      M->p = (int *) malloc((M->m+1) * sizeof(int));
      M->i = (int *) malloc(M->nzmax * sizeof(int));
      M->x = (double *) malloc(M->nzmax * sizeof(double));

      for(unsigned int i = 0; i < (M->m+1); i++)
      {
        M->p[i] = ((int *) array_data(array_indptr)) [i];
      }

      for(unsigned int i = 0; i< M->nzmax; i++)
      {
        M->i[i] = ((int *) array_data(array_indices)) [i];
      }
     
      memcpy(M->x, (double *) array_data(array_data), M->nzmax * sizeof(double));

      Py_DECREF(shape);
      Py_DECREF(nnz);
      Py_DECREF(data);
      Py_DECREF(indices);
      Py_DECREF(indptr);

      if (array_data && is_new_object1)
      {
        Py_DECREF(array_data);
      }

      if (array_indptr && is_new_object2)
      {
        Py_DECREF(array_indptr);
      }

      if (array_indices && is_new_object3)
      {
        Py_DECREF(array_indices);
      }
      

      return M;
    }
    catch (const std::invalid_argument& e)
    {
      Py_XDECREF(shape);
      Py_XDECREF(nnz);
      Py_XDECREF(data);
      Py_XDECREF(indices);
      Py_XDECREF(indptr);

      if (array_data && is_new_object1)
      {
        Py_DECREF(array_data);
      }
      
      if (array_indptr && is_new_object2)
      {
        Py_DECREF(array_indptr);
      }
      
      if (array_indices && is_new_object3)
      {
        Py_DECREF(array_indices);
      }
      
      assert(!M->p);
      assert(!M->i);
      assert(!M->x);
      assert(M);
      free(M);
      throw(e);
    }
  }

  ~cs_sparse()
  {
    freeSparse($self);
  }
}


%{

 
%}

#ifdef WITH_IO
%include picklable.i

%make_picklable(Callback, Numerics);
%make_picklable(_SolverOptions, Numerics);
%make_picklable(SolverOptions, Numerics);
%make_picklable(FrictionContactProblem, Numerics);
%make_picklable(NumericsMatrix, Numerics);
%make_picklable(SparseBlockStructuredMatrix, Numerics);
#endif
