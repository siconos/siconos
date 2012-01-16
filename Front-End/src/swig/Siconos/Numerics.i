// -*- C++ -*-
// Siconos-Front-End , Copyright INRIA 2005-2011.
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
%module Numerics

%{
#define SWIG_FILE_WITH_INIT
#include <sstream>
#if defined(Py_COMPLEXOBJECT_H)
#undef c_sum
#undef c_diff
#undef c_neg
#undef c_prod
#undef c_quot
#undef c_pow
#undef c_abs
#endif
#include "SiconosNumerics.h"
#include "NumericsConfig.h"
#include "SolverOptions.h"
#include "SparseBlockMatrix.h"
#include "FrictionContactProblem.h"
#include "FrictionContact3D_Solvers.h"
#include "Friction_cst.h"
#include "frictionContact_test_function.h"
#include "FrictionContact3D_AlartCurnier.h"
#include "FrictionContact3D_globalAlartCurnier.h"
#include "FrictionContact3D_compute_error.h"
#include "fclib_interface.h"

#include <boost/preprocessor/stringize.hpp>
%}

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
  }
%}
  
// needed macros
%include "NumericsConfig.h"


// more convenient
%rename (LCP) LinearComplementarityProblem;

%ignore lcp_compute_error_only;

// numpy macros
%include numpy.i 	

%init %{
  import_array();
%}

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

%include "Common.i"



%typemap(in) (LinearComplementarityProblem*) (npy_intp problem_size) {
  void *lcp;
  int res = SWIG_ConvertPtr($input, &lcp,SWIGTYPE_p_LinearComplementarityProblem, 0 |  0 );
  if (!SWIG_IsOK(res)) SWIG_fail;

  problem_size=((LinearComplementarityProblem *) lcp)->size;

  $1 = (LinearComplementarityProblem *) lcp;
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

// vectors of size problem_size from given *Problem as first input
// no conversion => inout array XXX FIX issue here
%typemap(in) (double *z) (PyArrayObject* array=NULL, int is_new_object = 0) {

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
%typemap(in) (double *blocklist3x3) (PyArrayObject* array=NULL, int is_new_object) {

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
%typemap(in) (double *blockarray3x3) (PyArrayObject* array=NULL, int is_new_object) {

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
      *p_problem_size1 = sqrt(array_size(array,0));
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
%typemap(in) (double *blocklist3) (PyArrayObject* array=NULL, int is_new_object) {

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
%typemap(in) (double *mu) (PyArrayObject* array=NULL, int is_new_object) {

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

  if (nummat->matrix0) free(nummat->matrix0);

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
  free($1);
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
%typemap(in) (FILE *file)
{
  // %typemap(in) (FILE *file)
  $1 = fopen(PyString_AsString($input), "r");
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
  fclose($1);
}

%typemap(in, numinputs=0) (SparseBlockStructuredMatrix* const M) 
{
  $1 = (SparseBlockStructuredMatrix*) malloc(sizeof(SparseBlockStructuredMatrix));
}

%typemap(argout) (SparseBlockStructuredMatrix* const M)
{
  $result = SWIG_Python_AppendOutput($result,
                                     SWIG_NewPointerObj(SWIG_as_voidptr($1), 
                                                        SWIGTYPE_p_SparseBlockStructuredMatrix, 0));
}

// signatures
%feature("autodoc", 1);

// generated docstrings from doxygen xml output
%include Numerics-docstrings.i
 

// LCP
%include "SparseBlockMatrix.h"
%include "NumericsMatrix.h"
%include "LinearComplementarityProblem.h"
%include "LCP_Solvers.h"
%include "lcp_cst.h"
%include "SolverOptions.h"
%include "NumericsOptions.h"
%include "frictionContact_test_function.h"

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

  SolverOptions* makeSolverOptions(FrictionContactProblem* fcp, FRICTION_SOLVER id)
  {
    SolverOptions *SO;
    SO = (SolverOptions *) malloc(sizeof(SolverOptions));
    switch(id)
    {
    case SICONOS_FRICTION_2D_NSGS:
    case SICONOS_FRICTION_2D_NLGS:
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

  SolverOptions(FRICTION_SOLVER id)
  {
    return SolverOptions_makeSolverOptions(NULL, NULL, id);
  }

  SolverOptions(FrictionContactProblem* fcp, FRICTION_SOLVER id)
  {
    return SolverOptions_makeSolverOptions(NULL, fcp, id);
  }

  ~SolverOptions() 
    { 
      deleteSolverOptions(self);
    }
};

%extend LinearComplementarityProblem
{
  LinearComplementarityProblem(PyObject *o1, PyObject *o2)
    {

      int is_new_object1, is_new_object2;
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

%include "FrictionContactProblem.h"
%include "FrictionContact3D_Solvers.h"
%include "Friction_cst.h"
%include "FrictionContact3D_AlartCurnier.h"
%include "FrictionContact3D_globalAlartCurnier.h"
%include "FrictionContact3D_compute_error.h"
%include "fclib_interface.h"

%extend FrictionContactProblem
{

  FrictionContactProblem(PyObject *dim, PyObject *o1, PyObject *o2, PyObject *o3)
    {

      int is_new_object1, is_new_object2, is_new_object3; 

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


// some extensions but numpy arrays should be used instead 
%extend NumericsMatrix
{

  NumericsMatrix(int nrows, int ncols, double *data)
    {
      NumericsMatrix *M;
      
      // return pointer : free by std swig destructor
      M = (NumericsMatrix *) malloc(sizeof(NumericsMatrix));
      M->storageType=0;
      M->size0 = nrows;
      M->size1 = ncols;
      M->matrix0=data;
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
    return PyString_FromString(result.str().c_str());
  }
  
}; 


