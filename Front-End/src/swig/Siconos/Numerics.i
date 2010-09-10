// -*- C++ -*-
// Siconos-Front-End version 3.2.0, Copyright INRIA 2005-2010.
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
#include "SolverOptions.h"
#include "frictionContact_test_function.h"
%}

%inline %{
#include <stdio.h>
  FrictionContactProblem* frictionContactProblemFromFile(const char * filename)
  {
    FILE * finput = fopen(filename, "r");
    FrictionContactProblem* problem = (FrictionContactProblem *) malloc(sizeof(FrictionContactProblem));
    if (frictionContact_newFromFile(problem,finput))
    {
      fprintf(stderr, "frictionContactProblemFromFile: cannot load %s\n",filename);
      return problem;
    }
    else
      return problem;
  }

%}




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

// std python sequence -> C array

// int
%{
  static int convert_iarray(PyObject *input, int *ptr) {
  int i;
  if (!PySequence_Check(input)) {
      PyErr_SetString(PyExc_TypeError,"Expecting a sequence");
      return NULL;
  }
  for (i =0; i <  PyObject_Length(input); i++) {
      PyObject *o = PySequence_GetItem(input,i);
      if (!PyInt_Check(o)) {
        Py_XDECREF(o);
        PyErr_SetString(PyExc_ValueError,"Expecting a sequence of ints");
        return NULL;
      }
      ptr[i] = (int) PyInt_AsLong(o);

      if (ptr[i] == -1 && PyErr_Occurred())
        return NULL;
      Py_DECREF(o);
  }
  return 1;
}
%}

// double
%{
static int convert_darray(PyObject *input, double *ptr) {
  int i;
  if (!PySequence_Check(input)) {
      PyErr_SetString(PyExc_TypeError,"Expecting a sequence");
      return NULL;
  }

  for (i =0; i < PyObject_Length(input); i++) {
      PyObject *o = PySequence_GetItem(input,i);
      if (!PyFloat_Check(o)) {
         Py_XDECREF(o);
         PyErr_SetString(PyExc_ValueError,"Expecting a sequence of floats");
         return NULL;
      }
      ptr[i] = PyFloat_AsDouble(o);

      if (ptr[i] == -1 && PyErr_Occurred())
        return NULL;

      Py_DECREF(o);
  }
  return 1;
  
}
%}

// conversions
%typemap(in) (int sizei, int *iparam) (int *temp) {
  
  temp = (int *) malloc(sizeof(int)*PyObject_Length($input));
    if (!convert_iarray($input,temp)) {
      SWIG_fail;
  }
  $1 = PyObject_Length($input);
  $2 = &temp[0];

  }

%typemap(in) (int sized, double *dparam) (double *temp) {

  temp = (double *) malloc(sizeof(double)*PyObject_Length($input));
    if (!convert_darray($input,temp)) {
      SWIG_fail;
    }
    $1 = PyObject_Length($input);
    $2 = &temp[0];
 }


// cleanup
%typemap(freearg) (int sized, double *dparam) {
  free($2);
 }

%typemap(freearg) (int sizei, int *iparam) {
  free($2);
 }

// info result
%typemap(in, numinputs=0) int *info (int temp) {
  $1 = &temp;
}

%typemap(argout) (int *info) {
  $result = PyInt_FromLong(*$1);
 }

%typemap(in) (LinearComplementarityProblem*) (npy_intp problem_size) {
  void *lcp;
  int res = SWIG_ConvertPtr($input, &lcp,SWIGTYPE_p_LinearComplementarityProblem, 0 |  0 );
  if (!SWIG_IsOK(res)) SWIG_fail;

  problem_size=((LinearComplementarityProblem *) lcp)->size;

  $1 = (LinearComplementarityProblem *) lcp;
}

%typemap(in, numinputs=0) (unsigned int problemSize) (npy_intp problem_size, npy_intp number_of_contacts)
{
  // the first array length sets problemSize
  problem_size = -1;
  number_of_contacts = -1;
}

%typemap(in) (FrictionContactProblem*) (npy_intp problem_size, npy_intp problem_dimension, npy_intp number_of_contacts) {
  void *fcp;
  int res = SWIG_ConvertPtr($input, &fcp,SWIGTYPE_p_FrictionContactProblem, 0 |  0 );
  if (!SWIG_IsOK(res)) SWIG_fail;

  problem_dimension=((FrictionContactProblem *) fcp)->dimension;
  number_of_contacts=((FrictionContactProblem *) fcp)->numberOfContacts;
  problem_size=((FrictionContactProblem *) fcp)->numberOfContacts * problem_dimension;

  $1 = (FrictionContactProblem*) fcp;
}

// vectors of size problem_size from given *Problem as first input
%typemap(in) (double *z) (PyArrayObject* array=NULL, int is_new_object) {

  array = obj_to_array_contiguous_allow_conversion($input, NPY_DOUBLE,&is_new_object);
  
  if (!array || !require_dimensions(array,1) ||
      !require_native(array) || !require_contiguous(array)
      || !require_size(array, &problem_size1, 1)) SWIG_fail;
  
  $1 = (double *) array_data(array);

 }

// vectors of size problem_size from problemSize as first input
%typemap(in) (double *vect3D) (PyArrayObject* array=NULL, int is_new_object) {

  array = obj_to_array_contiguous_allow_conversion($input, NPY_DOUBLE,&is_new_object);

  if (problem_size1 < 0)
  {
    if (array_numdims(array) > 0)
    {
      problem_size1 = array_size(array,0);
      arg1 = problem_size1;
    }
    
    if (problem_size1 % 3 != 0) SWIG_fail;
    
    if (problem_size1 / 3 == 0) SWIG_fail;
    
    
    number_of_contacts1 = problem_size1 / 3;
    
  }
  
  if (!array || !require_dimensions(array,1) ||
      !require_native(array) || !require_contiguous(array)
      || !require_size(array, &problem_size1, 1)) SWIG_fail;
  
  $1 = (double *) array_data(array);

 }



// vectors of size problem_size

// 1 : numinputs=0 mandatory to avoid arg
%typemap(in, numinputs=0) (double *result) (PyObject* array=NULL)
{
    // %typemap(in, numinputs=0)
    // we cannot get problem_size here as numinputs=0 => before
    // numinputs=1, how can we change this ??
}

// 2 : check must be done after in
%typemap(check) (double *result) 
{
  if (problem_size1 > 0)
  {
    
    npy_intp dims[1] = { problem_size1 };
    
    array$argnum = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (!array$argnum) SWIG_fail;
    $1 = ($1_ltype) array_data(array$argnum);
  }
  
}

// 3 : return arg
%typemap(argout) (double *result)
{
  if (problem_size1 > 0)
  {
     $result = SWIG_Python_AppendOutput($result,(PyObject *)array$argnum);
  }
  
}


// 3x3 matrices

// 1 : numinputs=0 mandatory to avoid arg
%typemap(in, numinputs=0) (double *result3x3) (PyObject* array=NULL)
{
    // %typemap(in, numinputs=0)
    // we cannot get problem_size here as numinputs=0 => before
    // numinputs=1, how can we change this ??
}

// 2 : check must be done after in
%typemap(check) (double *result3x3) 
{
  if (problem_size1 > 0)
  {
    
    npy_intp dims[1] = { problem_size1 * 3 };
    
    array$argnum = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (!array$argnum) SWIG_fail;
    $1 = ($1_ltype) array_data(array$argnum);
  }
  
}

// 3 : return arg
%typemap(argout) (double *result3x3)
{
  if (problem_size1 > 0)
  {
      $result = SWIG_Python_AppendOutput($result,(PyObject *)array$argnum);
  }
  
}


// vectors of size numberOfContacts
%typemap(in) (double *mu) (PyArrayObject* array=NULL, int is_new_object) {

  if (number_of_contacts1 > 0)

  {
    array = obj_to_array_contiguous_allow_conversion($input, NPY_DOUBLE,&is_new_object);
    
    if (!array || !require_dimensions(array,1) ||
        !require_native(array) || !require_contiguous(array)
        || !require_size(array, &number_of_contacts1, 1)) SWIG_fail;
    
    $1 = (double *) array_data(array);
    
  }
  
    
 }


// other names that must be transformed this way
%apply (double *z) { (double *w) }; 

%apply (double *z) { (double *zlem) };

%apply (double *z) { (double *wlem) };

%apply (double *z) { (double *reaction) };

%apply (double *z) { (double *velocity) };

%apply (double *vect3D) { (double *reaction3D) };

%apply (double *vect3D) { (double *velocity3D) };

%apply (double *vect3D) { (double *rho3D) };

%apply (double *result) { (double *result1) };

%apply (double *result) { (double *result2) };

%apply (double *result) { (double *result3) };

%apply (double *result) { (double *result4) };

%apply (double *result3x3) { (double *result3x3_1) };

%apply (double *result3x3) { (double *result3x3_2) };

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

%typemap(freearg) (NumericsMatrix* A) {
    free($1);
}

%typemap(out) (NumericsMatrix* M) {
  npy_intp dims[2];
  dims[0] = $1->size0;
  dims[1] = $1->size1;
  if ($1->matrix0)
  {
    $result = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, $1->matrix0);
  }
  else
    SWIG_fail;
 }

%typemap(out) (double* q) {
  npy_intp dims[1];

  dims[0] = arg1->size;
  if (arg1->q)
  {
    $result = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, arg1->q);
  }
  else
    SWIG_fail;
 }

%apply (NumericsMatrix *A) { (NumericsMatrix *m) };
%apply (NumericsMatrix *A) { (NumericsMatrix *M) };

%feature("autodoc", "1");

%typemap(in) (int iSize) ($1_type iparam_size) {
  iparam_size = PyInt_AsLong($input);
 }

%typemap(in) (int dSize) ($1_type dparam_size) {
  dparam_size = PyInt_AsLong($input);
 }

%typemap(in) (int *iparam) {
  
  $1_type temp;
  temp = ($1_type) malloc(sizeof($*1_type)*PyObject_Length($input));
  if (!convert_iarray($input,temp)) {
    SWIG_fail;
  }

  if ($1) { free($1); };

  $1 = &temp[0];

  // arg1 is *SolverOptions. May be version dependent, how to get
  // this?
  if (arg1) arg1->iSize = PyObject_Length($input);

 }

%typemap(in) (double *dparam) {
  
  $1_type temp;
  temp = ($1_type) malloc(sizeof($*1_type)*PyObject_Length($input));
  if (!convert_darray($input,temp)) {
    SWIG_fail;
  }

  if ($1) { free($1); };

  $1 = &temp[0];

  // arg1 is *SolverOptions. May be version dependent, how to get
  // this?
  if (arg1) arg1->dSize = PyObject_Length($input);

 }

// output lists

%typemap(out) (int *iparam) {
  $1_type piparam = $1;
  PyObject* return_pylist = PyList_New(0);
  for (int i=0; i<arg1->iSize; ++i)
  {
    PyList_Append(return_pylist, PyLong_FromLong(*piparam++));
  }
  $result = return_pylist;
 }


%typemap(out) (double *dparam) {
  $1_type pdparam = $1;
  PyObject* return_pylist = PyList_New(0);
  for (int i=0; i<arg1->dSize; ++i)
  {
    PyList_Append(return_pylist, PyFloat_FromDouble(*pdparam++));
  }
  $result = return_pylist;
 }

 

// LCP
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
  }
}

%extend SolverOptions
{
  SolverOptions()
    {
      SolverOptions *SO;
      SO = (SolverOptions *) malloc(sizeof(SolverOptions));
      return SO;
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
      M->matrix0 = (double *) array_data(array);
      LC->size = M->size0;
      LC->M = M;
      LC->q = (double *) array_data(vector);

      return LC;

    }


  ~LinearComplementarityProblem()
  {
    freeLinearComplementarity_problem(self);
  }

};



%typemap(out) (double* q) {
  npy_intp dims[1];

  dims[0] = arg1->dimension * arg1->numberOfContacts;
  if (arg1->q)
  {
    $result = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, arg1->q);
  }
  else
    SWIG_fail;
 }

// FrictionContact
%include "FrictionContactProblem.h"
%include "FrictionContact3D_Solvers.h"
%include "Friction_cst.h"
%include "FrictionContact3D_AlartCurnier.h"
%include "FrictionContact3D_globalAlartCurnier.h"

%extend FrictionContactProblem
{
  FrictionContactProblem(PyObject *dim, PyObject *o1, PyObject *o2, PyObject *o3)
    {

      int is_new_object1, is_new_object2, is_new_object3; // useless here ?

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
      M->matrix0 = (double *) array_data(array);
      FC->dimension = (int) PyInt_AsLong(dim);
      FC->numberOfContacts = M->size0 / FC->dimension;
      FC->M = M;
      FC->q = (double *) array_data(vector);
      FC->mu = (double *) array_data(mu_vector);

      return FC;
    }

  ~FrictionContactProblem()
  {
    freeFrictionContact_problem(self);
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


