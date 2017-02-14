// -*- C++ -*-
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.
//
// Copyright 2016 INRIA.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//

//	

// Common stuff

#if (SWIG_VERSION >= 0x020005)
#define FE_SWIG_INTERNAL_MEMBER _
#else
#define FE_SWIG_INTERNAL_MEMBER 
#endif

%{
#include "Newton_methods.h"
%}

#ifdef SWIGPYTHON
#define FPyArray_SimpleNewFromData(nd, dims, typenum, data)             \
  PyArray_New(&PyArray_Type, nd, dims, typenum, NULL,                   \
              data, 0, NPY_ARRAY_FARRAY, NULL)

// int
%{
  static int convert_iarray(PyObject *input, int *ptr) {
  int i = 0;
  if (!PySequence_Check(input)) {
      PyErr_SetString(PyExc_TypeError,"Expecting a sequence");
      return 0;
  }
  for (i =0; i <  PyObject_Length(input); i++) {
      PyObject *o = PySequence_GetItem(input,i);
      if (!PyInt_Check(o)) {
        Py_XDECREF(o);
        PyErr_SetString(PyExc_ValueError,"Expecting a sequence of ints");
        return 0;
      }
      ptr[i] = (int) PyInt_AsLong(o);

      if (ptr[i] == -1 && PyErr_Occurred())
        return 0;
      Py_DECREF(o);
  }
  return 1;
}
%}

// double
%{
static int convert_darray(PyObject *input, double *ptr) {
  int i = 0;
  if (!PySequence_Check(input)) {
      PyErr_SetString(PyExc_TypeError,"Expecting a sequence");
      return 0;
  }

  for (i =0; i < PyObject_Length(input); ++i) {
      PyObject *o = PySequence_GetItem(input,i);
      if (!PyFloat_Check(o)) {
         Py_XDECREF(o);
         PyErr_SetString(PyExc_ValueError,"Expecting a sequence of floats");
         return 0;
      }
      ptr[i] = PyFloat_AsDouble(o);

      if (ptr[i] == -1 && PyErr_Occurred())
        return 0;

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
#endif /* SWIGPYTHON */

// cleanup
%typemap(freearg) (int sized, double *dparam) {
  free($2);
 }

%typemap(freearg) (int sizei, int *iparam) {
  free($2);
 }

// info, error results
%typemap(in, numinputs=0) double *error (double temp) {
  $1 = &temp;
}

%typemap(argout) (int *info) {
#ifdef SWIGPYTHON
  target_mem_mgmt_instr($result);
#endif
  $result = SWIG_From_int(*$1);
 }

%typemap(argout) (double *error) {
#ifdef SWIGPYTHON
  target_mem_mgmt_instr($result);
#endif
  $result = SWIG_From_double(*$1);
 }

#ifdef SWIGPYTHON
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
#endif /* SWIGPYTHON */

#ifdef SWIGMATLAB
%typemap(in) (int* iparam) (SN_ARRAY_TYPE* array = NULL, int is_new_object = 0){
   array = obj_to_sn_vector_int($input, &is_new_object);

  if (!array)
  {
   SWIG_exception_fail(SWIG_TypeError, "Could not get a int from the object");
  }

  if (CHECK_ARRAY_VECTOR(array))
  {
   SWIG_exception_fail(SWIG_TypeError, "The given object does not have the right structure. We expect a vector (or list, tuple, ...)");
  }

  $1 = (int *) array_data(array);
}
 
 %typemap(memberin) (int* iparam) {
  if (arg1) arg1->iSize = array_size(array2, 1) == 1 ? array_size(array2, 0) : array_size(array2, 1);
  if (!$1) { $1 = (int*)malloc(arg1->iSize * sizeof(int)); }
  else { $1 = (int*)realloc($1, arg1->iSize * sizeof(int)); }

  memcpy($1, $input, arg1->iSize * sizeof(int));
}

%apply (double *z) { (double *dparam) };

%typemap(memberin) (double *dparam) {
  
  if (arg1) arg1->dSize = array_size(array2, 1) == 1 ? array_size(array2, 0) : array_size(array2, 1);
  if (!$1) { $1 = (double*)malloc(arg1->dSize * sizeof(double)); }
  else { $1 = (double*)realloc($1, arg1->dSize * sizeof(double)); }

  memcpy($1, $input, arg1->dSize * sizeof(double));
}
#endif /* SWIGMATLAB */

// output lists

%typemap(out) (int *iparam) {
  C_to_target_lang1_int($result, arg1->iSize, $1, SWIG_fail);
 }


%typemap(out) (double *dparam) {
  C_to_target_lang1($result, arg1->dSize, $1, SWIG_fail);
 }

#ifdef SWIGPYTHON
%fragment("NumPy_Fragments");

%{

  static void collectStatsIterationCallback(void *env,
                            int size, double* z,
                            double* Fz,
                            double error,
                            void* extra_data)
  {
    // A little bit of black magic
    PyObject* py_tuple;
    if (extra_data)
    {
      switch (*(int*)extra_data)
      {
        case NEWTON_STATS_ITERATION:
        {
          newton_stats* stats = (newton_stats*) extra_data;
          py_tuple = PyTuple_New(3);
          PyObject* py_merit_value = PyFloat_FromDouble(stats->merit_value);
          PyTuple_SetItem(py_tuple, 0, py_merit_value);
          PyObject* py_alpha =  PyFloat_FromDouble(stats->alpha);
          PyTuple_SetItem(py_tuple, 1, py_alpha);
          PyObject* py_status = PyTuple_New(2);
          PyObject* py_newton_step;
          if (stats->status & NEWTON_STATS_NEWTON_STEP)
          {
            Py_INCREF(Py_True);
            py_newton_step = Py_True;
          }
          else
          {
            Py_INCREF(Py_False);
            py_newton_step = Py_False;
          }
          PyTuple_SetItem(py_status, 0, py_newton_step);

          PyObject* py_desc_dir;
          if (stats->status & NEWTON_STATS_DESC_DIR)
          {
            Py_INCREF(Py_True);
            py_desc_dir = Py_True;
          }
          else
          {
            Py_INCREF(Py_False);
            py_desc_dir = Py_False;
          }
          PyTuple_SetItem(py_status, 1, py_desc_dir);

          PyTuple_SetItem(py_tuple, 2, py_status);
          break;
        }

        default:
          py_tuple = PyTuple_New(0);
      }
    }
    else
    {
      py_tuple = PyTuple_New(0);
    }
    if (PyCallable_Check((PyObject*) env))
    {

      npy_intp dim[1];
      dim[0] = size;

      PyObject* py_z = FPyArray_SimpleNewFromData(1, dim,
                                                  NPY_DOUBLE,
                                                  z);

      PyObject* py_Fz = FPyArray_SimpleNewFromData(1, dim,
                                                   NPY_DOUBLE,
                                                   Fz);

      PyObject* py_error = PyFloat_FromDouble(error);

      PyObject* py_args = PyTuple_New(4);
      PyTuple_SetItem(py_args, 0, py_z);
      PyTuple_SetItem(py_args, 1, py_Fz);
      PyTuple_SetItem(py_args, 2, py_error);
      PyTuple_SetItem(py_args, 3, py_tuple);

      PyGILState_STATE gstate;
      gstate = PyGILState_Ensure();

      PyObject* py_out = PyObject_CallObject((PyObject*) env, py_args);

      Py_DECREF(py_args);
      Py_XDECREF(py_out);

      PyGILState_Release(gstate);
    }
    else
    {
      PyErr_SetString(PyExc_TypeError,"Expecting a callable callback");
    }
  };

%}


%typemap(in) (Callback* callback) (Callback* pycallback) { 
  // %typemap(in) (Callback* callback) (Callback* pycallback)
  if ((arg1)->callback)
  {
    free(arg1->callback);
  }
  pycallback = (Callback*) malloc(sizeof(Callback)); // free in solver_options_delete
  pycallback->env = $input;
  pycallback->collectStatsIteration = &collectStatsIterationCallback;

  $1 = pycallback;

 }
#endif /* SWIGPYTHON */

%{
#include "SolverOptions.h"
%}

%include "SolverOptions.h"
