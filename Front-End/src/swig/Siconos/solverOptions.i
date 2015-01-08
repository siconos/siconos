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

// Common stuff

#if (SWIG_VERSION >= 0x020005)
#define FE_SWIG_INTERNAL_MEMBER _
#else
#define FE_SWIG_INTERNAL_MEMBER 
#endif

%{
#include "Newton_Methods.h"
%}

#define FPyArray_SimpleNewFromData(nd, dims, typenum, data)             \
  PyArray_New(&PyArray_Type, nd, dims, typenum, NULL,                   \
              data, 0, NPY_ARRAY_FARRAY, NULL)

// std python sequence -> C array

// int
%{
  static int convert_iarray(PyObject *input, int *ptr) {
  int i;
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
  int i;
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
  Py_DECREF($result);
  $result = PyInt_FromLong(*$1);
 }

%typemap(argout) (double *error) {
  $result = PyFloat_FromDouble(*$1);
 }


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

  npy_intp this_iparam_dim[1];
  this_iparam_dim[0] = arg1->iSize;

  $result = PyArray_SimpleNewFromData(1,this_iparam_dim,NPY_INT,piparam);
 }


%typemap(out) (double *dparam) {
  $1_type pdparam = $1;
  
  npy_intp this_dparam_dim[1];
  this_dparam_dim[0] = arg1->dSize;

  $result = PyArray_SimpleNewFromData(1, this_dparam_dim,NPY_DOUBLE,pdparam);
 }

%fragment("NumPy_Fragments");

%{

  void collectStatsIterationCallback(void *env,
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
  pycallback = new Callback(); // free in deleteSolverOptions
  pycallback->env = $input;
  pycallback->collectStatsIteration = &collectStatsIterationCallback;

  $1 = pycallback;

 }


%{
#include "SolverOptions.h"
%}

%include "SolverOptions.h"

