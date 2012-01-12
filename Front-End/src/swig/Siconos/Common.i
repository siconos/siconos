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

// Common stuff


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

  for (i =0; i < PyObject_Length(input); i++) {
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
%typemap(in, numinputs=0) int *info (int temp) {
  $1 = &temp;
}

%typemap(in, numinputs=0) double *error (double temp) {
  $1 = &temp;
}

%typemap(argout) (int *info) {
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

  $result = PyArray_SimpleNewFromData(1,this_dparam_dim,NPY_DOUBLE,pdparam);
 }


//PyArray_UpdateFlags does not seem to have any effect
//>>> r = K.FirstOrderLinearTIR()
//>>> r.setCPtr([[1,2,3],[4,5,6]])
//>>> C=r.C()
//>>> C.flags
//  C_CONTIGUOUS : True
//  F_CONTIGUOUS : False            <---- !!!
//  OWNDATA : False
//  WRITEABLE : True
//  ALIGNED : True
//  UPDATEIFCOPY : False
//
//
// with this macro : ok
#define FPyArray_SimpleNewFromData(nd, dims, typenum, data)             \
  PyArray_New(&PyArray_Type, nd, dims, typenum, NULL,                   \
              data, 0, NPY_FARRAY, NULL)


// copy shared ptr reference in a base PyCObject 
#define PYARRAY_FROM_SHARED_DATA(NDIM,DIMS,NAME,RESULT)                 \
  PyObject* pyarray = FPyArray_SimpleNewFromData(NDIM,                  \
                                                 DIMS,                  \
                                                 NPY_DOUBLE,            \
                                                 NAME->getArray());     \
  SharedPointerKeeper* savedSharedPointer = new                         \
    SharedPointerKeeper(boost::static_pointer_cast<void>(NAME));        \
  reinterpret_cast<PyArrayObject*>(pyarray)->base =                     \
    PyCObject_FromVoidPtr((void*) savedSharedPointer,                   \
                          &sharedPointerKeeperDelete);                  \
  RESULT = pyarray
