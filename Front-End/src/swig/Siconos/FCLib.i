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

%module FCLib

%include start.i 
%fragment("NumPy_Fragments");
%{
extern "C"
{
#include "fclib.h"
}
%}

 /* fclib_solution */
%typemap(in) (double *v) (PyArrayObject* array=NULL, int is_new_object) {
  array = obj_to_array_allow_conversion($input, NPY_DOUBLE, &is_new_object);
  $1 = (double *) array_data(array);
 }

%apply (double *v) { (double *u) } 
%apply (double *v) { (double *r) }
%apply (double *v) { (double *l) }

%{
  static int convert_fcsol_array(PyObject *input, struct fclib_solution *ptr) {
    void *argp1=0;
    int res1=0;
    int newmem=0;
    int i;
    if (!PySequence_Check(input)) {
      PyErr_SetString(PyExc_TypeError,"Expecting a sequence");
      return 0;
    }
    for (i =0; i <  PyObject_Length(input); i++) {
      PyObject *o = PySequence_GetItem(input,i);
      res1 = SWIG_ConvertPtrAndOwn(o, &argp1, SWIGTYPE_p_fclib_solution, 0 |  0 , &newmem);
      if (!SWIG_IsOK(res1)) {
        Py_XDECREF(o);
        PyErr_SetString(PyExc_ValueError,"Expecting a sequence of fclib_solution");
        return 0;
      }
      ptr[i] = *(reinterpret_cast<fclib_solution*>(argp1));

//      if (ptr[i] == -1 && PyErr_Occurred())
//        return 0;
      Py_DECREF(o);
  }
  return 1;
  }
%}

%typemap(in) (int number_of_guesses,  fclib_solution *guesses) (struct fclib_solution* temp) { 
  
  temp = (fclib_solution *) malloc(sizeof(fclib_solution)*PyObject_Length($input));
  if (!convert_fcsol_array($input,temp)) {
    SWIG_fail;
  }
  $1 = PyObject_Length($input);
  $2 = &temp[0];
 }


%include fclib.h



