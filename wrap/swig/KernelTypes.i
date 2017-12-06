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
// SWIG interface for Siconos Kernel types

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

// SimpleMatrix operators
%rename  (__add__) operator+;
%rename  (__less__) operator-;
%rename  (__mul__) operator*;
%rename  (__div__) operator/;
%rename  (__iadd__) operator+=;
%rename  (__iless__) operator-=;
%rename  (__imul__) operator*=;
%rename  (__idiv__) operator/=;
%rename  (__eq__) operator==;
%rename  (__ne__) operator!=;

// Default swig typemap checks only for shared_ptr, but in some cases
// it is a normal pointer.  We insert a normal pointer check before
// the usual shared_ptr check.
%define REF_PTR(TYPE)
%typemap(in) TYPE&
%{
  {
    int newmem = 0;
    void *argp$argnum;
    int res$argnum = SWIG_ConvertPtrAndOwn($input, &argp$argnum, $descriptor(TYPE *),  0 , &newmem);
    if (SWIG_IsOK(res$argnum)) {
      $1 = static_cast<TYPE*>(argp$argnum);
    }
    else
    {
      res$argnum = SWIG_ConvertPtrAndOwn($input, &argp$argnum, $descriptor(std11::shared_ptr< TYPE >*),  0 , &newmem);
      if (!SWIG_IsOK(res$argnum)) {
        %argument_fail(SWIG_ValueError, "$type", $symname, $argnum);
      }
      if (!argp$argnum) {
        %argument_fail(SWIG_ValueError, "$type", $symname, $argnum);
      }
      if (newmem & SWIG_CAST_NEW_MEMORY) {
        std11::shared_ptr< TYPE > tempshared$argnum
          (*reinterpret_cast< std11::shared_ptr<  TYPE > * >(argp$argnum));
        delete reinterpret_cast< std11::shared_ptr<  TYPE > * >(argp$argnum);
        $1 = const_cast< TYPE * >(tempshared$argnum.get());
      } else {
        $1 = const_cast< TYPE * >(reinterpret_cast< std11::shared_ptr<  TYPE > * >(argp$argnum)->get());
      }
    }
  }
%}
%typemap(in) (const TYPE&)
%{
  {
    int newmem = 0;
    void *argp$argnum;
    int res$argnum = SWIG_ConvertPtrAndOwn($input, &argp$argnum, $descriptor(TYPE *),  0 , &newmem);
    if (SWIG_IsOK(res$argnum)) {
      $1 = static_cast<TYPE*>(argp$argnum);
    }
    else
    {
      res$argnum = SWIG_ConvertPtrAndOwn($input, &argp$argnum, $descriptor(std11::shared_ptr< TYPE >*),  0 , &newmem);
      if (!SWIG_IsOK(res$argnum)) {
        %argument_fail(SWIG_ValueError, "$type", $symname, $argnum);
      }
      if (!argp$argnum) {
        %argument_fail(SWIG_ValueError, "$type", $symname, $argnum);
      }
      if (newmem & SWIG_CAST_NEW_MEMORY) {
        std11::shared_ptr< TYPE > tempshared$argnum
          (*reinterpret_cast< std11::shared_ptr< TYPE > * >(argp$argnum));
        delete reinterpret_cast< std11::shared_ptr< TYPE > * >(argp$argnum);
        $1 = const_cast< TYPE * >(tempshared$argnum.get());
      } else {
        $1 = const_cast< TYPE * >(reinterpret_cast< std11::shared_ptr< TYPE > * >(argp$argnum)->get());
      }
    }
  }
%}
%typecheck(SWIG_TYPECHECK_POINTER) (TYPE&)
%{
  {
    void *ptr;
    int res = SWIG_ConvertPtr($input, &ptr, $descriptor(TYPE *), 0);
    $1 = SWIG_CheckState(res);
    if (!$1) {
      res = SWIG_ConvertPtr($input, &ptr, $descriptor(std11::shared_ptr< TYPE > *), 0);
      $1 = SWIG_CheckState(res);
    }
  }
%}
%typecheck(SWIG_TYPECHECK_POINTER) (const TYPE&)
%{
  {
    void *ptr;
    int res = SWIG_ConvertPtr($input, &ptr, $descriptor(TYPE *), 0);
    $1 = SWIG_CheckState(res);
    if (!$1) {
      res = SWIG_ConvertPtr($input, &ptr, $descriptor(std11::shared_ptr< TYPE > *), 0);
      $1 = SWIG_CheckState(res);
    }
  }
%}
%enddef

#undef PY_REGISTER_WITHOUT_DIRECTOR

%define PY_REGISTER_WITHOUT_DIRECTOR(TYPE)
%inline
%{
#include "TYPE.hpp"
%}
%rename  (__getitem__) TYPE ## ::operator[];
%rename  (__add__) TYPE ## ::operator+;
%rename  (__mul__) TYPE ## ::operator*;
%rename  (__div__) TYPE ## ::operator/;
%rename  (__iadd__) TYPE ## ::operator+=;
%rename  (__imul__) TYPE ## ::operator*=;
%rename  (__idiv__) TYPE ## ::operator/=;
%rename  (__eq__) TYPE ## ::operator==;
%rename  (__ne__) TYPE ## ::operator!=;
%rename  (__copy__) TYPE ## ::operator=;
%ignore STD11::enable_shared_from_this<TYPE>;
%shared_ptr(STD11::enable_shared_from_this<TYPE>); // warning 520 suppression
%template (shared ## TYPE) STD11::enable_shared_from_this<TYPE>;
%shared_ptr(TYPE);
%make_picklable(TYPE, Kernel);
REF_PTR(TYPE);
%enddef


//%typemap(directorin) TYPE& ()
//%{
//  // %typemap(directorin) (TYPE&) ()
//  // swig issue shared pointer check in wrappers even if arg is a ref
//  {
//    SP::TYPE myptemp(createSPtr##TYPE($1));
//    $input = SWIG_NewPointerObj(SWIG_as_voidptr(&myptemp),
//                                $descriptor(SP::TYPE *), 0);
//  }
//%}
//
//%typemap(directorout) TYPE& ()
//%{
//  // %typemap(directorout) (TYPE&) ()
//  // swig issue shared pointer check in wrappers even if arg is a ref
//  {
//    SP::TYPE myptemp(createSPtr##TYPE($1));
//    $result = SWIG_NewPointerObj(SWIG_as_voidptr(&myptemp),
//                                 $descriptor(SP::TYPE *), 0);
//  }
//%}
//


PY_REGISTER_WITHOUT_DIRECTOR(SiconosMatrix);
PY_REGISTER_WITHOUT_DIRECTOR(SimpleMatrix);
PY_REGISTER_WITHOUT_DIRECTOR(SiconosVector);
PY_REGISTER_WITHOUT_DIRECTOR(BlockVector);

// set the base of the pyarray to a PyCapsule or PyCObject created from the shared_ptr
%{
static inline void fillBasePyarray(PyObject* pyarray, SharedPointerKeeper* savedSharedPointer)
{
  PyObject* cap =
#ifdef SWIGPY_USE_CAPSULE
    PyCapsule_New((void*)( savedSharedPointer), SWIGPY_CAPSULE_NAME, sharedPointerKeeperDeleteCap);
#else
    PyCObject_FromVoidPtr((void*)(savedSharedPointer), sharedPointerKeeperDelete);
#endif

#if NPY_API_VERSION < 0x00000007
  PyArray_BASE((PyArrayObject*)pyarray) = cap;
#else
  PyArray_SetBaseObject((PyArrayObject*) pyarray,cap);
#endif
}
%}

// copy shared ptr reference in a base PyCObject || PyCapsule
#define PYARRAY_FROM_SHARED_SICONOS_DATA(TYPE,NDIM,DIMS,NAME,RESULT)\
  PyObject* pyarray = FPyArray_SimpleNewFromData(NDIM,              \
                                                 DIMS,              \
                                                 TYPE,              \
                                                 NAME->getArray()); \
  SharedPointerKeeper* savedSharedPointer = new                     \
    SharedPointerKeeper(std11::static_pointer_cast<void>(NAME));    \
  fillBasePyarray(pyarray, savedSharedPointer);                         \
  RESULT = pyarray

#define PYARRAY_FROM_SHARED_STL_VECTOR(TYPE,NDIM,DIMS,NAME,RESULT)      \
  PyObject* pyarray = FPyArray_SimpleNewFromData(NDIM,                  \
                                                 DIMS,                  \
                                                 TYPE,                  \
                                                 &(*NAME)[0]);          \
  SharedPointerKeeper* savedSharedPointer = new                         \
    SharedPointerKeeper(std11::static_pointer_cast<void>(NAME));        \
  fillBasePyarray(pyarray, savedSharedPointer);                             \
  RESULT = pyarray

#define PYARRAY_FROM_SHARED_SICONOS_DATA_REF(TYPE,NDIM,DIMS,NAME,RESULT)\
  PyObject* pyarray = FPyArray_SimpleNewFromData(NDIM,              \
                                                 DIMS,              \
                                                 TYPE,              \
                                                 NAME.getArray()); \
  RESULT = pyarray



%typemap(in) PyArrayObject* {
   $1 = (PyArrayObject*) $input;
}

//////////////////////////////////////////////////////////////////////////////
// allow unsigned int be numpy types
%typecheck(SWIG_TYPECHECK_INTEGER,fragment="SWIG_AsVal_int") unsigned int {
  int ecode = SWIG_AsVal_int($input, NULL);
  $1 = SWIG_IsOK(ecode);
  if (!$1) {
    if (PyArray_CheckAnyScalar($input)) {
      int x = PyArray_PyIntAsInt($input);
      $1 = !(x == -1 && PyErr_Occurred());
      PyErr_Clear();
    }
  }
}
%typemap(in,fragment="SWIG_AsVal_int") unsigned int {
  int x;
  int ecode = SWIG_AsVal_int($input, &x);
  if (SWIG_IsOK(ecode)) {
    if (x < 0)
      %argument_fail(ecode, "$type", $symname, $argnum);
  } else if (PyArray_CheckScalar($input)) {
    // TODO: By treating as int, we lose half the range of the
    // unsigned int!  Okay for current uses of unsigned int in
    // Siconos.
    x = PyArray_PyIntAsInt($input);
    if (x == -1 && PyErr_Occurred())
      %argument_fail(SWIG_TypeError, "$type", $symname, $argnum);
    if (x < 0) {
      %argument_fail(SWIG_ValueError, "$type", $symname, $argnum);
    }
  }
  $1 = (unsigned int) x;
}

//////////////////////////////////////////////////////////////////////////////
// allow integers to be numpy types
%typecheck(SWIG_TYPECHECK_INTEGER,fragment="SWIG_AsVal_int") int {
  int ecode = SWIG_AsVal_int($input, NULL);
  $1 = SWIG_IsOK(ecode);
  if (!$1) {
    if (PyArray_CheckAnyScalar($input)) {
      int x = PyArray_PyIntAsInt($input);
      $1 = !(x == -1 && PyErr_Occurred());
      PyErr_Clear();
    }
  }
}
%typemap(in,fragment="SWIG_AsVal_int") int {
  int ecode = SWIG_AsVal_int($input, &$1);
  if (!SWIG_IsOK(ecode) && PyArray_CheckAnyScalar($input))
  {
    $1 = PyArray_PyIntAsInt($input);
    if ($1 == -1 && PyErr_Occurred())
      %argument_fail(SWIG_TypeError, "$type", $symname, $argnum);
  }
}


//////////////////////////////////////////////////////////////////////////////
// allow double to be numpy types
%typecheck(SWIG_TYPECHECK_DOUBLE,fragment="SWIG_AsVal_double") double {
  $1 = PyArray_CheckAnyScalar($input);
}
%typemap(in,fragment="SWIG_AsVal_double") double {
  int ecode = SWIG_AsVal_double($input, &$1);
  if (!SWIG_IsOK(ecode) && PyArray_CheckAnyScalar($input)) {
    $1 = PyFloat_AsDouble($input);
    if ($1 == -1.0 && PyErr_Occurred()) {
      PyObject *o = PyNumber_Float($input);
      if (o) {
        $1 = PyFloat_AsDouble(o);
        if ($1 == -1.0 && PyErr_Occurred()) {
          Py_DECREF(o);
          %argument_fail(SWIG_TypeError, "$type", $symname, $argnum);
        }
        Py_DECREF(o);
      }
      else
        %argument_fail(SWIG_TypeError, "$type", $symname, $argnum);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////
// check on input : a numpy array or a TYPE
%define TYPECHECK(TYPE)
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY)
(TYPE)
{
  // %typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY) (TYPE)
  TYPE * ptr;

  int res = SWIG_ConvertPtr($input, (void **) (&ptr), $descriptor(TYPE *), 0);
  int state = SWIG_CheckState(res);
  if (SWIG_IsNewObj(res)) { delete ptr; };
  $1 = is_array($input) || PySequence_Check($input) || state;
}
%enddef
%define TYPECHECK_REF(TYPE)
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY)
(TYPE&)
{
  // %typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY) (TYPE&)
  TYPE * ptr;

  int res = SWIG_ConvertPtr($input, (void **) (&ptr), $descriptor(TYPE *), 0);
  int state = SWIG_CheckState(res);
  if (SWIG_IsNewObj(res)) { delete ptr; };
  if (!state)
  {
    int swig_res = SWIG_ConvertPtr($input, (void **) (&ptr), $descriptor(SP::TYPE *), 0);
    state = SWIG_CheckState(swig_res);
    if (SWIG_IsNewObj(swig_res)) { delete ptr; };
  }
  $1 = is_array($input) || PySequence_Check($input) || state;
}
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY)
(const TYPE&)
{
  // %typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY) (const TYPE&)
  TYPE * ptr;

  int res = SWIG_ConvertPtr($input, (void **) (&ptr), $descriptor(const TYPE *), 0);
  int state = SWIG_CheckState(res);
  if (SWIG_IsNewObj(res)) { delete ptr; };
  if (!state)
  {
    int swig_res = SWIG_ConvertPtr($input, (void **) (&ptr), $descriptor(TYPE *), 0);
    state = SWIG_CheckState(swig_res);
    if (SWIG_IsNewObj(swig_res)) { delete ptr; };
  }
  if (!state)
  {
    int swig_res = SWIG_ConvertPtr($input, (void **) (&ptr), $descriptor(SP::TYPE *), 0);
    state = SWIG_CheckState(swig_res);
    if (SWIG_IsNewObj(swig_res)) { delete ptr; };
  }
  $1 = is_array($input) || PySequence_Check($input) || state;
}
%enddef
//////////////////////////////////////////////////////////////////////////////

%{
#include <Question.hpp>
#include <SiconosVector.hpp>
struct IsDense : public Question<bool>
{
  using SiconosVisitor::visit;

  void visit(const SiconosVector& v)
  {
    answer = v._dense;
  }

  void visit(const BlockVector& v)
  {
    answer = false;
  }
};
%}

%fragment("SiconosVector", "header", fragment="NumPy_Fragments")
{
  PyObject * SiconosVector_to_numpy(SiconosVector const & v)
  {
    npy_intp this_vector_dim[1] = { v.size() };

    PyObject * larray;
    PYARRAY_FROM_SHARED_SICONOS_DATA_REF(NPY_DOUBLE, 1, this_vector_dim, v, larray);

    return larray;
  }

  PyObject * SP_SiconosVector_to_numpy(SP::SiconosVector v)
  {
    npy_intp this_vector_dim[1] = { v->size() };
    PyObject* lresult;
    PYARRAY_FROM_SHARED_SICONOS_DATA(NPY_DOUBLE, 1, this_vector_dim, v, lresult);
    return lresult;
  }

  SP::SiconosVector SP_SiconosVector_from_numpy(PyObject* vec, PyArrayObject** array_p, int* is_new_object)
  {
    if (vec==Py_None)
      return SP::SiconosVector();

    PyArrayObject* array = obj_to_array_fortran_allow_conversion(vec, NPY_DOUBLE, is_new_object);

    if (!array)
    {
     PyErr_SetString(PyExc_TypeError, "Could not get array obj from the python object");
     PyObject_Print(vec, stderr, 0);
     return std11::shared_ptr<SiconosVector>();
    }
    if (!require_dimensions(array,1) ||
      !require_native(array) || !require_fortran(array))
    {
      PyErr_SetString(PyExc_TypeError, "The given object does not have the right structure. We expect a 1 dimensional array (or list, tuple, ...)");
      PyObject_Print(vec, stderr, 0);
      return std11::shared_ptr<SiconosVector>();
    }

    SP::SiconosVector tmp;
    tmp.reset(new SiconosVector(array_size(array,0)));
    // copy : with SiconosVector based on resizable std::vector there is
    // no other way
    memcpy(tmp->getArray(),array_data(array),array_size(array,0)*sizeof(double));

    // for cleanup
    *array_p = array;
    return tmp;
  }

  SP::SiconosVector SP_SiconosVector_in(PyObject* vec, PyArrayObject** array_p, int* is_new_object)
  {
    void *argp1=0;
    int res1=0;
    int newmem = 0;
    std11::shared_ptr<SiconosVector> tempshared1 ;
    std11::shared_ptr<SiconosVector> *smartarg1 = NULL;

    // try a conversion from std11::shared_ptr<SiconosVector>
    res1 = SWIG_ConvertPtrAndOwn(vec, &argp1, $descriptor(std11::shared_ptr<SiconosVector> *), 0 |  0 , &newmem);
    if (SWIG_IsOK(res1) && argp1)
    {
      if (newmem & SWIG_CAST_NEW_MEMORY)
      {
        // taken from generated code
        tempshared1 = *reinterpret_cast< std11::shared_ptr<SiconosVector> * >(argp1);
        delete reinterpret_cast< std11::shared_ptr<SiconosVector> * >(argp1);
        return tempshared1;
      }
      else {
        smartarg1 = reinterpret_cast< std11::shared_ptr<SiconosVector> * >(argp1);
        return  *smartarg1;
      }
    }
    else
    {
      return SP_SiconosVector_from_numpy(vec, array_p, is_new_object);
    }
  }

  SiconosVector* SiconosVector_in(PyObject* vec, PyArrayObject** array_p, int* is_new_object, std::vector<SP::SiconosVector>& keeper)
  {
    void *argp1=0;
    int res0=0;
    int res1=0;
    std11::shared_ptr<SiconosVector> tempshared1 ;
    std11::shared_ptr<SiconosVector> *smartarg1 = NULL;
    SiconosVector* smartarg2 = NULL;

    // try a conversion from std11::shared_ptr<SiconosVector>
    res0 = SWIG_ConvertPtr(vec, &argp1, $descriptor(std11::shared_ptr<SiconosVector> *), 0 |  0);
    if (SWIG_IsOK(res0) && argp1)
    {
      // return the raw pointer -- assuming it will be used
      // temporarily, therefore no ownership or keeper needed.
      smartarg1 = reinterpret_cast< std11::shared_ptr<SiconosVector> * >(argp1);
      return smartarg1->get();
    }

    // try a conversion from SiconosVector
    res1 = SWIG_ConvertPtr(vec, &argp1, $descriptor(SiconosVector *), 0 |  0);
    if (SWIG_IsOK(res1) && argp1)
    {
      // no newmem & SWIG_CAST_NEW_MEMORY + tempshared for non SP
      smartarg2 = reinterpret_cast< SiconosVector * >(argp1);
      return  smartarg2;
    }
    else
    {
      SP::SiconosVector tmp = SP_SiconosVector_from_numpy(vec, array_p, is_new_object);
      keeper.push_back(tmp);
      return tmp.get();
    }
  }

  PyObject * SP_SiconosVector_directorin(SP::SiconosVector v)
  {
    if(v)
    {
      if (ask<IsDense>(*v))
      {
        return SP_SiconosVector_to_numpy(v);
      }
      else
      {
        // not a dense vector : no conversion
        return SWIG_NewPointerObj(SWIG_as_voidptr(&v), $descriptor(SP::SiconosVector *),  0 );
      }
    }
    else
    {
      Py_INCREF(Py_None);
      return Py_None;
    }
  }

  PyObject * SP_SiconosVector_out(SP::SiconosVector result, bool l_upcall, PyObject* obj0)
  {
    // call from director?
    if (l_upcall)
    {
      // result from C++ method, return the pointer
      return SWIG_NewPointerObj(SWIG_as_voidptr(&result), $descriptor(SP::SiconosVector *),  0 );
    }
    // call from python : return numpy from SiconosVector
    else
    {
      if (result)
      {
        // /!\ need check for a dense vector!
        return SP_SiconosVector_to_numpy(result);

      }
      else
      {
        Py_INCREF(Py_None);
        return Py_None;
      }
    }
  }

  SP::SiconosVector SiconosVector_from_python(PyObject* obj)
  {
    void * swig_argp;
    SP::SiconosVector c_result;

    // try a conversion from SP::SiconosVector
    int swig_res = SWIG_ConvertPtr(obj, &swig_argp, $descriptor(SP::SiconosVector *),  0  | 0);

    if (!SWIG_IsOK(swig_res))
    {
      // try a conversion from numpy
      PyArrayObject* array = NULL;
      int is_new_object = 0;
      c_result = SP_SiconosVector_from_numpy(obj, &array, &is_new_object);
      if (array && is_new_object) { Py_DECREF(array); }
    }
    else if (swig_argp)
    {
      c_result = *(reinterpret_cast< SP::SiconosVector * >(swig_argp));
      if (SWIG_IsNewObj(swig_res)) delete reinterpret_cast< SP::SiconosVector * >(swig_argp);
    }
    return c_result;
  }

}

%fragment("SiconosMatrix", "header", fragment="NumPy_Fragments")
{
  PyObject* SiconosMatrix_to_numpy(SiconosMatrix& m)
  {
    if (m.num() == 1)
    {
      npy_intp this_matrix_dim[2];
      this_matrix_dim[0] = m.size(0);
      this_matrix_dim[1] = m.size(1);

      PyObject * linput;
      PYARRAY_FROM_SHARED_SICONOS_DATA_REF(NPY_DOUBLE,2,this_matrix_dim, m, linput);
      return linput;
    }
    else
    {
      // not a dense matrix : no conversion
      return SWIG_NewPointerObj(SWIG_as_voidptr(&m), $descriptor(SiconosMatrix *),  0 );
    }
  }

  PyObject* SiconosMatrix_to_numpy(SimpleMatrix& m)
  {
    if (m.num() == 1)
    {
      npy_intp this_matrix_dim[2];
      this_matrix_dim[0] = m.size(0);
      this_matrix_dim[1] = m.size(1);

      PyObject * linput;
      PYARRAY_FROM_SHARED_SICONOS_DATA_REF(NPY_DOUBLE,2,this_matrix_dim, m, linput);
      return linput;
    }
    else
    {
      // not a dense matrix : no conversion
      return SWIG_NewPointerObj(SWIG_as_voidptr(&m), $descriptor(SimpleMatrix *),  0 );
    }
  }

  PyObject* SiconosMatrix_to_numpy(SP::SimpleMatrix m)
  {
    if (m && m->size(0) > 0 && m->size(1) > 0)
    {
      if (m->num() == 1)
      {
        npy_intp this_matrix_dim[2];
        this_matrix_dim[0] = m->size(0);
        this_matrix_dim[1] = m->size(1);

        PyObject * linput;
        PYARRAY_FROM_SHARED_SICONOS_DATA(NPY_DOUBLE,2,this_matrix_dim, m, linput);
        return linput;
      }
      else
      {
        // not a dense matrix : no conversion
        return SWIG_NewPointerObj(SWIG_as_voidptr(&m), $descriptor(SP::SimpleMatrix *),  0 );
      }
    }
    else
    {
      Py_INCREF(Py_None);
      return Py_None;
    }
  }

  PyObject* SiconosMatrix_to_numpy(SP::SiconosMatrix m)
  {
    if (m && m->size(0) > 0 && m->size(1) > 0)
    {
      if (m->num() == 1)
      {
        npy_intp this_matrix_dim[2];
        this_matrix_dim[0] = m->size(0);
        this_matrix_dim[1] = m->size(1);

        PyObject * linput;
        PYARRAY_FROM_SHARED_SICONOS_DATA(NPY_DOUBLE,2,this_matrix_dim, m, linput);
        return linput;
      }
      else
      {
        // not a dense matrix : no conversion
        return SWIG_NewPointerObj(SWIG_as_voidptr(&m), $descriptor(SP::SiconosMatrix *),  0 );
      }
    }
    else
    {
      Py_INCREF(Py_None);
      return Py_None;
    }
  }

  PyObject* SiconosMatrix_to_numpy(SP::SimpleMatrix m, bool l_upcall)
  {
    if (l_upcall)
    {
      // result from C++ method, return the pointer
      return SWIG_NewPointerObj(SWIG_as_voidptr(&m), $descriptor(SP::SimpleMatrix *),  0 );
    }
    // call from python : return numpy from SiconosMatrix
    else
    {
      return SiconosMatrix_to_numpy(m);
    }
  }

  PyObject* SiconosMatrix_to_numpy(SP::SiconosMatrix m, bool l_upcall)
  {
    if (l_upcall)
    {
      // result from C++ method, return the pointer
      return SWIG_NewPointerObj(SWIG_as_voidptr(&m), $descriptor(SP::SiconosMatrix *),  0 );
    }
    // call from python : return numpy from SiconosMatrix
    else
    {
      return SiconosMatrix_to_numpy(m);
    }
  }

  SP::SimpleMatrix SimpleMatrix_from_numpy(PyObject* obj, PyArrayObject** array_p, int* is_new_object)
  {
    if (obj==Py_None)
      return SP::SimpleMatrix();

    PyArrayObject* array = obj_to_array_fortran_allow_conversion(obj, NPY_DOUBLE, is_new_object);
    if (!array)
    {
      PyErr_SetString(PyExc_TypeError, "Could not get array obj from the python object");
      PyObject_Print(obj, stderr, 0);
      return std11::shared_ptr<SimpleMatrix>();
    }

    if (!require_dimensions(array,2) ||
        !require_native(array) || !require_fortran(array))
    {
      PyErr_SetString(PyExc_TypeError, "The given object does not have the right structure. We expect a 2 dimensional array (or list, tuple, ...)");
      PyObject_Print(obj, stderr, 0);
      return std11::shared_ptr<SimpleMatrix>();
    }

    SP::SimpleMatrix result = SP::SimpleMatrix(new SimpleMatrix(array_size(array,0), array_size(array,1)));
    // copy this is due to SimpleMatrix based on resizable std::vector
    memcpy(result->getArray(), array_data(array), array_size(array,0)*array_size(array,1)*sizeof(double));
    // for cleanup
    *array_p = array;
    return result;
  }

  bool SiconosMatrix_from_python(PyObject* obj, PyArrayObject** array_p, int* is_new_object, SP::SimpleMatrix* c_result)
  {
    void * swig_argp;

    int swig_res = SWIG_ConvertPtr(obj, &swig_argp, $descriptor(SP::SimpleMatrix *),  0  | 0);

    if (!SWIG_IsOK(swig_res))
    {
      // try a conversion from numpy
      *c_result = SimpleMatrix_from_numpy(obj, array_p, is_new_object);
      if (!c_result) { return false; }
    }
    else if (swig_argp)
    {
      *c_result = *(reinterpret_cast< SP::SimpleMatrix * >(swig_argp));
      if (SWIG_IsNewObj(swig_res)) delete reinterpret_cast< SP::SimpleMatrix * >(swig_argp);
    }
    return true;
  }

  bool SiconosMatrix_from_python(PyObject* obj, PyArrayObject** array_p, int* is_new_object, SP::SiconosMatrix* c_result)
  {
    void * swig_argp;
    int swig_res = SWIG_ConvertPtr(obj, &swig_argp, $descriptor(SP::SiconosMatrix *),  0  | 0);

    if (!SWIG_IsOK(swig_res))
    {
      *c_result = SimpleMatrix_from_numpy(obj, array_p, is_new_object);
      if (!c_result) { return false; }
    }
    else if (swig_argp)
    {
      *c_result = *(reinterpret_cast< SP::SiconosMatrix * >(swig_argp));
      if (SWIG_IsNewObj(swig_res)) delete reinterpret_cast< SP::SiconosMatrix * >(swig_argp);
    }
    return true;
  }

  bool SiconosMatrix_from_python(PyObject* obj, PyArrayObject** array_p, int* is_new_object, SiconosMatrix** c_result, std::vector<SP::SiconosMatrix>& keeper)
  {
    void * swig_argp;
    int swig_res = SWIG_ConvertPtr(obj, &swig_argp, $descriptor(SiconosMatrix *),  0  | 0);

    if (!SWIG_IsOK(swig_res))
    {
      swig_res = SWIG_ConvertPtr(obj, &swig_argp, $descriptor(SP::SiconosMatrix *),  0  | 0);

      if (!SWIG_IsOK(swig_res))
      {
        SP::SiconosMatrix tmp = SimpleMatrix_from_numpy(obj, array_p, is_new_object);
        if (!tmp) { return false; }
        keeper.push_back(tmp);
        *c_result = tmp.get();
      }
      else if (swig_argp)
      {
        *c_result = (reinterpret_cast< SP::SiconosMatrix * >(swig_argp))->get();
        if (SWIG_IsNewObj(swig_res)) delete reinterpret_cast< SP::SiconosMatrix * >(swig_argp);
      }
    }
    else if (swig_argp)
    {
      *c_result = (reinterpret_cast< SiconosMatrix * >(swig_argp));
      if (SWIG_IsNewObj(swig_res)) delete reinterpret_cast< SiconosMatrix * >(swig_argp);
    }
    return true;
  }

 bool SiconosMatrix_from_python(PyObject* obj, PyArrayObject** array_p, int* is_new_object, SimpleMatrix** c_result, std::vector<SP::SimpleMatrix>& keeper)
  {
    void * swig_argp;
    int swig_res = SWIG_ConvertPtr(obj, &swig_argp, $descriptor(SimpleMatrix *),  0  | 0);

    if (!SWIG_IsOK(swig_res))
    {
      swig_res = SWIG_ConvertPtr(obj, &swig_argp, $descriptor(SP::SimpleMatrix *),  0  | 0);

      if (!SWIG_IsOK(swig_res))
      {
        SP::SimpleMatrix tmp = SimpleMatrix_from_numpy(obj, array_p, is_new_object);
        if (!tmp) { return false; }
        keeper.push_back(tmp);
        *c_result = tmp.get();
      }
      else if (swig_argp)
      {
        *c_result = (reinterpret_cast< SP::SimpleMatrix * >(swig_argp))->get();
        if (SWIG_IsNewObj(swig_res)) delete reinterpret_cast< SP::SimpleMatrix * >(swig_argp);
      }
    }
    else if (swig_argp)
    {
      *c_result = (reinterpret_cast< SimpleMatrix * >(swig_argp));
      if (SWIG_IsNewObj(swig_res)) delete reinterpret_cast< SimpleMatrix * >(swig_argp);
    }
    return true;
  }
}

//////////////////////////////////////////////////////////////////////////////
%typemap(in,fragment="SiconosVector") (std11::shared_ptr<SiconosVector>) (PyArrayObject* array = NULL, int is_new_object = 0)
{
  // %typemap(in,fragment="SiconosVector") (std11::shared_ptr<SiconosVector>)
  $1 = SP_SiconosVector_in($input, &array, &is_new_object);
}

%typemap(in,fragment="SiconosVector")
  const SiconosVector &
  (PyArrayObject* array=NULL, int is_new_object = 0, std::vector<SP::SiconosVector> keeper)
{
  // %typemap(in,fragment="NumPy_Fragments")
  // %TYPE (PyArrayObject* array=NULL, int
  // %is_new_object)
  $1 = SiconosVector_in($input, &array, &is_new_object, keeper);
  if (!$1)
  {
    void * swig_argp;
    int swig_res = SWIG_ConvertPtr($input, &swig_argp, $descriptor(SP::SiconosVector *), 0);
    if (SWIG_IsOK(swig_res) && swig_argp)
    {
      $1 = (%reinterpret_cast(swig_argp, SP::SiconosVector *))->get();
      if (!$1) { SWIG_exception_fail(SWIG_ValueError, "expected vector"); }
    }
    else { SWIG_exception_fail(SWIG_ValueError, "expected vector"); }
  }
}

%typemap(in,fragment="SiconosVector") SiconosVector &
(PyArrayObject* array=NULL, int is_new_object = 0, std::vector<SP::SiconosVector> keeper)
{
  // %typemap(in,fragment="NumPy_Fragments")
  // %TYPE (PyArrayObject* array=NULL, int
  // %is_new_object)
  $1 = SiconosVector_in($input, &array, &is_new_object, keeper);
  if (!$1)
  {
    void * swig_argp;
    int swig_res = SWIG_ConvertPtr($input, &swig_argp, $descriptor(SP::SiconosVector *), 0);
    if (SWIG_IsOK(swig_res) && swig_argp)
    {
      $1 = (%reinterpret_cast(swig_argp, SP::SiconosVector *))->get();
      if (!$1) { SWIG_exception_fail(SWIG_ValueError, "expected vector"); }
    }
    else { SWIG_exception_fail(SWIG_ValueError, "expected vector"); }
  }
}

%typemap(freearg) (std11::shared_ptr<SiconosVector>)
{
  if (is_new_object$argnum && array$argnum)
    { Py_DECREF(array$argnum); }
}

%typemap(freearg) (SiconosVector &)
{
  if (is_new_object$argnum && array$argnum)
    { Py_DECREF(array$argnum); }
}

// director input : TYPE -> numpy
%typemap(directorin, fragment="SiconosVector") (std11::shared_ptr<SiconosVector>) ()
{
  // %typemap(directorin, fragment="SiconosVector") (std11::shared_ptr<SiconosVector>) ()
  $input = SP_SiconosVector_directorin($1_name);
}

%typemap(out, fragment="SiconosVector") (std11::shared_ptr<SiconosVector>)
{
  // %typemap(out) TYPE
  // compile time test to reproduce swig director test. swig does not
  // seem to provides facilities to customize this
  // arg1 is the class instantiation (swig 2.0) => no way to get this as a swig
  // variable.
  typedef BOOST_TYPEOF(arg1) self_type;

  typedef boost::mpl::eval_if<boost::is_polymorphic<self_type >,
    DirectorCast<self_type >,
    DirectorNoCast<self_type > >::type CastMaybe;

  CastMaybe cast_maybe;

  Swig::Director* l_director = cast_maybe.value(arg1);
  bool l_upcall = (l_director && (l_director->swig_get_self()==obj0));

  $result = SP_SiconosVector_out($1, l_upcall, obj0);
}

// director output : PyObject -> SP::SiconosVector
%typemap(directorout, fragment="SiconosVector") std11::shared_ptr<SiconosVector> ()
{
  // %typemap(directorout, fragment="SiconosVector") std11::shared_ptr<SiconosVector> ()
  c_result = SiconosVector_from_python($input);
  if (!c_result) throw Swig::DirectorMethodException();
}

//////////////////////////////////////////////////////////////////////////////

// director input : TYPE -> numpy
%typemap(out, fragment="SiconosVector") SiconosVector & ()
{
  // %typemap(out, fragment="SiconosVector") SiconosVector& ()
  $result = SiconosVector_to_numpy(*$1);
}

// director input : TYPE -> numpy
%typemap(directorin, fragment="SiconosVector") SiconosVector & ()
{
  //%typemap(directorin, fragment="SiconosVector") SiconosVector & ()
  $input = SiconosVector_to_numpy($1_name);
}

//%typemap(directorout, fragment="NumPy_Fragments") (SiconosMatrix&) ()
//{
//  // %typemap(directorout, fragment="NumPy_Fragments") SiconosMatrix& ()
//  void * swig_argp;
//  int swig_res = SWIG_ConvertPtr(result,&swig_argp,SWIGTYPE_p_SiconosMatrix,  0  | 0);
//
//  if (!SWIG_IsOK(swig_res))
//  {
//    // try a conversion from numpy
//    PyArrayObject* array = NULL;
//    int is_new_object = 0;
//    array = obj_to_array_fortran_allow_conversion($input, NPY_DOUBLE,&is_new_object);
//    if (!require_dimensions(array,2) ||
//        !require_native(array) || !require_fortran(array)) throw Swig::DirectorMethodException();
//
//
//    SP::SimpleMatrix tmp;
//    tmp.reset(new SimpleMatrix(array_size(array,0), array_size(array,1)));
//    // copy this is due to SimpleMatrix based on resizable std::vector
//    memcpy(&*tmp->getArray(),array_data(array),array_size(array,0)*array_size(array,1)*sizeof(double));
//    return *tmp;
//  }
//
//  if (!swig_argp)
//  {
//    return (SiconosMatrix) c_result;
//  }
//  else
//  {
//    c_result = *(reinterpret_cast< SP::SiconosMatrix * >(swig_argp));
//    if (SWIG_IsNewObj(swig_res)) delete reinterpret_cast< SP::SiconosMatrix * >(swig_argp);
//    return (SP::SiconosMatrix) c_result;
//  }
//}
//

//%typemap(directorin) VectorOfSMatrices& ()
//{
//  // %typemap(directorin) (VectorOfSMatrices&) ()
//  // swig issue shared pointer check in wrappers even if arg is a ref
//  SP::VectorOfSMatrices myptemp(createSPtrVectorOfSMatrices($1));
//  $input = SWIG_NewPointerObj(SWIG_as_voidptr(&myptemp),
//                              SWIGTYPE_p_std__vectorT_std11__shared_ptrT_SimpleMatrix_t_std__allocatorT_std11__shared_ptrT_SimpleMatrix_t_t_t, 0);
//}
//



//////////////////////////////////////////////////////////////////////////////
%define TYPEMAP_MATRIX(TYPE)
// numpy or TYPE on input -> TYPE
%typemap(in, fragment="SiconosMatrix") (std11::shared_ptr<TYPE>) (PyArrayObject* array=NULL, int is_new_object = 0)
{

  void *argp1=0;
  int res1=0;
  int newmem = 0;
  SP::TYPE tempshared1 ;
  SP::TYPE *smartarg1 = 0 ;

   // try a conversion from a SiconosMatrix
  res1 = SWIG_ConvertPtrAndOwn($input, &argp1, $descriptor(SP::TYPE *), 0 |  0 , &newmem);
  if (SWIG_IsOK(res1) && argp1)
  {
    if (newmem & SWIG_CAST_NEW_MEMORY)
    {
      // taken from generated code
      tempshared1 = *reinterpret_cast< SP::TYPE * >(argp1);
      delete reinterpret_cast< SP::TYPE * >(argp1);
      $1 = tempshared1;
    }
    else {
      smartarg1 = reinterpret_cast< SP::TYPE * >(argp1);
      $1 = *smartarg1;
    }
  }
  else
  {
    bool ok = SiconosMatrix_from_python($input, &array, &is_new_object, &$1);
    if (!ok) { SWIG_exception_fail(SWIG_ValueError, "expected matrix"); }
  }
}

%typemap(in, fragment="SiconosMatrix")
  const TYPE &
  (PyArrayObject* array = NULL, int is_new_object = 0, std::vector<SP::TYPE> keeper)
{
   bool ok = SiconosMatrix_from_python($input, &array, &is_new_object, &$1, keeper);
   if (!ok)
   {
     void* swig_argp=0;
     int swig_res = SWIG_ConvertPtr($input, &swig_argp, $descriptor(SP::TYPE *), 0);
     if (SWIG_IsOK(swig_res) && swig_argp)
     {
       $1 = (%reinterpret_cast(swig_argp, SP::TYPE *))->get();
     }
     else { SWIG_exception_fail(SWIG_ValueError, "expected matrix"); }
   }
   if (!$1) { SWIG_exception_fail(SWIG_ValueError, "expected matrix"); }
}

%typemap(in, fragment="SiconosMatrix")
  TYPE&
  (PyArrayObject* array = NULL, int is_new_object = 0, std::vector<SP::TYPE> keeper)
{
   bool ok = SiconosMatrix_from_python($input, &array, &is_new_object, &$1, keeper);
   if (!ok)
   {
     void* swig_argp=0;
     int swig_res = SWIG_ConvertPtr($input, &swig_argp, $descriptor(SP::TYPE *), 0);
     if (SWIG_IsOK(swig_res) && swig_argp)
     {
       $1 = (%reinterpret_cast(swig_argp, SP::TYPE *))->get();
     }
     else { SWIG_exception_fail(SWIG_ValueError, "expected matrix"); }
   }
   if (!$1) { SWIG_exception_fail(SWIG_ValueError, "expected matrix"); }
}

%typemap(freearg) (std11::shared_ptr<TYPE>)
{
  if (is_new_object$argnum && array$argnum)
    { Py_DECREF(array$argnum); }
}

%typemap(freearg) (TYPE)
{
  if (is_new_object$argnum && array$argnum)
    { Py_DECREF(array$argnum); }
}

%typemap(freearg) (TYPE&)
{
  if (is_new_object$argnum && array$argnum)
    { Py_DECREF(array$argnum); }
}

// director input : TYPE -> numpy
%typemap(directorin, fragment="SiconosMatrix") (std11::shared_ptr<TYPE>) ()
{
  // %typemap(directorin, fragment="SiconosMatrix") std11::shared_ptr<TYPE> ()
  $input = SiconosMatrix_to_numpy($1_name);
 }

// director input : TYPE -> numpy
%typemap(directorin, fragment="SiconosMatrix") TYPE& ()
{
  // %typemap(directorin, fragment="SiconosMatrix") TYPE& ()
  $input = SiconosMatrix_to_numpy($1_name);
 }

// director output : PyObject -> SP::SiconosMatrix
%typemap(directorout, fragment="SiconosMatrix") std11::shared_ptr<TYPE> ()
{
  // %typemap(directorout, fragment="NumPy_Fragments") std11::shared_ptr<SiconosMatrix> ()
  PyArrayObject* array_dout = NULL;
  int is_new_object_dout;
  bool ok = SiconosMatrix_from_python($input, &array_dout, &is_new_object_dout, &c_result);
  if (!ok) throw Swig::DirectorMethodException();
}

%typemap(out, fragment="SiconosMatrix") (std11::shared_ptr<TYPE>)
{
  // %typemap(out) TYPE

  // compile time test to reproduce swig director test. swig does not
  // seem to provides facilities to customize this
  // arg1 is the class instantiation (swig 2.0) => no way to get this as a swig
  // variable.
  typedef BOOST_TYPEOF(arg1) self_type;

  typedef boost::mpl::eval_if<boost::is_polymorphic<self_type >,
    DirectorCast<self_type >,
    DirectorNoCast<self_type > >::type CastMaybe;

  CastMaybe cast_maybe;

  Swig::Director* l_director = cast_maybe.value(arg1);
  bool l_upcall = (l_director && (l_director->swig_get_self()==obj0));

  // call from director?
  $result = SiconosMatrix_to_numpy($1, l_upcall);
}
%enddef
//////////////////////////////////////////////////////////////////////////////

%inline %{
  // SWIG_DIRECTOR_CAST cannot be use on non polymorphic data type
  // (i.e. classes without at least one virtual method) so we have to
  // switch at compile time. This is done with boost::mpl::eval_if

  template<typename T>
    struct DirectorCast
  {
    typedef DirectorCast type;
    Swig::Director* value(T& p)
    {
      return SWIG_DIRECTOR_CAST(p);
    }
  };

  template<typename T>
    struct DirectorNoCast
  {
    typedef DirectorNoCast type;
    Swig::Director* value(T& p)
    {
      return 0;
    }
  };
%}


%typemap(out) std11::shared_ptr<std::vector<unsigned int> >
{
  // %typemap(out) std11::shared_ptr<SiconosVector>

  // compile time test to reproduce swig director test. swig does not
  // seem to provides facilities to customize this
  // arg1 is the class instantiation (swig 2.0) => no way to get this as a swig
  // variable.
  typedef BOOST_TYPEOF(arg1) self_type;

  typedef boost::mpl::eval_if<boost::is_polymorphic<self_type >,
    DirectorCast<self_type >,
    DirectorNoCast<self_type > >::type CastMaybe;

  CastMaybe cast_maybe;

  Swig::Director* l_director = cast_maybe.value(arg1);
  bool l_upcall = (l_director && (l_director->swig_get_self()==obj0));

  // call from director?
  if (l_upcall)
  {
    // result from C++ method, return the pointer
    $result = SWIG_NewPointerObj(SWIG_as_voidptr(&$1), SWIGTYPE_p_std11__shared_ptrT_std__vectorT_unsigned_int_std__allocatorT_unsigned_int_t_t_const_t,  0 );
  }
  // call from python : return numpy from SiconosVector
  else
  {
    if ($1)
    {
      // /!\ need check for a dense vector!

      npy_intp this_vector_dim[1];
      this_vector_dim[0]=$1->size();

      PYARRAY_FROM_SHARED_STL_VECTOR(NPY_UINT,1,this_vector_dim,$1,$result);
    }
    else
    {
      Py_INCREF(Py_None);
      $result = Py_None;
    }
  }
}

// check on input : a python sequence
%typecheck(SWIG_TYPECHECK_INTEGER) (std11::shared_ptr<std::vector<unsigned int> >)
%{
  $1 = PySequence_Check($input);
%}

// python int sequence => std::vector<unsigned int>
%fragment("sequenceToUnsignedIntVector","header",fragment="SWIG_AsVal_int")
{
  static inline int sequenceToUnsignedIntVector(
    PyObject *input,
    std11::shared_ptr<std::vector<unsigned int> >& ptr)
  {
    // Return success but leave ptr as nullptr if input is None
    if (input == Py_None) {
      ptr.reset();
      return 1;
    }

    if (!PySequence_Check(input)) {
      PyErr_SetString(PyExc_TypeError,"Expecting a sequence");
      return 0;
    }

    ptr.reset(new std::vector<unsigned int>());
    assert(ptr);

    unsigned int i = 0;
    PyObject *o = NULL;
    for (i = 0; i < PyObject_Length(input); i++)
    {
      o = PySequence_GetItem(input,i);
      int x;
      int ecode = SWIG_AsVal_int(o, &x);
      if (!SWIG_IsOK(ecode) && PyArray_CheckScalar(o)) {
        // TODO: By treating as int, we lose half the range of the
        // unsigned int!  Okay for current uses of unsigned int in
        // Siconos.
        x = PyArray_PyIntAsInt(o);
      }
      if (x < 0)
        break;
      ptr->push_back((unsigned int) x);
      Py_DECREF(o);
    }
    if (i < PyObject_Length(input)) {
      Py_XDECREF(o);
      PyErr_SetString(PyExc_ValueError,"Expecting a sequence of unsigned int");
      return 0;
    }
    return 1;
  }
}

// int sequence => std::vector<unsigned int>
%typemap(in,fragment="sequenceToUnsignedIntVector") std11::shared_ptr<std::vector<unsigned int> > (std11::shared_ptr<std::vector<unsigned int> > temp)
{
  if (!sequenceToUnsignedIntVector($input, $1))
  {
    SWIG_fail;
  }
}

//////////////////////////
// VectorOf..
%template(VectorOfVectors) std::vector< std11::shared_ptr<SiconosVector> >;
%template(VectorOfBlockVectors) std::vector< std11::shared_ptr<BlockVector> >;
%template(VectorOfMatrices) std::vector< std11::shared_ptr<SiconosMatrix> >;
%template(VectorOfSMatrices) std::vector< std11::shared_ptr<SimpleMatrix> >;
%template(VectorOfMemories) std::vector< std11::shared_ptr<SiconosMemory> >;

// Other vector types
%template(UnsignedIntVector) std::vector<unsigned int>;

//////////////////////////

TYPECHECK(std11::shared_ptr<SiconosVector>);
TYPECHECK(std11::shared_ptr<SiconosMatrix>);
TYPECHECK(std11::shared_ptr<SimpleMatrix>);
TYPECHECK_REF(SiconosVector);
TYPECHECK_REF(SiconosMatrix);
TYPECHECK_REF(SimpleMatrix);

TYPEMAP_MATRIX(SiconosMatrix);
TYPEMAP_MATRIX(SimpleMatrix);

%apply (std11::shared_ptr<SiconosVector>) { (SP::SiconosVector) };
//%apply (SiconosVector &) { (const SiconosVector &) };

%apply (std11::shared_ptr<SiconosMatrix>) { (SP::SiconosMatrix) };
%apply (std11::shared_ptr<SimpleMatrix>) { (SP::SimpleMatrix) };

//%apply (SiconosMatrix &) { (SiconosMatrix const &) };
//%apply (SimpleMatrix &) { (SimpleMatrix const &) };
//%apply (SiconosMatrix &) { (const SiconosMatrix &) };
//%apply (SimpleMatrix &) { (const SimpleMatrix &) };

%apply (std11::shared_ptr<std::vector<unsigned int> >) { (SP::UnsignedIntVector) };

// cast to get the right class in Python

%define %_factory_dispatch_SP(TYPE)
if (!dcast) {
  SP::TYPE dobj = std11::dynamic_pointer_cast<TYPE>(ds);
  if (dobj) {
    dcast = 1;
    SP::TYPE *dobjp = new SP::TYPE(dobj);
    *_tmpobj = SWIG_NewPointerObj(%as_voidptr(dobjp),$descriptor(SP::TYPE *), SWIG_POINTER_OWN);
  }
}%enddef

%define %convert_sp_ds(mds, resultobj)
  int dcast = 0;
  SP::DynamicalSystem ds = mds;
  PyObject** _tmpobj = &resultobj;
  %formacro(%_factory_dispatch_SP, FirstOrderLinearTIDS, FirstOrderLinearDS, FirstOrderNonLinearDS)
  %formacro(%_factory_dispatch_SP, LagrangianLinearTIDS, LagrangianDS)
  %formacro(%_factory_dispatch_SP, NewtonEulerDS)
  if (!dcast) {
    SP::DynamicalSystem *pds = new SP::DynamicalSystem(ds);
    resultobj = SWIG_NewPointerObj(%as_voidptr(pds), $descriptor(SP::DynamicalSystem *), SWIG_POINTER_OWN);
  }
%enddef

// from factory.swg, but with dynamic_pointer_cast instead of dynamic_cast
%define %_factory_dispatch_SP_output(TYPE)
if (!dcast) {
  SP::TYPE dobj = std11::dynamic_pointer_cast<TYPE>($1);
  if (dobj) {
    dcast = 1;
    SP::TYPE *dobjp = new SP::TYPE(dobj);
    %set_output(SWIG_NewPointerObj(%as_voidptr(dobjp),$descriptor(SP::TYPE *), SWIG_POINTER_OWN));
  }
}%enddef

%typemap(out) SP::DynamicalSystem
{
  int dcast = 0;
  %formacro(%_factory_dispatch_SP_output, FirstOrderLinearTIDS, FirstOrderLinearDS, FirstOrderNonLinearDS)
  %formacro(%_factory_dispatch_SP_output, LagrangianLinearTIDS, LagrangianDS)
  %formacro(%_factory_dispatch_SP_output, NewtonEulerDS)
  if (!dcast) {
    %set_output(SWIG_NewPointerObj(%as_voidptr(&$1),$descriptor(SP::DynamicalSystem *), $owner | %newpointer_flags));
  }
}

%typemap(out) SP::Simulation
{
   int dcast = 0;
   %formacro(%_factory_dispatch_SP_output, TimeSteppingCombinedProjection, TimeSteppingDirectProjection, TimeStepping, TimeSteppingD1Minus, EventDriven)
   if (!dcast) {
      %set_output(SWIG_NewPointerObj(%as_voidptr(&$1),$descriptor(SP::Type *), $owner | %newpointer_flags));
   }
}

%typemap(out) SP::Relation
{
   int dcast = 0;
   %formacro(%_factory_dispatch_SP_output, FirstOrderLinearTIR, FirstOrderLinearR, FirstOrderType1R, FirstOrderType2R, FirstOrderNonLinearR, FirstOrderR)
   %formacro(%_factory_dispatch_SP_output, LagrangianLinearTIR, LagrangianCompliantR, LagrangianRheonomousR, LagrangianScleronomousR, LagrangianR)
   %formacro(%_factory_dispatch_SP_output, NewtonEulerFrom3DLocalFrameR, NewtonEulerFrom1DLocalFrameR, NewtonEulerR)
   if (!dcast) {
      %set_output(SWIG_NewPointerObj(%as_voidptr(&$1),$descriptor(SP::Type *), $owner | %newpointer_flags));
   }
}
