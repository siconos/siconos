// -*- c++ -*-

// without include
%define PY_REGISTER(TYPE)
%inline
%{
#include <TYPE.hpp>
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
%feature("director") TYPE;
%ignore STD11::enable_shared_from_this<TYPE>;
%shared_ptr(STD11::enable_shared_from_this<TYPE>); // warning 520 suppression
%template (shared ## TYPE) STD11::enable_shared_from_this<TYPE>;
%typemap(directorin) TYPE& ()
{
  // %typemap(directorin) (TYPE&) ()
  // swig issue shared pointer check in wrappers even if arg is a ref
  SP::TYPE myptemp(createSPtr##TYPE($1));
  $input = SWIG_NewPointerObj(SWIG_as_voidptr(&myptemp),
                              SWIGTYPE_p_std11__shared_ptrT_##TYPE##_t, 0);
}
%typemap(directorin) const TYPE& ()
{
  // %typemap(directorin) (TYPE&) ()
  // swig issue shared pointer check in wrappers even if arg is a ref
  SPC::TYPE myptemp(createSPtrConst##TYPE($1));
  $input = SWIG_NewPointerObj(SWIG_as_voidptr(&myptemp),
                              SWIGTYPE_p_std11__shared_ptrT_##TYPE##_const_t, 0);
}
%shared_ptr(TYPE);
%make_picklable(TYPE, Kernel);
%enddef

%define PY_FULL_REGISTER(TYPE)
%inline
%{
#include <TYPE.hpp>
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
%feature("director") TYPE;
%ignore STD11::enable_shared_from_this<TYPE>;
%shared_ptr(STD11::enable_shared_from_this<TYPE>); // warning 520 suppression
%template (shared ## TYPE) STD11::enable_shared_from_this<TYPE>;
%typemap(directorin) TYPE& ()
{
  // %typemap(directorin) (TYPE&) ()
  // swig issue shared pointer check in wrappers even if arg is a ref
  SP::TYPE myptemp(createSPtr##TYPE($1));
  $input = SWIG_NewPointerObj(SWIG_as_voidptr(&myptemp),
                              SWIGTYPE_p_std11__shared_ptrT_##TYPE##_t, 0);
}
%typemap(directorin) const TYPE& ()
{
  // %typemap(directorin) (TYPE&) ()
  // swig issue shared pointer check in wrappers even if arg is a ref
  SP::TYPE myptemp(createSPtrConst##TYPE($1));
  $input = SWIG_NewPointerObj(SWIG_as_voidptr(&myptemp),
                              SWIGTYPE_p_std11__shared_ptrT_##TYPE##_const_t, 0);
}
%shared_ptr(TYPE);
%include TYPE.hpp
%make_picklable(TYPE, Kernel);
%enddef

%define PY_REGISTER_WITHOUT_DIRECTOR_REF_ONLY(TYPE)
PY_REGISTER_WITHOUT_DIRECTOR(TYPE)
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY)
(const TYPE&)
{
  // %typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY) (const TYPE&)
  $1 = 1;
}
%typemap(directorin) TYPE& ()
{
  // %typemap(directorin) (TYPE&) ()
  // swig issue shared pointer check in wrappers even if arg is a ref
  SP::TYPE myptemp(createSPtr##TYPE($1));
  $input = SWIG_NewPointerObj(SWIG_as_voidptr(&myptemp),
                              SWIGTYPE_p_std11__shared_ptrT_##TYPE##_t, 0);
}

%typemap(directorout) TYPE& ()
{
  // %typemap(directorout) (TYPE&) ()
  // swig issue shared pointer check in wrappers even if arg is a ref
  SP::TYPE myptemp(createSPtr##TYPE($1));
  $result = SWIG_NewPointerObj(SWIG_as_voidptr(&myptemp),
                              SWIGTYPE_p_std11__shared_ptrT_##TYPE##_t, 0);
}


%enddef

%define PY_REGISTER_SIMPLEMATRIX(TYPE)
PY_REGISTER_WITHOUT_DIRECTOR(TYPE)
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY)
(const TYPE&)
{
  // %typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY) (const TYPE&)
  $1 = 1;
}
// numpy or TYPE on input -> TYPE 
%typemap(in,fragment="NumPy_Fragments") 
  TYPE &
(PyArrayObject* array=NULL, int is_new_object, std::vector<SP::SimpleMatrix> keeper)
{
  // %typemap(in,fragment="NumPy_Fragments")
  // %TYPE (PyArrayObject* array=NULL, int
  // %is_new_object)
  void *argp1=0;
  int res1=0;
  int newmem = 0;
  TYPE temp1 ;
  TYPE *smartarg1 = 0 ;
 
  // try a conversion from TYPE
  res1 = SWIG_ConvertPtrAndOwn($input, &argp1, $descriptor(TYPE *), 0 |  0 , &newmem);
  if (SWIG_IsOK(res1)) 
  {
    if (newmem & SWIG_CAST_NEW_MEMORY) 
    {
      // taken from generated code
      temp1 = *reinterpret_cast< TYPE * >(argp1);
      delete reinterpret_cast< TYPE * >(argp1);
      $1 = &temp1;
    } 
    else {
      smartarg1 = reinterpret_cast< TYPE * >(argp1);
      $1 = smartarg1;
    }
  }
  else
  {
    array = obj_to_array_fortran_allow_conversion($input, NPY_DOUBLE,&is_new_object);

    if (!array)
    {
      void *argp;
      SWIG_fail; // not implemented : $1 = type_conv($input) (type check done above)
    }
    else
    {
      if (!require_dimensions(array,2) ||
          !require_native(array) || !require_fortran(array)) SWIG_fail;

      SP::SimpleMatrix tmp(new SimpleMatrix(array_size(array,0), array_size(array,1)));
      // copy : with SiconosVector based on resizable std::vector there is
      // no other way
      memcpy(&*tmp->getArray(),array_data(array),array_size(array,0)*array_size(array,1)*sizeof(double));
      $1 = &*tmp;
      keeper.push_back(tmp);
    }
  }
}
%enddef

%define PY_REGISTER_WITHOUT_DIRECTOR_REF_BASEMAT(TYPE)
PY_REGISTER_WITHOUT_DIRECTOR(TYPE)
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY)
(const TYPE&)
{
  // %typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY) (const TYPE&)
  $1 = 1;
}
%typemap(in,fragment="NumPy_Fragments") 
  const TYPE &
(PyArrayObject* array=NULL, int is_new_object, std::vector<SP::SiconosMatrix> keeper)
{
  // %typemap(in,fragment="NumPy_Fragments")
  // %TYPE (PyArrayObject* array=NULL, int
  // %is_new_object)
  void *argp1=0;
  int res1=0;
  int newmem = 0;
  TYPE temp1 ;
  TYPE *smartarg1 = 0 ;
 
  // try a conversion from TYPE
  res1 = SWIG_ConvertPtrAndOwn($input, &argp1, SWIGTYPE_p_std11__shared_ptrT_##TYPE##_t, 0 |  0 , &newmem);
  if (SWIG_IsOK(res1)) 
  {
    if (newmem & SWIG_CAST_NEW_MEMORY) 
    {
      // taken from generated code
      temp1 = *reinterpret_cast< TYPE * >(argp1);
      delete reinterpret_cast< TYPE * >(argp1);
      $1 = &temp1;
    } 
    else {
      smartarg1 = reinterpret_cast< TYPE * >(argp1);
      $1 = smartarg1;
    }
  }
  else
  {
    array = obj_to_array_fortran_allow_conversion($input, NPY_DOUBLE,&is_new_object);

    if (!array)
    {
      void *argp;
      SWIG_fail; // not implemented : $1 = type_conv($input) (type check done above)
    }
    else
    {
      if (!require_dimensions(array,1) ||
          !require_native(array) || !require_fortran(array)) SWIG_fail;
      
      SP::SimpleMatrix tmp(new SimpleMatrix(array_size(array,0), array_size(array,1)));
      // copy : with SiconosVector based on resizable std::vector there is
      // no other way
      memcpy(&*tmp->getArray(),array_data(array),array_size(array,0)*sizeof(double)*array_size(array,1));
      $1 = &*tmp;
      keeper.push_back(tmp);
    }
  }
}

%enddef

%define PY_REGISTER_WITHOUT_DIRECTOR_REF(TYPE)
PY_REGISTER_WITHOUT_DIRECTOR(TYPE)
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY)
(const TYPE&)
{
  // %typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY) (const TYPE&)
  $1 = 1;
}

// numpy or TYPE on input -> TYPE 
%typemap(in,fragment="NumPy_Fragments") 
  const TYPE &
(PyArrayObject* array=NULL, int is_new_object, std::vector<SP::SiconosVector> keeper)
{
  // %typemap(in,fragment="NumPy_Fragments")
  // %TYPE (PyArrayObject* array=NULL, int
  // %is_new_object)
  void *argp1=0;
  int res1=0;
  int newmem = 0;
  TYPE temp1 ;
  TYPE *smartarg1 = 0 ;
 
  // try a conversion from TYPE
  res1 = SWIG_ConvertPtrAndOwn($input, &argp1, SWIGTYPE_p_std11__shared_ptrT_##TYPE##_t, 0 |  0 , &newmem);
  if (SWIG_IsOK(res1)) 
  {
    if (newmem & SWIG_CAST_NEW_MEMORY) 
    {
      // taken from generated code
      temp1 = *reinterpret_cast< TYPE * >(argp1);
      delete reinterpret_cast< TYPE * >(argp1);
      $1 = &temp1;
    } 
    else {
      smartarg1 = reinterpret_cast< TYPE * >(argp1);
      $1 = smartarg1;
    }
  }
  else
  {
    array = obj_to_array_fortran_allow_conversion($input, NPY_DOUBLE,&is_new_object);

    if (!array)
    {
      void *argp;
      SWIG_fail; // not implemented : $1 = type_conv($input) (type check done above)
    }
    else
    {
      if (!require_dimensions(array,1) ||
          !require_native(array) || !require_fortran(array)) SWIG_fail;
      
      SP::SiconosVector tmp(new SiconosVector(array_size(array,0)));
      // copy : with SiconosVector based on resizable std::vector there is
      // no other way
      memcpy(&*tmp->getArray(),array_data(array),array_size(array,0)*sizeof(double));
      $1 = &*tmp;
      keeper.push_back(tmp);
    }
  }
}

// numpy or TYPE on input -> TYPE 
%typemap(in,fragment="NumPy_Fragments") 
  TYPE &
(PyArrayObject* array=NULL, int is_new_object, std::vector<SP::SiconosVector> keeper)
{
  // %typemap(in,fragment="NumPy_Fragments")
  // %TYPE (PyArrayObject* array=NULL, int
  // %is_new_object)
  void *argp1=0;
  int res1=0;
  int newmem = 0;
  TYPE temp1 ;
  TYPE *smartarg1 = 0 ;
 
  // try a conversion from TYPE
  res1 = SWIG_ConvertPtrAndOwn($input, &argp1, SWIGTYPE_p_std11__shared_ptrT_##TYPE##_t, 0 |  0 , &newmem);
  if (SWIG_IsOK(res1)) 
  {
    if (newmem & SWIG_CAST_NEW_MEMORY) 
    {
      // taken from generated code
      temp1 = *reinterpret_cast< TYPE * >(argp1);
      delete reinterpret_cast< TYPE * >(argp1);
      $1 = &temp1;
    } 
    else {
      smartarg1 = reinterpret_cast< TYPE * >(argp1);
      $1 = smartarg1;
    }
  }
  else
  {
    array = obj_to_array_fortran_allow_conversion($input, NPY_DOUBLE,&is_new_object);

    if (!array)
    {
      void *argp;
      SWIG_fail; // not implemented : $1 = type_conv($input) (type check done above)
    }
    else
    {
      if (!require_dimensions(array,1) ||
          !require_native(array) || !require_fortran(array)) SWIG_fail;
      
      SP::SiconosVector tmp(new SiconosVector(array_size(array,0)));
      // copy : with SiconosVector based on resizable std::vector there is
      // no other way
      memcpy(&*tmp->getArray(),array_data(array),array_size(array,0)*sizeof(double));
      $1 = &*tmp;
      keeper.push_back(tmp);
    }
  }
}

%typemap(freearg) (const TYPE &)
{
  if (is_new_object$argnum && array$argnum)
    { Py_DECREF(array$argnum); }
}
%enddef

%define PY_REGISTER_WITHOUT_DIRECTOR(TYPE)
%inline
%{
#include <TYPE.hpp>
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
%enddef

