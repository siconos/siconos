// -*- c++ -*-

#undef PY_REGISTER
#undef PY_FULL_REGISTER
#undef PY_REGISTER_UNSHARED
#undef PY_REGISTER_WITHOUT_DIRECTOR


// without include
%define PY_REGISTER_WITHOUT_HEADER(TYPE)
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
%typemap(directorin) (TYPE&) ()
%{
  // %typemap(directorin) (TYPE&) ()
  // swig issue shared pointer check in wrappers even if arg is a ref
  SP::TYPE myptemp = createSPtr##TYPE($1);
  $input = SWIG_NewPointerObj(SWIG_as_voidptr(&myptemp),
                              SWIGTYPE_p_std11__shared_ptrT_##TYPE##_t, 0);
%}

%typemap(directorin) (const TYPE&) ()
%{
  // %typemap(directorin) (const TYPE&) ()
  // swig issue shared pointer check in wrappers even if arg is a ref
  SPC::TYPE myptemp = createSPtrConst##TYPE($1);
  $input = SWIG_NewPointerObj(SWIG_as_voidptr(&myptemp),
                              SWIGTYPE_p_std11__shared_ptrT_##TYPE##_const_t, 0);
%}

%shared_ptr(TYPE);
%make_picklable(TYPE, Kernel);
REF_PTR(TYPE)
%enddef

%define PY_REGISTER(TYPE)
%inline
%{
#include <TYPE.hpp>
%}
PY_REGISTER_WITHOUT_HEADER(TYPE)
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
%typemap(directorin) (TYPE&) ()
%{
  // %typemap(directorin) (TYPE&) ()
  // swig issue shared pointer check in wrappers even if arg is a ref
  SP::TYPE myptemp = createSPtr##TYPE($1);
  $input = SWIG_NewPointerObj(SWIG_as_voidptr(&myptemp),
                              SWIGTYPE_p_std11__shared_ptrT_##TYPE##_t, 0);
%}

%typemap(directorin) (const TYPE&) ()
%{
  // %typemap(directorin) (const TYPE&) ()
  // swig issue shared pointer check in wrappers even if arg is a ref
  SPC::TYPE myptemp = createSPtrConst##TYPE($1);
  $input = SWIG_NewPointerObj(SWIG_as_voidptr(&myptemp),
                              SWIGTYPE_p_std11__shared_ptrT_##TYPE##_const_t, 0);
%}

%shared_ptr(TYPE);
%include TYPE.hpp
%make_picklable(TYPE, Kernel);
REF_PTR(TYPE)
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
REF_PTR(TYPE)
%enddef
