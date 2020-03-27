// -*- c++ -*-

#undef PY_REGISTER
#undef PY_FULL_REGISTER
#undef PY_REGISTER_UNSHARED
#undef PY_REGISTER_WITHOUT_DIRECTOR


// without include
%define PY_REGISTER_WITHOUT_HEADER(TYPE, COMPONENT)
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
%ignore std::enable_shared_from_this<TYPE>;
%shared_ptr(std::enable_shared_from_this<TYPE>); // warning 520 suppression
%template (shared ## TYPE) std::enable_shared_from_this<TYPE>;
%shared_ptr(TYPE);
FIX_DIRECTOR_TYPEMAPS(TYPE)
%make_picklable(TYPE, COMPONENT);
REF_PTR(TYPE)
%enddef

%define PY_REGISTER(TYPE, COMPONENT)
%inline
%{
#include <TYPE.hpp>
%}
PY_REGISTER_WITHOUT_HEADER(TYPE, COMPONENT)
%enddef

%define PY_FULL_REGISTER(TYPE, COMPONENT)
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
%ignore std::enable_shared_from_this<TYPE>;
%shared_ptr(std::enable_shared_from_this<TYPE>); // warning 520 suppression
%template (shared ## TYPE) std::enable_shared_from_this<TYPE>;
%shared_ptr(TYPE);
FIX_DIRECTOR_TYPEMAPS(TYPE)
%include TYPE.hpp
%make_picklable(TYPE, COMPONENT);
REF_PTR(TYPE)
%enddef

%define PY_REGISTER_WITHOUT_DIRECTOR(TYPE, COMPONENT)
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
%ignore std::enable_shared_from_this<TYPE>;
%shared_ptr(std::enable_shared_from_this<TYPE>); // warning 520 suppression
%template (shared ## TYPE) std::enable_shared_from_this<TYPE>;
%shared_ptr(TYPE);
%make_picklable(TYPE, COMPONENT);
REF_PTR(TYPE)
%enddef
