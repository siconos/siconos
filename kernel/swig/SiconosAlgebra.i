// The SiconosAlgebra classes are handled seperately from the other Kernel classes
// This is because we define typemaps for them. It look like swig forgets some typemaps
// after being told more about a class (more specifically after applying PY_REGISTER_WITHOUT_DIRECTOR)
// Hence, we declare them fully here, and just after we define the typemaps (note the %include KernelTypes.i at the end)


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

%include KernelTypes.i

%include SiconosMatrix.hpp
%include SimpleMatrix.hpp
%include SiconosVector.hpp
%include BlockVector.hpp
