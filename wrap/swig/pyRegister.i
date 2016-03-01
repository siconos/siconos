// -*- c++ -*-

#undef PY_REGISTER
#undef PY_FULL_REGISTER
#undef PY_REGISTER_UNSHARED
#undef PY_REGISTER_WITHOUT_DIRECTOR
#undef PY_REGISTER_WITHOUT_DIRECTOR_REF
#undef PY_REGISTER_WITHOUT_DIRECTOR_REF_ONLY
#undef PY_REGISTER_SIMPLEMATRIX


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
%typemap(directorin) (TYPE&) ()
%{
  // %typemap(directorin) (TYPE&) ()
  // swig issue shared pointer check in wrappers even if arg is a ref
  {
    SP::TYPE myptemp = createSPtr##TYPE($1);
      $input = SWIG_NewPointerObj(SWIG_as_voidptr(&myptemp),
                                  SWIGTYPE_p_std11__shared_ptrT_##TYPE##_t, 0);
  }
%}

%typemap(directorin) (const TYPE&) ()
%{
  // %typemap(directorin) (const TYPE&) ()
  // swig issue shared pointer check in wrappers even if arg is a ref
  {
    SPC::TYPE myptemp = createSPtrConst##TYPE($1);
      $input = SWIG_NewPointerObj(SWIG_as_voidptr(&myptemp),
                                  SWIGTYPE_p_std11__shared_ptrT_##TYPE##_const_t, 0);
  }
%}

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
%typemap(directorin) (TYPE&) ()
%{
  // %typemap(directorin) (TYPE&) ()
  // swig issue shared pointer check in wrappers even if arg is a ref
  {
    SP::TYPE myptemp = createSPtr##TYPE($1);
      $input = SWIG_NewPointerObj(SWIG_as_voidptr(&myptemp),
                                  SWIGTYPE_p_std11__shared_ptrT_##TYPE##_t, 0);
  }
%}

%typemap(directorin) (const TYPE&) ()
%{
  // %typemap(directorin) (const TYPE&) ()
  // swig issue shared pointer check in wrappers even if arg is a ref
  {
    SPC::TYPE myptemp = createSPtrConst##TYPE($1);
      $input = SWIG_NewPointerObj(SWIG_as_voidptr(&myptemp),
                                  SWIGTYPE_p_std11__shared_ptrT_##TYPE##_const_t, 0);
  }
%}

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
// swig does not apply properly the typemap
// TODO fix this
%typemap(in, fragment="SiconosMatrix")
  const TYPE &
  (PyArrayObject* array = NULL, int is_new_object, std::vector<SP::TYPE> keeper)
{
   bool ok = SiconosMatrix_from_python($input, &array, &is_new_object, &$1, keeper);
   if (!ok) SWIG_fail;
}

// swig does not apply properly the typemap
// TODO fix this
%typemap(in, fragment="SiconosMatrix") 
  TYPE &
  (PyArrayObject* array = NULL, int is_new_object, std::vector<SP::TYPE> keeper)
{
   bool ok = SiconosMatrix_from_python($input, &array, &is_new_object, &$1, keeper);
   if (!ok) SWIG_fail;
}

%typemap(directorin) TYPE& ()
%{
  // %typemap(directorin) (TYPE&) ()
  // swig issue shared pointer check in wrappers even if arg is a ref
  {
    SP::TYPE myptemp(createSPtr##TYPE($1));
    $input = SWIG_NewPointerObj(SWIG_as_voidptr(&myptemp),
                                $descriptor(SP::TYPE *), 0);
  }
%}

%typemap(directorout) TYPE& ()
%{
  // %typemap(directorout) (TYPE&) ()
  // swig issue shared pointer check in wrappers even if arg is a ref
  {
    SP::TYPE myptemp(createSPtr##TYPE($1));
    $result = SWIG_NewPointerObj(SWIG_as_voidptr(&myptemp),
                                 $descriptor(SP::TYPE *), 0);
  }
%}

%enddef

%define PY_REGISTER_SIMPLEMATRIX(TYPE)
PY_REGISTER_WITHOUT_DIRECTOR(TYPE)
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY)
(const TYPE&)
{
  // %typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY) (const TYPE&)
  $1 = 1;
}
%typemap(in, fragment="SiconosMatrix")
  const TYPE &
  (PyArrayObject* array = NULL, int is_new_object, std::vector<SP::TYPE> keeper)
{
   bool ok = SiconosMatrix_from_python($input, &array, &is_new_object, &$1, keeper);
   if (!ok) SWIG_fail;
}

%typemap(in, fragment="SiconosMatrix") 
  TYPE &
  (PyArrayObject* array = NULL, int is_new_object, std::vector<SP::TYPE> keeper)
{
   bool ok = SiconosMatrix_from_python($input, &array, &is_new_object, &$1, keeper);
   if (!ok) SWIG_fail;
}

%enddef

%define PY_REGISTER_WITHOUT_DIRECTOR_REF(TYPE)
PY_REGISTER_WITHOUT_DIRECTOR(TYPE)

// swig does not apply properly the typemap
// TODO fix this
%typemap(in,fragment="SiconosVector")
  const SiconosVector &
  (PyArrayObject* array=NULL, int is_new_object, std::vector<SP::SiconosVector> keeper)
{
  // %typemap(in,fragment="NumPy_Fragments")
  // %TYPE (PyArrayObject* array=NULL, int
  // %is_new_object)
  $1 = SiconosVector_in($input, &array, &is_new_object, keeper);
  if (!$1) SWIG_fail;
}

// swig does not apply properly the typemap
// TODO fix this
%typemap(in,fragment="SiconosVector") SiconosVector &
(PyArrayObject* array=NULL, int is_new_object, std::vector<SP::SiconosVector> keeper)
{
  // %typemap(in,fragment="NumPy_Fragments")
  // %TYPE (PyArrayObject* array=NULL, int
  // %is_new_object)
  $1 = SiconosVector_in($input, &array, &is_new_object, keeper);
  if (!$1) SWIG_fail;
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

