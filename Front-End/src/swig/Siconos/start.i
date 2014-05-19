// named parameters
%feature("kwargs");


// head declarations

// to avoid name conflicts
%{
#define SWIG_FILE_WITH_INIT
#include <sstream>
#if defined(Py_COMPLEXOBJECT_H)
#undef c_sum
#undef c_diff
#undef c_neg
#undef c_prod
#undef c_quot
#undef c_pow
#undef c_abs
#endif
#include "FrontEndConfig.h"

#include <assert.h>
//#define DEBUG_MESSAGES 1
#include <debug.h>
#include <misc.h>
%}

%include "FrontEndConfig.h"

// numpy macros
%include numpy.i 	

%include <misc.h>

%init %{
  import_array();
%}

%inline %{

#define FPyArray_SimpleNewFromData(nd, dims, typenum, data)             \
  PyArray_New(&PyArray_Type, nd, dims, typenum, NULL,                   \
              data, 0, NPY_ARRAY_FARRAY, NULL)

%}
