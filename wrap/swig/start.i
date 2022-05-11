// signatures
// docstrings are generated from doxygen using swig option -doxygen.
// Do not use %feature("autodoc", 3) !!

// named parameters (broken with overloaded function)
// %feature("kwargs");

// head declarations

// to avoid name conflicts
%{
// Numpy stuff - https://numpy.org/doc/stable/reference/swig.interface-file.html
#define SWIG_FILE_WITH_INIT

#ifdef __cplusplus
extern "C"
#endif
SWIGEXPORT PyObject* SWIG_init(void);
 

#ifdef __cplusplus
#include <sstream>
#endif


#if defined(Py_COMPLEXOBJECT_H)
#undef c_sum
#undef c_diff
#undef c_neg
#undef c_prod
#undef c_quot
#undef c_pow
#undef c_abs
#endif

#include <assert.h>
//#define DEBUG_MESSAGES 1
#include "siconos_debug.h"
#include "SiconosConfig.h"
#include "numerics_verbose.h"

%}

%include "SiconosConfig.h"

%include target_datatypes.i

#ifdef SWIGPYTHON
// numpy macros
%include numpy.i

%init %{
  import_array();
%}

%inline %{

#define FPyArray_SimpleNewFromData(nd, dims, typenum, data)             \
  PyArray_New(&PyArray_Type, nd, dims, typenum, NULL,                   \
              data, 0, NPY_ARRAY_FARRAY, NULL)

%}

// mandatory !
%rename (lambda_) lambda;
#endif /* SWIGPYTHON */

#ifdef SWIGMATLAB
%include numpy_matlab.i
#endif /* SWIGMATLAB */

#ifdef __cplusplus

%include ignored_functions.i

 // swig / STL. http://www.swig.org/Doc4.0/Library.html#Library_stl_cpp_library
%include stl.i

%feature("doxygen:alias:rst") "\verbatim embed:rst"
%feature("doxygen:alias:endrst") "\endverbatim"

#endif
