// signatures
%feature("autodoc", 2);

// named parameters (broken with overloaded function)
// %feature("kwargs");

// head declarations

// to avoid name conflicts
%{
#define SWIG_FILE_WITH_INIT

#ifdef __cplusplus
extern "C"
#endif

SWIGEXPORT
#if PY_VERSION_HEX >= 0x03000000
PyObject*
#else
void
#endif
SWIG_init(void);

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
#include <debug.h>
#include <numerics_verbose.h>

#if PY_VERSION_HEX < 0x02070000
#ifndef PyBytes_Check
#define PyBytes_Check PyString_Check
#define PyBytes_Size PyString_Size
#define PyBytes_AsString PyString_AsString
#define PyBytes_AS_STRING PyString_AS_STRING
#define PyBytes_GET_SIZE PyString_GET_SIZE
#endif

#ifndef PyCapsule_New
#define PyCapsule_New PyCObject_FromVoidPtrAndDesc
#define PyCapsule_CheckExact PyCObject_Check
#define PyCapsule_GetPointer(o, n) PyCObject_GetDesc((o))
#endif

#endif
%}

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

%include stl.i

#endif
