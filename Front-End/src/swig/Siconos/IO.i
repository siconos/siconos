

// SWIG interface for Siconos IO
%module(directors="1", allprotected="1") IO

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
#include <SiconosKernel.hpp>
#include <SiconosRestart.hpp>
%}

%include "FrontEndConfig.h";

%include HandleException.i

%import Kernel.i

// shared ptr management
#define SWIG_SHARED_PTR_NAMESPACE std11
%include "boost_shared_ptr.i"

%include "SiconosRestart.hpp"
