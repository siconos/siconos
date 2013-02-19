

// SWIG interface for Siconos IO
%module(directors="1", allprotected="1") IO

%import Kernel.i

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
#include <SiconosVisitor.hpp>
#include "SiconosPointers.hpp"
#include <SiconosRestart.hpp>
%}

%include "FrontEndConfig.h";

// common declarations

%include Common.i

// mandatory !
%rename (lambda_) lambda;

// shared ptr management
#define SWIG_SHARED_PTR_NAMESPACE std11
%include "boost_shared_ptr.i"

// handle standard exceptions
%{
 static void handle_exception(void) {
    try {
      throw;
    } 
    catch (const std::invalid_argument& e)
    {
      PyErr_SetString(PyExc_ValueError, e.what());
    }
    catch (const std::out_of_range& e)
    {
      PyErr_SetString(PyExc_IndexError, e.what());
    }
    catch (const SiconosException& e)
    {
      PyErr_SetString(PyExc_Exception, e.report().c_str());
    }
    catch (const Swig::DirectorException& e)
    {
      PyErr_SetString(PyExc_ValueError, e.getMessage());
    }
  }
%} 


%include "exception.i"
%exception
{
  try
  {
    $action;
  }
  catch (...) {
    if (!PyErr_Occurred()) {
      handle_exception();
    }
    SWIG_fail;
  } 
}

%include "SiconosRestart.hpp"
