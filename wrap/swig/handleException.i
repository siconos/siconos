%{
#include <SiconosException.hpp>
%}


%feature("director:except") {
  if ($error != NULL) {
    throw Swig::DirectorMethodException();
  }
 }

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
    catch (const siconos::exception& e)
    {
      auto strout = boost::diagnostic_information(e, true);
      PyErr_SetString(PyExc_Exception, strout.data());
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
