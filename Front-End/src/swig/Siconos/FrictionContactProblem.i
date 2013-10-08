// -*- c++ -*- 
%{
#include "FrictionContactProblem.h"
#ifdef WITH_FCLIB
#include "fclib_interface.h"
#endif
%}

%include "FrictionContactProblem.h"
#ifdef WITH_FCLIB
%include fclib_interface.h
#endif
