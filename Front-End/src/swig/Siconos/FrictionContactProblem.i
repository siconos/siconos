// -*- c++ -*- 
%{
#include "FrictionContactProblem.h"
#include "GlobalFrictionContactProblem.h"
#ifdef WITH_FCLIB
#include "fclib_interface.h"
#endif
%}

%include "FrictionContactProblem.h"
%include "GlobalFrictionContactProblem.h"
#ifdef WITH_FCLIB
%include fclib_interface.h
#endif
