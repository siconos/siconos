%{
#include "rolling_fc_Solvers.h"
#include "Friction_cst.h"
#include "RollingFrictionContactProblem.h"
#include "rolling_fc3d_compute_error.h"
#ifdef WITH_FCLIB
// avoid a conflict with old csparse.h in case fclib.h includes it
#define _CS_H
#include "fclib_interface.h"
#endif
%}

%include "RollingFrictionContactProblem.h"
#ifdef WITH_FCLIB
// avoid a conflict with old csparse.h in case fclib.h includes it
#define _CS_H
%include fclib_interface.h
#endif

%include "fclib_interface.h"
%include "rolling_fc_Solvers.h"
%include "Friction_cst.h"
%include "rolling_fc3d_compute_error.h"

