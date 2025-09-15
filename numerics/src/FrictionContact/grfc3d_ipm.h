
#ifndef GLOBALFRICTIONCONTACT3DSOLVERS_H
#define GLOBALFRICTIONCONTACT3DSOLVERS_H

#include "Friction_cst.h"
#include "GlobalRollingFrictionContactProblem.h"
#include "SolverOptions.h"
#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C" {
#endif

void grfc3d_IPM(GlobalRollingFrictionContactProblem* restrict problem,
                double* restrict reaction, double* restrict velocity,
                double* restrict globalVelocity, int* restrict info,
                SolverOptions* restrict options);
/** @} */

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
