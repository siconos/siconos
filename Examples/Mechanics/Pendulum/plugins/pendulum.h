#include <stdlib.h>
#include <math.h>

extern "C" void Inertia(double* Inertia, const double* Q);

extern "C" void JacobianQFGyr(double* Jac, const double* Q, const double* QP);

extern "C" void JacobianVFGyr(double* Jac, const double* Q, const double* QP);

extern "C" void FGyrEffects(double* FGyr, const double* Q, const double* QP);


