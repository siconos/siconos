#include <stdlib.h>
#include <math.h>

extern "C" void Inertia(double* Inertia, const double* Q);

extern "C" void JacobianQNNL(double* Jac, const double* Q, const double* QP);

extern "C" void JacobianVNNL(double* Jac, const double* Q, const double* QP);

extern "C" void NNLEffects(double* NNL, const double* Q, const double* QP);


