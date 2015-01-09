#include <math.h>
void JacobianVFGyr(double* Jac, const double* Q, const double* QP)
{

double Q1 = Q[0];
double QP1 = QP[0];
double Q2 = Q[1];
double QP2 = QP[1];

double tmp0 = 2.0*sin(Q2);
double tmp1 = -QP2*tmp0;
double tmp2 = QP1*tmp0;

Jac[0]=tmp1;
Jac[1]=tmp1 - tmp2;
Jac[2]=tmp2;
Jac[3]=0;

return;

}

