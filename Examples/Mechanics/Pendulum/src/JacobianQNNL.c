#include <math.h>
void JacobianQFGyr(double* Jac, const double* Q, const double* QP)
{

double Q1 = Q[0];
double QP1 = QP[0];
double Q2 = Q[1];
double QP2 = QP[1];

double tmp0 = 10.0*cos(Q1 + Q2);
double tmp1 = cos(Q2);
double tmp2 = 1.0*tmp1;

Jac[0]=tmp0 + 20.0*cos(Q1);
Jac[1]=-2.0*QP1*QP2*tmp1 - pow(QP2, 2)*tmp2 + tmp0;
Jac[2]=tmp0;
Jac[3]=pow(QP1, 2)*tmp2 + tmp0;

return;

}

