#include <math.h>
void FGyrEffects(double* FGyr, const double* Q, const double* QP)
{

double Q1 = Q[0];
double QP1 = QP[0];
double Q2 = Q[1];
double QP2 = QP[1];

double tmp0 = sin(Q2);
double tmp1 = 10.0*sin(Q1 + Q2);
double tmp2 = 1.0*tmp0;

FGyr[0]=-2.0*QP1*QP2*tmp0 - pow(QP2, 2)*tmp2 + tmp1 + 20.0*sin(Q1);
FGyr[1]=pow(QP1, 2)*tmp2 + tmp1;

return;

}

