#include <math.h>
void Inertia(double* I, const double* Q)
{

double Q1 = Q[0];
double Q2 = Q[1];

double tmp0 = cos(Q2);
double tmp1 = 1.0*tmp0 + 1.0;

I[0]=2.0*tmp0 + 3.0;
I[1]=tmp1;
I[2]=tmp1;
I[3]=1.00000000000000;

return;

}

