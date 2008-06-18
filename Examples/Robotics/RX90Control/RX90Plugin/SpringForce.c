#include "LagrangianModel.h"
#include <math.h>
void SpringForce(S, q)
double S[6];
double q[6];
{
  double t1;
  double t10;
  double t2;
  double t3;
  double t5;
  double t8;
  {
    S[0] = 0.0;
    t1 = q[1];
    t2 = sin(t1);
    t3 = t2 * t2;
    t5 = cos(t1);
    t8 = pow(0.45 - 0.6E-1 * t5, 2.0);
    t10 = sqrt(0.36E-2 * t3 + t8);
    S[1] = -0.27E-1 * (20295.0 * t10 - 0.625775E4) * t2 / t10;
    S[2] = 0.0;
    S[3] = 0.0;
    S[4] = 0.0;
    S[5] = 0.0;
    return;
  }
}

