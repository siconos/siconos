#include <math.h>
void NLEffects(N, q, qdot)
double N[2];
double q[2];
double qdot[2];
{
  double t1;
  double t10;
  double t11;
  double t12;
  double t14;
  double t2;
  double t20;
  double t23;
  double t4;
  double t5;
  double t6;
  double t8;
  double t9;
  {
    t1 = q[0];
    t2 = cos(t1);
    t4 = q[1];
    t5 = sin(t4);
    t6 = sin(t1);
    t8 = qdot[0];
    t9 = t8 * t8;
    t10 = 0.981E1 * t6 - t9;
    t11 = t5 * t10;
    t12 = cos(t4);
    t14 = 0.981E1 * t2 * t12;
    t20 = pow(t8 + qdot[1], 2.0);
    t23 = -t11 + t14;
    N[0] = 0.981E1 * t2 - t11 + t14 + t5 * (t12 * t10 + 0.981E1 * t2 * t5 - t20) + t12 * t23;
    N[1] = t23;
    return;
  }
}

