#include <math.h>
void Contact(double * CC, double * q)
{
  double t1;
  double t10;
  double t11;
  double t14;
  double t16;
  double t18;
  double t2;
  double t20;
  double t21;
  double t3;
  double t4;
  double t6;
  double t7;
  double t8;
  {
    t1 = q[0];
    t2 = cos(t1);
    t3 = q[1];
    t4 = sin(t3);
    CC[0] = 0.45 * t2 * t4;
    t6 = cos(t3);
    t7 = q[2];
    t8 = cos(t7);
    t10 = q[3];
    t11 = sin(t10);
    t14 = cos(t10);
    t16 = -0.45 - 0.48 * t14;
    t18 = 0.48 * t6 * t8 * t11 - t4 * t16;
    t20 = sin(t1);
    t21 = sin(t7);
    CC[1] = t2 * t18 - 0.48 * t20 * t21 * t11;
    CC[2] = 0.45 * t20 * t4;
    CC[3] = t20 * t18 + 0.48 * t2 * t21 * t11;
    CC[4] = 0.45 * t6;
    CC[5] = -0.48 * t4 * t8 * t11 - t6 * t16;
    return;
  }
}

