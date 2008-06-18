#include <math.h>
void Contact(CC, q)
double CC[6];
double q[7];
{
  double t1;
  double t12;
  double t13;
  double t15;
  double t17;
  double t19;
  double t2;
  double t20;
  double t3;
  double t4;
  double t5;
  double t6;
  double t8;
  double t9;
  {
    CC[0] = 0.0;
    t1 = q[0];
    t2 = cos(t1);
    t3 = q[1];
    t4 = cos(t3);
    t5 = q[2];
    t6 = cos(t5);
    t8 = q[3];
    t9 = sin(t8);
    t12 = sin(t3);
    t13 = cos(t8);
    t15 = -0.45 - 0.48 * t13;
    t17 = 0.48 * t4 * t6 * t9 - t12 * t15;
    t19 = sin(t1);
    t20 = sin(t5);
    CC[1] = t2 * t17 - 0.48 * t19 * t20 * t9;
    CC[2] = 0.0;
    CC[3] = t19 * t17 + 0.48 * t2 * t20 * t9;
    CC[4] = 0.0;
    CC[5] = -0.48 * t12 * t6 * t9 - t4 * t15;
    return;
  }
}

#include <math.h>
void Contact(CC, q)
double CC[6];
double q[2];
{
  double t1;
  double t2;
  double t3;
  double t4;
  double t5;
  double t7;
  double t8;
  {
    CC[0] = 0.0;
    t1 = q[0];
    t2 = cos(t1);
    t3 = q[1];
    t4 = cos(t3);
    t5 = t4 + 1.0;
    t7 = sin(t1);
    t8 = sin(t3);
    CC[1] = t2 * t5 - t7 * t8;
    CC[2] = 0.0;
    CC[3] = t7 * t5 + t2 * t8;
    CC[4] = 0.0;
    CC[5] = 0.0;
    return;
  }
}

