#include <math.h>
void ContactJacobian(CJ, q)
double CJ[42];
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
  double t21;
  double t25;
  double t28;
  double t3;
  double t32;
  double t34;
  double t4;
  double t5;
  double t52;
  double t6;
  double t7;
  double t8;
  double t9;
  {
    CJ[0] = 0.0;
    t1 = q[0];
    t2 = sin(t1);
    t3 = q[1];
    t4 = cos(t3);
    t5 = q[2];
    t6 = cos(t5);
    t7 = t4 * t6;
    t8 = q[3];
    t9 = sin(t8);
    t12 = sin(t3);
    t13 = cos(t8);
    t15 = -0.45 - 0.48 * t13;
    t17 = 0.48 * t7 * t9 - t12 * t15;
    t19 = cos(t1);
    t20 = sin(t5);
    t21 = t19 * t20;
    CJ[1] = -t2 * t17 - 0.48 * t21 * t9;
    CJ[2] = 0.0;
    t25 = t2 * t20;
    CJ[3] = t19 * t17 - 0.48 * t25 * t9;
    CJ[4] = 0.0;
    CJ[5] = 0.0;
    CJ[6] = 0.0;
    t28 = t12 * t6;
    t32 = -0.48 * t28 * t9 - t4 * t15;
    CJ[7] = t19 * t32;
    CJ[8] = 0.0;
    CJ[9] = t2 * t32;
    CJ[10] = 0.0;
    CJ[11] = -t17;
    CJ[12] = 0.0;
    t34 = t20 * t9;
    CJ[13] = -0.48 * t19 * t4 * t34 - 0.48 * t2 * t6 * t9;
    CJ[14] = 0.0;
    CJ[15] = -0.48 * t2 * t4 * t34 + 0.48 * t19 * t6 * t9;
    CJ[16] = 0.0;
    CJ[17] = 0.48 * t12 * t20 * t9;
    CJ[18] = 0.0;
    t52 = 0.48 * t7 * t13 - 0.48 * t9 * t12;
    CJ[19] = t19 * t52 - 0.48 * t25 * t13;
    CJ[20] = 0.0;
    CJ[21] = t2 * t52 + 0.48 * t21 * t13;
    CJ[22] = 0.0;
    CJ[23] = -0.48 * t28 * t13 - 0.48 * t4 * t9;
    CJ[24] = 0.0;
    CJ[25] = 0.0;
    CJ[26] = 0.0;
    CJ[27] = 0.0;
    CJ[28] = 0.0;
    CJ[29] = 0.0;
    CJ[30] = 0.0;
    CJ[31] = 0.0;
    CJ[32] = 0.0;
    CJ[33] = 0.0;
    CJ[34] = 0.0;
    CJ[35] = 0.0;
    CJ[36] = 0.0;
    CJ[37] = 0.0;
    CJ[38] = 0.0;
    CJ[39] = 0.0;
    CJ[40] = 0.0;
    CJ[41] = 0.0;
    return;
  }
}

#include <math.h>
void ContactJacobian(CJ, q)
double CJ[12];
double q[2];
{
  double t1;
  double t11;
  double t2;
  double t3;
  double t4;
  double t5;
  double t7;
  double t8;
  double t9;
  {
    CJ[0] = 0.0;
    t1 = q[0];
    t2 = sin(t1);
    t3 = q[1];
    t4 = cos(t3);
    t5 = t4 + 1.0;
    t7 = cos(t1);
    t8 = sin(t3);
    t9 = t7 * t8;
    CJ[1] = -t2 * t5 - t9;
    CJ[2] = 0.0;
    t11 = t2 * t8;
    CJ[3] = t7 * t5 - t11;
    CJ[4] = 0.0;
    CJ[5] = 0.0;
    CJ[6] = 0.0;
    CJ[7] = -t9 - t2 * t4;
    CJ[8] = 0.0;
    CJ[9] = -t11 + t7 * t4;
    CJ[10] = 0.0;
    CJ[11] = 0.0;
    return;
  }
}

