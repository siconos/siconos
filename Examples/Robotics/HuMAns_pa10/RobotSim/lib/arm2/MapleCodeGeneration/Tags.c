#include <math.h>
void Tags(T, q)
double T[24];
double q[7];
{
  double t1;
  double t10;
  double t11;
  double t14;
  double t16;
  double t18;
  double t19;
  double t2;
  double t20;
  double t21;
  double t24;
  double t25;
  double t26;
  double t27;
  double t29;
  double t3;
  double t30;
  double t33;
  double t35;
  double t37;
  double t39;
  double t4;
  double t43;
  double t49;
  double t51;
  double t58;
  double t6;
  double t62;
  double t65;
  double t7;
  double t72;
  double t73;
  double t8;
  {
    T[0] = 0.0;
    T[1] = 0.0;
    t1 = q[0];
    t2 = cos(t1);
    t3 = q[1];
    t4 = sin(t3);
    T[2] = 0.45 * t2 * t4;
    T[3] = T[2];
    t6 = cos(t3);
    t7 = q[2];
    t8 = cos(t7);
    t10 = q[3];
    t11 = sin(t10);
    t14 = cos(t10);
    t16 = -0.45 - 0.48 * t14;
    t18 = 0.48 * t6 * t8 * t11 - t4 * t16;
    t19 = t2 * t18;
    t20 = sin(t1);
    t21 = sin(t7);
    t24 = 0.48 * t20 * t21 * t11;
    t25 = &*(t, tag_5);
    T[4] = t19 - t24 + t25;
    T[5] = t19 - t24;
    T[6] = T[5];
    t26 = q[4];
    t27 = cos(t26);
    t29 = q[5];
    t30 = sin(t29);
    t33 = cos(t29);
    t35 = -0.276576E1 - 0.11004 * t33;
    t37 = 0.11004 * t14 * t27 * t30 - t11 * t35;
    t39 = sin(t26);
    t43 = t8 * t37 - 0.11004 * t21 * t39 * t30;
    t49 = -0.67391E1 + 0.11004 * t11 * t27 * t30 + t14 * t35;
    t51 = t6 * t43 - t4 * t49;
    t58 = t21 * t37 + 0.11004 * t8 * t39 * t30;
    T[7] = 0.3684598379E-1 * t2 * t51 - 0.3684598379E-1 * t20 * t58;
    T[8] = 0.0;
    T[9] = 0.0;
    T[10] = 0.45 * t20 * t4;
    T[11] = T[10];
    t62 = t20 * t18;
    t65 = 0.48 * t2 * t21 * t11;
    T[12] = t62 + t65 + t25;
    T[13] = t62 + t65;
    T[14] = T[13];
    T[15] = 0.3684598379E-1 * t20 * t51 + 0.3684598379E-1 * t2 * t58;
    T[16] = 0.0;
    T[17] = 0.0;
    T[18] = 0.45 * t6;
    T[19] = T[18];
    t72 = 0.48 * t4 * t8 * t11;
    t73 = t6 * t16;
    T[20] = -t72 - t73 + t25;
    T[21] = -t72 - t73;
    T[22] = T[21];
    T[23] = -0.3397199705E-2 - 0.3684598379E-1 * t4 * t43 - 0.3684598379E-1 * t6 * t49;
    return;
  }
}

#include <math.h>
void Tags(T, q)
double T[12];
double q[2];
{
  double t1;
  double t13;
  double t2;
  double t3;
  double t4;
  double t6;
  double t7;
  double t8;
  double t9;
  {
    T[0] = 0.0;
    t1 = q[0];
    T[1] = cos(t1);
    t2 = q[1];
    t3 = cos(t2);
    t4 = t3 + 1.0;
    t6 = sin(t1);
    t7 = sin(t2);
    t8 = t6 * t7;
    T[2] = T[1] * t4 - t8;
    t9 = 2.0 + t3;
    T[3] = T[1] * t9 / 2.0 - t8 / 2.0;
    T[4] = 0.0;
    T[5] = t6;
    t13 = T[1] * t7;
    T[6] = T[5] * t4 + t13;
    T[7] = T[5] * t9 / 2.0 + t13 / 2.0;
    T[8] = 0.0;
    T[9] = 0.0;
    T[10] = 0.0;
    T[11] = 0.0;
    return;
  }
}

#include <math.h>
void Tags(T, q)
double T[12];
double q[2];
{
  double t1;
  double t13;
  double t2;
  double t3;
  double t4;
  double t6;
  double t7;
  double t8;
  double t9;
  {
    T[0] = 0.0;
    t1 = q[0];
    T[1] = cos(t1);
    t2 = q[1];
    t3 = cos(t2);
    t4 = t3 + 1.0;
    t6 = sin(t1);
    t7 = sin(t2);
    t8 = t6 * t7;
    T[2] = T[1] * t4 - t8;
    t9 = 2.0 + t3;
    T[3] = T[1] * t9 / 2.0 - t8 / 2.0;
    T[4] = 0.0;
    T[5] = t6;
    t13 = T[1] * t7;
    T[6] = T[5] * t4 + t13;
    T[7] = T[5] * t9 / 2.0 + t13 / 2.0;
    T[8] = 0.0;
    T[9] = 0.0;
    T[10] = 0.0;
    T[11] = 0.0;
    return;
  }
}

