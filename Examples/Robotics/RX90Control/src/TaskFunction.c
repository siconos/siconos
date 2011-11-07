#include "TaskFunctionDefinition.h"
#include <math.h>
void TaskFunction(T, q)
double T[6];
double q[6];
{
  double t1;
  double t10;
  double t12;
  double t14;
  double t18;
  double t2;
  double t20;
  double t21;
  double t23;
  double t24;
  double t26;
  double t3;
  double t35;
  double t36;
  double t38;
  double t4;
  double t40;
  double t45;
  double t46;
  double t48;
  double t5;
  double t50;
  double t6;
  double t9;
  {
    t1 = q[0];
    t2 = cos(t1);
    t3 = q[1];
    t4 = cos(t3);
    t5 = q[2];
    t6 = sin(t5);
    t9 = sin(t3);
    t10 = cos(t5);
    t12 = 0.45 + 0.45 * t10;
    t14 = 0.45 * t4 * t6 + t9 * t12;
    T[0] = t2 * t14;
    T[1] = 0.42 - 0.45 * t9 * t6 + t4 * t12;
    t18 = sin(t1);
    T[2] = -t18 * t14;
    t20 = q[3];
    t21 = cos(t20);
    t23 = q[4];
    t24 = sin(t23);
    t26 = cos(t23);
    t35 = t4 * (t10 * t21 * t24 + t6 * t26) + t9 * (-t6 * t21 * t24 + t10 * t26);
    t36 = t2 * t35;
    t38 = sin(t20);
    t40 = t18 * t38 * t24;
    T[3] = 0.1 * t36 - 0.1 * t40;
    t45 = q[5];
    t46 = sin(t45);
    t48 = cos(t45);
    t50 = t21 * t26 * t46 + t38 * t48;
    T[4] = 0.15 * t36 - 0.15 * t40 + 0.5E-1 * t2 * (t4 * (t10 * t50 - t6 * t24 * t46) + t9 * (-t6 * t50 - t10
           * t24 * t46)) + 0.5E-1 * t18 * (-t38 * t26 * t46 + t21 * t48);
    T[5] = -0.1 * t18 * t35 - 0.1 * t2 * t38 * t24;
    return;
  }
}

