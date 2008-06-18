#include <math.h>
void ContactHessian(CH, q, qdot)
double CH[6];
double q[7];
double qdot[7];
{
  double t1;
  double t10;
  double t11;
  double t112;
  double t12;
  double t126;
  double t15;
  double t16;
  double t18;
  double t2;
  double t20;
  double t22;
  double t23;
  double t24;
  double t25;
  double t26;
  double t29;
  double t3;
  double t30;
  double t31;
  double t35;
  double t39;
  double t4;
  double t44;
  double t47;
  double t48;
  double t5;
  double t52;
  double t56;
  double t59;
  double t6;
  double t60;
  double t64;
  double t67;
  double t68;
  double t7;
  double t71;
  double t73;
  double t79;
  double t8;
  double t81;
  double t82;
  double t85;
  double t87;
  double t9;
  double t91;
  double t94;
  {
    CH[0] = 0.0;
    t1 = q[0];
    t2 = cos(t1);
    t3 = qdot[0];
    t4 = t3 * t3;
    t5 = t2 * t4;
    t6 = q[1];
    t7 = cos(t6);
    t8 = q[2];
    t9 = cos(t8);
    t10 = t7 * t9;
    t11 = q[3];
    t12 = sin(t11);
    t15 = sin(t6);
    t16 = cos(t11);
    t18 = -0.45 - 0.48 * t16;
    t20 = 0.48 * t10 * t12 - t15 * t18;
    t22 = sin(t1);
    t23 = t22 * t3;
    t24 = qdot[1];
    t25 = t15 * t24;
    t26 = t9 * t12;
    t29 = sin(t8);
    t30 = t7 * t29;
    t31 = qdot[2];
    t35 = qdot[3];
    t39 = t7 * t24;
    t44 = -0.48 * t25 * t26 - 0.48 * t30 * t31 * t12 + 0.48 * t10 * t16 * t35 - t39 * t18 - 0.48 * t15 * t12 *
          t35;
    t47 = t24 * t24;
    t48 = t7 * t47;
    t52 = t29 * t31 * t12;
    t56 = t9 * t16 * t35;
    t59 = t31 * t31;
    t60 = t59 * t12;
    t64 = t31 * t16 * t35;
    t67 = t35 * t35;
    t68 = t12 * t67;
    t71 = t15 * t47;
    t73 = t12 * t35;
    t79 = -0.48 * t48 * t26 + 0.96 * t52 * t25 - 0.96 * t25 * t56 - 0.48 * t10 * t60 - 0.96 * t30 * t64
          - 0.48 * t10 * t68 + t71 * t18 - 0.96 * t39 * t73 - 0.48 * t15 * t16 * t67;
    t81 = t22 * t4;
    t82 = t29 * t12;
    t85 = t2 * t3;
    t87 = t9 * t31 * t12;
    t91 = t29 * t16 * t35;
    t94 = t22 * t29;
    CH[1] = -t5 * t20 - 2.0 * t23 * t44 + t2 * t79 + 0.48 * t81 * t82 - 0.96 * t85 * t87 - 0.96 * t85 * t91 +
            0.48 * t94 * t60 - 0.96 * t22 * t9 * t64 + 0.48 * t94 * t68;
    CH[2] = 0.0;
    t112 = t2 * t29;
    CH[3] = -t81 * t20 + 2.0 * t85 * t44 + t22 * t79 - 0.48 * t5 * t82 - 0.96 * t23 * t87 - 0.96 * t23 * t91
            - 0.48 * t112 * t60 + 0.96 * t2 * t9 * t64 - 0.48 * t112 * t68;
    CH[4] = 0.0;
    t126 = t15 * t9;
    CH[5] = 0.48 * t71 * t26 + 0.96 * t39 * t52 - 0.96 * t39 * t56 + 0.48 * t126 * t60 + 0.96 * t15 * t29 *
            t64 + 0.48 * t126 * t68 + t48 * t18 + 0.96 * t25 * t73 - 0.48 * t7 * t16 * t67;
    return;
  }
}

#include <math.h>
void ContactHessian(CH, q, qdot)
double CH[6];
double q[2];
double qdot[2];
{
  double t1;
  double t10;
  double t11;
  double t12;
  double t13;
  double t14;
  double t18;
  double t2;
  double t20;
  double t22;
  double t23;
  double t3;
  double t4;
  double t5;
  double t6;
  double t7;
  double t8;
  {
    CH[0] = 0.0;
    t1 = q[0];
    t2 = cos(t1);
    t3 = qdot[0];
    t4 = t3 * t3;
    t5 = t2 * t4;
    t6 = q[1];
    t7 = cos(t6);
    t8 = t7 + 1.0;
    t10 = sin(t1);
    t11 = t10 * t3;
    t12 = sin(t6);
    t13 = qdot[1];
    t14 = t12 * t13;
    t18 = t13 * t13;
    t20 = t10 * t4;
    t22 = t2 * t3;
    t23 = t7 * t13;
    CH[1] = -t5 * t8 + 2.0 * t11 * t14 - t2 * t7 * t18 + t20 * t12 - 2.0 * t22 * t23 + t10 * t12 * t18;
    CH[2] = 0.0;
    CH[3] = -t20 * t8 - 2.0 * t22 * t14 - t10 * t7 * t18 - t5 * t12 - 2.0 * t11 * t23 - t2 * t12 * t18;
    CH[4] = 0.0;
    CH[5] = 0.0;
    return;
  }
}

