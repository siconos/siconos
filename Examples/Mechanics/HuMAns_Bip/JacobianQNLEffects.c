#include "JacobianQNLEffects.h"

void JacobianQNLEffects(double *jac, const double *q, const double *qdot)
{
  int i;
  double JC[63];

  JacobianQNLEffects1(JC, q, qdot);
  for (i = 0; i < 63; i++)
    jac[i] = JC[i];

  JacobianQNLEffects2(JC, q, qdot);
  for (i = 0; i < 63; i++)
    jac[63 + i] = JC[i];

  JacobianQNLEffects3(JC, q, qdot);
  for (i = 0; i < 63; i++)
    jac[2 * 63 + i] = JC[i];

  JacobianQNLEffects4(JC, q, qdot);
  for (i = 0; i < 63; i++)
    jac[3 * 63 + i] = JC[i];

  JacobianQNLEffects5(JC, q, qdot);
  for (i = 0; i < 63; i++)
    jac[4 * 63 + i] = JC[i];

  JacobianQNLEffects6(JC, q, qdot);
  for (i = 0; i < 63; i++)
    jac[5 * 63 + i] = JC[i];

  JacobianQNLEffects7(JC, q, qdot);
  for (i = 0; i < 63; i++)
    jac[6 * 63 + i] = JC[i];
}
