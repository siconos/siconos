#include <math.h>
#include "robotmodel.h"

void modele_coriolis(const double q[NDDL], const double qdot[NDDL], double N[NDDL*NDDL])
{
  double s2;
  // Compute cos and sin
  s2 = sin(q[1]);


  N[0] = -2 * M2 * L1 * L2 * s2 * qdot[1];
  N[1] = -M2 * L1 * L2 * s2 * qdot[1];
  N[2] = M2 * L1 * L2 * qdot[0];
  N[3] = 0.0;
}

