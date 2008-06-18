#include <math.h>

#include "robotmodel.h"

void modele_inertie(const double q[NDDL], double M[NDDL*NDDL])
{
  double c2;
  // Compute cos and sin
  c2 = cos(q[1]);


  M[0] = L2 * L2 * M2 + 2 * L1 * L2 * M2 * c2 + L1 * L1 * (M1 + M2);
  M[1] = L2 * L2 * M2 + L1 * L2 * M2 * c2;
  M[2] = M[1];
  M[3] = L2 * L2 * M2;
}
