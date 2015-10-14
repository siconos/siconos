#include <math.h>
#include "robotmodel.h"

void modele_gravite(const double q[NDDL], double G[NDDL])
{

  double c1, c12;
  // Compute cos and sin
  c1 = cos(q[0]);
  c12 = cos(q[0] + q[1]);

  G[0] = M2 * L2 * GRAV * c12 + (M1 + M2) * L1 * GRAV * c1;
  G[1] = M2 * L2 * GRAV * c12;
}
