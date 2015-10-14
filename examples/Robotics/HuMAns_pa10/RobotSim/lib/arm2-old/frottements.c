#include <math.h>
#include "robotmodel.h"

void modele_frottements(const double q[NDDL], const double qdot[NDDL], double F[NDDL])
{
  // Friction  (viscosity model)
  F[0] = FRIC1 * qdot[0];
  F[1] = FRIC2 * qdot[1];
}

