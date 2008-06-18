#include <math.h>
#include "robotmodel.h"

void modele_frottements(const double q[N_DOF], const double qdot[N_DOF], double F[N_DOF])
{
  // Friction  (viscosity model)
  F[0] = FRIC1 * qdot[0];
  F[1] = FRIC2 * qdot[1];
  F[2] = FRIC3 * qdot[2];
  F[3] = FRIC4 * qdot[3];
  F[4] = FRIC5 * qdot[4];
  F[5] = FRIC6 * qdot[5];
  F[6] = FRIC7 * qdot[6];
}

