#include <math.h>
#include "RobotModel.h"

void Friction(double F[NDOF], const double q[NDOF], const double qdot[NDOF])
{
  // Friction  (viscosity model)
  F[0] = FRIC1 * qdot[0];
  F[1] = FRIC2 * qdot[1];
}

