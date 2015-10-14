#include <math.h>
#include "RobotModel.h"

#define EPS 0.0001
#define SIGN(a) ((a) < -EPS ? (-1): ((a) > EPS ?(1): (a/EPS)))

void Friction(double F[NDOF], const double q[NDOF], const double qdot[NDOF])
{
  // Friction  (viscosity model)
  F[0] = SIGN(qdot[0]) * SFRIC1 + VFRIC1 * qdot[0];
  F[1] = SIGN(qdot[1]) * SFRIC2 + VFRIC2 * qdot[1];
  F[2] = SIGN(qdot[2]) * SFRIC3 + VFRIC3 * qdot[2];
  F[3] = SIGN(qdot[3]) * SFRIC4 + VFRIC4 * qdot[3];
  F[4] = SIGN(qdot[4]) * SFRIC5 + VFRIC5 * qdot[4];
  F[5] = SIGN(qdot[5]) * SFRIC6 + VFRIC6 * qdot[5];
  F[6] = SIGN(qdot[6]) * SFRIC7 + VFRIC7 * qdot[6];
}



