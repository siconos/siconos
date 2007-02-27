#include <math.h>
#include "robotmodel.h"

void modele_coriolis(const double q[N_DOF], const double qdot[N_DOF], double N[N_DOF*N_DOF])
{
  int il, ic;

  for (ic = 0; ic < N_DOF; ic++)
    for (il = 0; il < N_DOF; il++)
      N[N_DOF * ic + il] = 0.0;
}


