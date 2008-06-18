// Copyright (C) INRIA 1999-2008
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 2 as published
// by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
// Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
//%
// @file ActuationModel/NoDynamics/TaskFunctionControl/Trajectory.scilab
// @author Florence Billet
//
// Affiliation(s): INRIA, team BIPOP
//
// Email(s): Florence.Billet@inria.fr
//
// @brief Compute the position, velocity and acceleration desired at a given time t
//


// Description:
//


#include "stdio.h"
#include "math.h"
#include "string.h"

#include "ActuationModelSomeDefinitions.hpp"


void trajectory(double * t, double * position, double * velocity, double * acceleration, int * contacts)
{
  int i;
  double a, r, ry, rz;
  char trajectoryName[20] = "";

  a = M_PI / 3.0;
  r = 0.2;

  for (i = 0; i < 6; i++)
  {
    position[i] = 0;
    velocity[i] = 0;
    acceleration[i] = 0;
  }

  position[0] = (0.45 + 0.07) * cos(M_PI / 2 - M_PI / 3) + 0.38 * cos(M_PI / 2 - M_PI / 3 + M_PI / 6);
  position[1] = r * cos(a * (*t)) + 0.42 + (0.45 + 0.07) * sin(M_PI / 2 - M_PI / 3) + 0.38 * sin(M_PI / 2 - M_PI / 3 + M_PI / 6) - r;
  position[2] = r * sin(a * (*t));
  position[3] = 0.1 * cos(M_PI / 2 - M_PI / 3);
  position[4] = 0.15 * cos(M_PI / 2 - M_PI / 3);

  velocity[1] = -r * a * sin(a * (*t));
  velocity[2] = r * a * cos(a * (*t));

  acceleration[1] = -r * a * a * cos(a * (*t));
  acceleration[2] = -r * a * a * sin(a * (*t));

  *contacts = 0;
}
