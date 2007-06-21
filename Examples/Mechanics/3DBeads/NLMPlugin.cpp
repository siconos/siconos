/* Siconos-sample version 2.0.0, Copyright INRIA 2005-2006.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
using namespace std;

const double R = 0.1; // beads radius
const double m = 1; // beads mass
const double g = 9.81; // gravity

extern "C" void gravity(double time, unsigned int sizeOfq, double *gravity, unsigned int sizeZ, double* z)
{
  for (unsigned int i = 0; i < sizeOfq; i++)
    gravity[i] = 0.0;

  gravity[2] = -2 * m * g;
}

// // Momentum inertia
const double I1 = 16. / 5 * m * R * R;
const double I3 = 6. / 5 * m * R * R;

extern "C" void mass(unsigned int sizeOfq, const double *q, double *mass, double* param)
{
  mass[0]  = mass[7]  = mass[14] = 2 * m;
  mass[21] = I1;
  mass[28] = I1 * sin(q[3]) * sin(q[3]) + I3 * cos(q[3]) * cos(q[3]);
  mass[29] = mass[34] = I3 * cos(q[3]);
  mass[35] = I3;
}

extern "C" void NNL(unsigned int sizeOfq, const double *q, const double *velocity, double *NNL, double* param)
{
  NNL[3] = -((I1 - I3) * velocity[4] * cos(q[3]) - I3 * velocity[5]) * velocity[4] * sin(q[3]);
  NNL[4] = -(2 * (I3 - I1) * velocity[4] * cos(q[3]) + I3 * velocity[5]) * velocity[3] * sin(q[3]);
  NNL[5] = -I3 * velocity[3] * velocity[4] * sin(q[3]);
}

extern "C" void jacobianQNNL(unsigned int sizeOfq, const double *q, const double *velocity, double *jacob, double* param)
{
  jacob[21] = (I1 - I3) * velocity[4] * velocity[4] * sin(q[3]) * sin(q[3]) - ((I1 - I3) * velocity[4] * cos(q[3]) - I3 * velocity[5]) * velocity[4] * cos(q[3]);
  jacob[22] = 2 * (I3 - I1) * velocity[3] * velocity[4] * sin(q[3]) * sin(q[3]) - (2 * (I3 - I1) * velocity[4] * cos(q[3]) + I3 * velocity[5]) * velocity[3] * cos(q[3]);
  jacob[23] = -I3 * velocity[3] * velocity[4] * cos(q[3]);
}

extern "C" void jacobianVNNL(unsigned int sizeOfq, const double *q, const  double *velocity, double *jacob, double* param)
{
  jacob[22] = -2 * (I1 - I3) * sin(q[3]) * cos(q[3]) * velocity[4] + I3 * velocity[5] * sin(q[3]);
  jacob[23] =  I3 * velocity[4] * sin(q[3]);
  jacob[27] = -(2 * (I3 - I1) * velocity[4] * cos(q[3]) + I3 * velocity[5]) * sin(q[3]);
  jacob[28] = -2 * (I3 - I1) * velocity[3] * sin(q[3]) * cos(q[3]);
  jacob[29] = jacob[34] = -I3 * velocity[3] * sin(q[3]);
  jacob[33] = -I3 * velocity[4] * sin(q[3]);
}

extern "C" void h1(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* y, unsigned int sizeZ, double* param)
{
  double epsilon1 = -1;
  y[0] = q[2] + epsilon1 * R * cos(q[3]) - R;
  //cout << "h1=" << y[0] << endl;
}

extern "C" void G1(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* G, unsigned int sizeZ, double* param)
{
  double epsilon1 = -1;
  G[1] =  G[5] = G[6] = 1.;
  G[9] =  -epsilon1 * R * sin(q[3]);
  G[11] = -epsilon1 * R * (cos(q[3]) - 1);
  G[13] = epsilon1 * R * (cos(q[3]) - 1);
  G[16] = epsilon1 * R * sin(q[3]) * cos(q[4]);
  G[17] = epsilon1 * R * sin(q[3]) * sin(q[4]);
}


extern "C" void h2(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* y, unsigned int sizeZ, double* param)
{
  double epsilon2 = 1;
  y[0] = q[2] + epsilon2 * R * cos(q[3]) - R;
  //cout << "h2=" <<  y[0] << endl;
}

extern "C" void G2(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* G, unsigned int sizeZ, double* param)
{
  double epsilon2 = 1;
  G[1] =  G[5] = G[6] = 1.;
  G[9] =  -epsilon2 * R * sin(q[3]);
  G[11] = -epsilon2 * R * (cos(q[3]) - 1);
  G[13] = epsilon2 * R * (cos(q[3]) - 1);
  G[16] = epsilon2 * R * sin(q[3]) * cos(q[4]);
  G[17] = epsilon2 * R * sin(q[3]) * sin(q[4]);
}


extern "C" void h22(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* y, unsigned int sizeZ, double* param)
{
  double epsilon1 = 1;
  double epsilon2 = 1;

  double d = sqrt(
               ((q[0] + epsilon1 * R * sin(q[3]) * sin(q[4])) - (q[6] + epsilon2 * R * sin(q[9]) * sin(q[10]))) * ((q[0] + epsilon1 * R * sin(q[3]) * sin(q[4])) - (q[6] + epsilon2 * R * sin(q[9]) * sin(q[10])))
               + ((q[1] - epsilon1 * R * sin(q[3]) * cos(q[4])) - (q[7] - epsilon2 * R * sin(q[9]) * cos(q[10]))) * ((q[1] - epsilon1 * R * sin(q[3]) * cos(q[4])) - (q[7] - epsilon2 * R * sin(q[9]) * cos(q[10])))
               + ((q[2] + epsilon1 * R * cos(q[3])) - (q[8] + epsilon2 * R * cos(q[9]))) * ((q[2] + epsilon1 * R * cos(q[3])) - (q[8] + epsilon2 * R * cos(q[9])))
             );
  y[0] =  d - 2 * R;
  cout << "h22 = " << y[0] << endl;
  cout << "d22 = " << d << endl;
}

/* For contact between couples of beads*/ /* A FAIRE */

extern "C" void G22(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* G, unsigned int sizeZ, double* param)
{
  double epsilon1 = 1;
  double epsilon2 = 1;

  double d = sqrt(
               ((q[0] + epsilon1 * R * sin(q[3]) * sin(q[4])) - (q[6] + epsilon2 * R * sin(q[9]) * sin(q[10]))) * ((q[0] + epsilon1 * R * sin(q[3]) * sin(q[4])) - (q[6] + epsilon2 * R * sin(q[9]) * sin(q[10])))
               + ((q[1] - epsilon1 * R * sin(q[3]) * cos(q[4])) - (q[7] - epsilon2 * R * sin(q[9]) * cos(q[10]))) * ((q[1] - epsilon1 * R * sin(q[3]) * cos(q[4])) - (q[7] - epsilon2 * R * sin(q[9]) * cos(q[10])))
               + ((q[2] + epsilon1 * R * cos(q[3])) - (q[8] + epsilon2 * R * cos(q[9]))) * ((q[2] + epsilon1 * R * cos(q[3])) - (q[8] + epsilon2 * R * cos(q[9])))
             );

  G[0] = ((q[0] + epsilon1 * R * sin(q[3]) * sin(q[4])) - (q[6] + epsilon2 * R * sin(q[9]) * sin(q[10]))) / d;
  G[3] = ((q[1] - epsilon1 * R * sin(q[3]) * cos(q[4])) - (q[7] - epsilon2 * R * sin(q[9]) * cos(q[10]))) / d;
  G[6] = ((q[2] + epsilon1 * R * cos(q[3])) - (q[8] + epsilon2 * R * cos(q[9]))) / d;

  G[9] = epsilon1 * R * (((q[0] + epsilon1 * R * sin(q[3]) * sin(q[4])) - (q[6] + epsilon2 * R * sin(q[9]) * sin(q[10]))) * cos(q[3]) * sin(q[4])
                         - ((q[1] - epsilon1 * R * sin(q[3]) * cos(q[4])) - (q[7] - epsilon2 * R * sin(q[9]) * cos(q[10]))) * cos(q[3]) * cos(q[4])
                         - ((q[2] + epsilon1 * R * cos(q[3])) - (q[8] + epsilon2 * R * cos(q[9]))) * sin(q[3])) / d;

  G[12] = epsilon1 * R * (((q[0] + epsilon1 * R * sin(q[3]) * sin(q[4])) - (q[6] + epsilon2 * R * sin(q[9]) * sin(q[10]))) * sin(q[3]) * cos(q[4])
                          + ((q[1] - epsilon1 * R * sin(q[3]) * cos(q[4])) - (q[7] - epsilon2 * R * sin(q[9]) * cos(q[10]))) * sin(q[3]) * sin(q[4])) / d;

  G[18] = -G[0];
  G[21] = -G[1];
  G[24] = -G[2];

  G[27] = -epsilon2 * R * (((q[0] + epsilon1 * R * sin(q[3]) * sin(q[4])) - (q[6] + epsilon2 * R * sin(q[9]) * sin(q[10]))) * cos(q[9]) * sin(q[10])
                           - ((q[1] - epsilon1 * R * sin(q[3]) * cos(q[4])) - (q[7] - epsilon2 * R * sin(q[9]) * cos(q[10]))) * cos(q[9]) * cos(q[10])
                           - ((q[2] + epsilon1 * R * cos(q[3])) - (q[8] + epsilon2 * R * cos(q[9]))) * sin(q[9])) / d;

  G[30] = -epsilon2 * R * (((q[0] + epsilon1 * R * sin(q[3]) * sin(q[4])) - (q[6] + epsilon2 * R * sin(q[9]) * sin(q[10]))) * sin(q[9]) * cos(q[10])
                           + ((q[1] - epsilon1 * R * sin(q[3]) * cos(q[4])) - (q[7] - epsilon2 * R * sin(q[9]) * cos(q[10]))) * sin(q[9]) * sin(q[10])) / d;

}

extern "C" void h21(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* y, unsigned int sizeZ, double* param)
{
  double epsilon1 = 1;
  double epsilon2 = -1;

  double d = sqrt(
               ((q[0] + epsilon1 * R * sin(q[3]) * sin(q[4])) - (q[6] + epsilon2 * R * sin(q[9]) * sin(q[10]))) * ((q[0] + epsilon1 * R * sin(q[3]) * sin(q[4])) - (q[6] + epsilon2 * R * sin(q[9]) * sin(q[10])))
               + ((q[1] - epsilon1 * R * sin(q[3]) * cos(q[4])) - (q[7] - epsilon2 * R * sin(q[9]) * cos(q[10]))) * ((q[1] - epsilon1 * R * sin(q[3]) * cos(q[4])) - (q[7] - epsilon2 * R * sin(q[9]) * cos(q[10])))
               + ((q[2] + epsilon1 * R * cos(q[3])) - (q[8] + epsilon2 * R * cos(q[9]))) * ((q[2] + epsilon1 * R * cos(q[3])) - (q[8] + epsilon2 * R * cos(q[9])))
             );
  y[0] =  d - 2 * R;
  cout << "h21 = " << y[0] << endl;
  cout << "d21 = " << d << endl;
}

extern "C" void G21(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* G, unsigned int sizeZ, double* param)
{
  double epsilon1 = 1;
  double epsilon2 = -1;

  double d = sqrt(
               ((q[0] + epsilon1 * R * sin(q[3]) * sin(q[4])) - (q[6] + epsilon2 * R * sin(q[9]) * sin(q[10]))) * ((q[0] + epsilon1 * R * sin(q[3]) * sin(q[4])) - (q[6] + epsilon2 * R * sin(q[9]) * sin(q[10])))
               + ((q[1] - epsilon1 * R * sin(q[3]) * cos(q[4])) - (q[7] - epsilon2 * R * sin(q[9]) * cos(q[10]))) * ((q[1] - epsilon1 * R * sin(q[3]) * cos(q[4])) - (q[7] - epsilon2 * R * sin(q[9]) * cos(q[10])))
               + ((q[2] + epsilon1 * R * cos(q[3])) - (q[8] + epsilon2 * R * cos(q[9]))) * ((q[2] + epsilon1 * R * cos(q[3])) - (q[8] + epsilon2 * R * cos(q[9])))
             );

  G[0] = ((q[0] + epsilon1 * R * sin(q[3]) * sin(q[4])) - (q[6] + epsilon2 * R * sin(q[9]) * sin(q[10]))) / d;
  G[3] = ((q[1] - epsilon1 * R * sin(q[3]) * cos(q[4])) - (q[7] - epsilon2 * R * sin(q[9]) * cos(q[10]))) / d;
  G[6] = ((q[2] + epsilon1 * R * cos(q[3])) - (q[8] + epsilon2 * R * cos(q[9]))) / d;

  G[9] = epsilon1 * R * (((q[0] + epsilon1 * R * sin(q[3]) * sin(q[4])) - (q[6] + epsilon2 * R * sin(q[9]) * sin(q[10]))) * cos(q[3]) * sin(q[4])
                         - ((q[1] - epsilon1 * R * sin(q[3]) * cos(q[4])) - (q[7] - epsilon2 * R * sin(q[9]) * cos(q[10]))) * cos(q[3]) * cos(q[4])
                         - ((q[2] + epsilon1 * R * cos(q[3])) - (q[8] + epsilon2 * R * cos(q[9]))) * sin(q[3])) / d;

  G[12] = epsilon1 * R * (((q[0] + epsilon1 * R * sin(q[3]) * sin(q[4])) - (q[6] + epsilon2 * R * sin(q[9]) * sin(q[10]))) * sin(q[3]) * cos(q[4])
                          + ((q[1] - epsilon1 * R * sin(q[3]) * cos(q[4])) - (q[7] - epsilon2 * R * sin(q[9]) * cos(q[10]))) * sin(q[3]) * sin(q[4])) / d;

  G[18] = -G[0];
  G[21] = -G[1];
  G[24] = -G[2];

  G[27] = -epsilon2 * R * (((q[0] + epsilon1 * R * sin(q[3]) * sin(q[4])) - (q[6] + epsilon2 * R * sin(q[9]) * sin(q[10]))) * cos(q[9]) * sin(q[10])
                           - ((q[1] - epsilon1 * R * sin(q[3]) * cos(q[4])) - (q[7] - epsilon2 * R * sin(q[9]) * cos(q[10]))) * cos(q[9]) * cos(q[10])
                           - ((q[2] + epsilon1 * R * cos(q[3])) - (q[8] + epsilon2 * R * cos(q[9]))) * sin(q[9])) / d;

  G[30] = -epsilon2 * R * (((q[0] + epsilon1 * R * sin(q[3]) * sin(q[4])) - (q[6] + epsilon2 * R * sin(q[9]) * sin(q[10]))) * sin(q[9]) * cos(q[10])
                           + ((q[1] - epsilon1 * R * sin(q[3]) * cos(q[4])) - (q[7] - epsilon2 * R * sin(q[9]) * cos(q[10]))) * sin(q[9]) * sin(q[10])) / d;


}

extern "C" void h12(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* y, unsigned int sizeZ, double* param)
{
  double epsilon1 = -1;
  double epsilon2 = 1;

  double d = sqrt(
               ((q[0] + epsilon1 * R * sin(q[3]) * sin(q[4])) - (q[6] + epsilon2 * R * sin(q[9]) * sin(q[10]))) * ((q[0] + epsilon1 * R * sin(q[3]) * sin(q[4])) - (q[6] + epsilon2 * R * sin(q[9]) * sin(q[10])))
               + ((q[1] - epsilon1 * R * sin(q[3]) * cos(q[4])) - (q[7] - epsilon2 * R * sin(q[9]) * cos(q[10]))) * ((q[1] - epsilon1 * R * sin(q[3]) * cos(q[4])) - (q[7] - epsilon2 * R * sin(q[9]) * cos(q[10])))
               + ((q[2] + epsilon1 * R * cos(q[3])) - (q[8] + epsilon2 * R * cos(q[9]))) * ((q[2] + epsilon1 * R * cos(q[3])) - (q[8] + epsilon2 * R * cos(q[9])))
             );

  y[0] =  d - 2 * R;
  cout << "h12 = " << y[0] << endl;
  cout << "d12 = " << d << endl;
}

/* For contact between couples of beads*/ /* A FAIRE */

extern "C" void G12(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* G, unsigned int sizeZ, double* param)
{
  double epsilon1 = -1;
  double epsilon2 = 1;

  double d = sqrt(
               ((q[0] + epsilon1 * R * sin(q[3]) * sin(q[4])) - (q[6] + epsilon2 * R * sin(q[9]) * sin(q[10]))) * ((q[0] + epsilon1 * R * sin(q[3]) * sin(q[4])) - (q[6] + epsilon2 * R * sin(q[9]) * sin(q[10])))
               + ((q[1] - epsilon1 * R * sin(q[3]) * cos(q[4])) - (q[7] - epsilon2 * R * sin(q[9]) * cos(q[10]))) * ((q[1] - epsilon1 * R * sin(q[3]) * cos(q[4])) - (q[7] - epsilon2 * R * sin(q[9]) * cos(q[10])))
               + ((q[2] + epsilon1 * R * cos(q[3])) - (q[8] + epsilon2 * R * cos(q[9]))) * ((q[2] + epsilon1 * R * cos(q[3])) - (q[8] + epsilon2 * R * cos(q[9])))
             );

  G[0] = ((q[0] + epsilon1 * R * sin(q[3]) * sin(q[4])) - (q[6] + epsilon2 * R * sin(q[9]) * sin(q[10]))) / d;
  G[3] = ((q[1] - epsilon1 * R * sin(q[3]) * cos(q[4])) - (q[7] - epsilon2 * R * sin(q[9]) * cos(q[10]))) / d;
  G[6] = ((q[2] + epsilon1 * R * cos(q[3])) - (q[8] + epsilon2 * R * cos(q[9]))) / d;

  G[9] = epsilon1 * R * (((q[0] + epsilon1 * R * sin(q[3]) * sin(q[4])) - (q[6] + epsilon2 * R * sin(q[9]) * sin(q[10]))) * cos(q[3]) * sin(q[4])
                         - ((q[1] - epsilon1 * R * sin(q[3]) * cos(q[4])) - (q[7] - epsilon2 * R * sin(q[9]) * cos(q[10]))) * cos(q[3]) * cos(q[4])
                         - ((q[2] + epsilon1 * R * cos(q[3])) - (q[8] + epsilon2 * R * cos(q[9]))) * sin(q[3])) / d;

  G[12] = epsilon1 * R * (((q[0] + epsilon1 * R * sin(q[3]) * sin(q[4])) - (q[6] + epsilon2 * R * sin(q[9]) * sin(q[10]))) * sin(q[3]) * cos(q[4])
                          + ((q[1] - epsilon1 * R * sin(q[3]) * cos(q[4])) - (q[7] - epsilon2 * R * sin(q[9]) * cos(q[10]))) * sin(q[3]) * sin(q[4])) / d;

  G[18] = -G[0];
  G[21] = -G[1];
  G[24] = -G[2];

  G[27] = -epsilon2 * R * (((q[0] + epsilon1 * R * sin(q[3]) * sin(q[4])) - (q[6] + epsilon2 * R * sin(q[9]) * sin(q[10]))) * cos(q[9]) * sin(q[10])
                           - ((q[1] - epsilon1 * R * sin(q[3]) * cos(q[4])) - (q[7] - epsilon2 * R * sin(q[9]) * cos(q[10]))) * cos(q[9]) * cos(q[10])
                           - ((q[2] + epsilon1 * R * cos(q[3])) - (q[8] + epsilon2 * R * cos(q[9]))) * sin(q[9])) / d;

  G[30] = -epsilon2 * R * (((q[0] + epsilon1 * R * sin(q[3]) * sin(q[4])) - (q[6] + epsilon2 * R * sin(q[9]) * sin(q[10]))) * sin(q[9]) * cos(q[10])
                           + ((q[1] - epsilon1 * R * sin(q[3]) * cos(q[4])) - (q[7] - epsilon2 * R * sin(q[9]) * cos(q[10]))) * sin(q[9]) * sin(q[10])) / d;

}

extern "C" void h11(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* y, unsigned int sizeZ, double* param)
{
  double epsilon1 = -1;
  double epsilon2 = -1;

  double d = sqrt(
               ((q[0] + epsilon1 * R * sin(q[3]) * sin(q[4])) - (q[6] + epsilon2 * R * sin(q[9]) * sin(q[10]))) * ((q[0] + epsilon1 * R * sin(q[3]) * sin(q[4])) - (q[6] + epsilon2 * R * sin(q[9]) * sin(q[10])))
               + ((q[1] - epsilon1 * R * sin(q[3]) * cos(q[4])) - (q[7] - epsilon2 * R * sin(q[9]) * cos(q[10]))) * ((q[1] - epsilon1 * R * sin(q[3]) * cos(q[4])) - (q[7] - epsilon2 * R * sin(q[9]) * cos(q[10])))
               + ((q[2] + epsilon1 * R * cos(q[3])) - (q[8] + epsilon2 * R * cos(q[9]))) * ((q[2] + epsilon1 * R * cos(q[3])) - (q[8] + epsilon2 * R * cos(q[9])))
             );

  y[0] =  d - 2 * R;
  cout << "h11 = " << y[0] << endl;
  cout << "d11 = " << d << endl;
  cout << " " << endl;
  cout << " " << endl;
}

extern "C" void G11(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* G, unsigned int sizeZ, double* param)
{
  double epsilon1 = -1;
  double epsilon2 = -1;

  double d = sqrt(
               ((q[0] + epsilon1 * R * sin(q[3]) * sin(q[4])) - (q[6] + epsilon2 * R * sin(q[9]) * sin(q[10]))) * ((q[0] + epsilon1 * R * sin(q[3]) * sin(q[4])) - (q[6] + epsilon2 * R * sin(q[9]) * sin(q[10])))
               + ((q[1] - epsilon1 * R * sin(q[3]) * cos(q[4])) - (q[7] - epsilon2 * R * sin(q[9]) * cos(q[10]))) * ((q[1] - epsilon1 * R * sin(q[3]) * cos(q[4])) - (q[7] - epsilon2 * R * sin(q[9]) * cos(q[10])))
               + ((q[2] + epsilon1 * R * cos(q[3])) - (q[8] + epsilon2 * R * cos(q[9]))) * ((q[2] + epsilon1 * R * cos(q[3])) - (q[8] + epsilon2 * R * cos(q[9])))
             );
  G[0] = ((q[0] + epsilon1 * R * sin(q[3]) * sin(q[4])) - (q[6] + epsilon2 * R * sin(q[9]) * sin(q[10]))) / d;
  G[3] = ((q[1] - epsilon1 * R * sin(q[3]) * cos(q[4])) - (q[7] - epsilon2 * R * sin(q[9]) * cos(q[10]))) / d;
  G[6] = ((q[2] + epsilon1 * R * cos(q[3])) - (q[8] + epsilon2 * R * cos(q[9]))) / d;

  G[9] = epsilon1 * R * (((q[0] + epsilon1 * R * sin(q[3]) * sin(q[4])) - (q[6] + epsilon2 * R * sin(q[9]) * sin(q[10]))) * cos(q[3]) * sin(q[4])
                         - ((q[1] - epsilon1 * R * sin(q[3]) * cos(q[4])) - (q[7] - epsilon2 * R * sin(q[9]) * cos(q[10]))) * cos(q[3]) * cos(q[4])
                         - ((q[2] + epsilon1 * R * cos(q[3])) - (q[8] + epsilon2 * R * cos(q[9]))) * sin(q[3])) / d;

  G[12] = epsilon1 * R * (((q[0] + epsilon1 * R * sin(q[3]) * sin(q[4])) - (q[6] + epsilon2 * R * sin(q[9]) * sin(q[10]))) * sin(q[3]) * cos(q[4])
                          + ((q[1] - epsilon1 * R * sin(q[3]) * cos(q[4])) - (q[7] - epsilon2 * R * sin(q[9]) * cos(q[10]))) * sin(q[3]) * sin(q[4])) / d;

  G[18] = -G[0];
  G[21] = -G[1];
  G[24] = -G[2];

  G[27] = -epsilon2 * R * (((q[0] + epsilon1 * R * sin(q[3]) * sin(q[4])) - (q[6] + epsilon2 * R * sin(q[9]) * sin(q[10]))) * cos(q[9]) * sin(q[10])
                           - ((q[1] - epsilon1 * R * sin(q[3]) * cos(q[4])) - (q[7] - epsilon2 * R * sin(q[9]) * cos(q[10]))) * cos(q[9]) * cos(q[10])
                           - ((q[2] + epsilon1 * R * cos(q[3])) - (q[8] + epsilon2 * R * cos(q[9]))) * sin(q[9])) / d;

  G[30] = -epsilon2 * R * (((q[0] + epsilon1 * R * sin(q[3]) * sin(q[4])) - (q[6] + epsilon2 * R * sin(q[9]) * sin(q[10]))) * sin(q[9]) * cos(q[10])
                           + ((q[1] - epsilon1 * R * sin(q[3]) * cos(q[4])) - (q[7] - epsilon2 * R * sin(q[9]) * cos(q[10]))) * sin(q[9]) * sin(q[10])) / d;

}
