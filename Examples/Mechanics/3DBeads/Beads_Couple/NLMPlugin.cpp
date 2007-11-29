/* Siconos-sample version 2.0.0, Copyright INRIA 2005-2007.
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

const double R      = 0.1; // beads radius
const double m      = 1.; // beads mass
const double g      = 9.81; // gravity
const double Wall   = 1.; // Wall Position
const double Top    = 2.2; // Top Wall Position
const double Ground = 0.; // Ground Wall Position

extern "C" void gravity(double time, unsigned int sizeOfq, double *gravity, unsigned int sizeZ, double* z)
{
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

/* Interactions with walls with respect to z axis
   We have two interactions for each couple concerning each bead and for each couple you have also interactions with top and ground
*/

extern "C" void h1z(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* y, unsigned int sizeZ, double* param)
{
  double epsilon1 = -1;
  y[0] = q[2] + epsilon1 * R * cos(q[3]) - R - Ground;
  //cout << "h1=" << y[0] << endl;
}

extern "C" void G1z(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* G, unsigned int sizeZ, double* param)
{
  double epsilon1 = -1;

  G[6]  = 1.;
  G[9]  = -epsilon1 * R * sin(q[3]);

  G[1]  =  G[5] = 1.;
  G[11] =  epsilon1 * R;
  G[13] = -epsilon1 * R;

}


extern "C" void h2z(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* y, unsigned int sizeZ, double* param)
{
  double epsilon2 = 1;
  y[0] = q[2] + epsilon2 * R * cos(q[3]) - R - Ground;
  //cout << "h2=" <<  y[0] << endl;
}

extern "C" void G2z(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* G, unsigned int sizeZ, double* param)
{
  double epsilon2 = 1;

  G[6]  = 1.;
  G[9]  = -epsilon2 * R * sin(q[3]);

  G[1]  = G[5] = 1.;
  G[11] = epsilon2 * R;
  G[13] = -epsilon2 * R;
}

extern "C" void h1z_(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* y, unsigned int sizeZ, double* param)
{
  double epsilon1 = -1;
  y[0] = Top - (q[2] + epsilon1 * R * cos(q[3])) - R;
  //cout << "h1=" << y[0] << endl;
}

extern "C" void G1z_(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* G, unsigned int sizeZ, double* param)
{
  double epsilon1 = -1;

  G[6]  = -1.;
  G[9]  = epsilon1 * R * sin(q[3]);

  G[1]  = G[5] = 1.;
  G[11] = epsilon1 * R;
  G[13] = -epsilon1 * R;
}

extern "C" void h2z_(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* y, unsigned int sizeZ, double* param)
{
  double epsilon2 = 1;
  y[0] = Top - (q[2] + epsilon2 * R * cos(q[3])) - R;
  //cout << "h1=" << y[0] << endl;
}

extern "C" void G2z_(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* G, unsigned int sizeZ, double* param)
{
  double epsilon2 = 1;

  G[6]  = -1.;
  G[9]  = epsilon2 * R * sin(q[3]);

  G[1]  = G[5] = 1.;
  G[11] = epsilon2 * R;
  G[13] = -epsilon2 * R;
}


/* Interactions with walls with respect to y axis
   We have two interactions for each couple concerning each bead and for each couple you have also interactions with right an left wall
*/

extern "C" void h1y(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* y, unsigned int sizeZ, double* param)
{
  double epsilon1 = -1;
  y[0] = q[1] - epsilon1 * R * sin(q[3]) * cos(q[4]) - R + Wall;
  //cout << "h1=" << y[0] << endl;
}

extern "C" void G1y(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* G, unsigned int sizeZ, double* param)
{
  double epsilon1 = -1;

  G[3]  = 1.;
  G[9]  = -epsilon1 * R * cos(q[3]) * cos(q[4]);
  G[12] =  epsilon1 * R * sin(q[3]) * sin(q[4]);

  G[1]  =  G[8] = 1.;
  G[11] =  R;
  G[16] = -R;
}

extern "C" void h2y(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* y, unsigned int sizeZ, double* param)
{
  double epsilon2 = 1;
  y[0] = (q[1] - epsilon2 * R * sin(q[3]) * cos(q[4]) - R) + Wall;
  //cout << "h1=" << y[0] << endl;
}

extern "C" void G2y(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* G, unsigned int sizeZ, double* param)
{
  double epsilon2 = 1;

  G[3]  = 1.;
  G[9]  = -epsilon2 * R * cos(q[3]) * cos(q[4]);
  G[12] =  epsilon2 * R * sin(q[3]) * sin(q[4]);

  G[1]  =  G[8] = 1.;
  G[11] =  R;
  G[16] = -R;
}

extern "C" void h1y_(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* y, unsigned int sizeZ, double* param)
{
  double epsilon1 = -1;
  y[0] = Wall - (q[1] - epsilon1 * R * sin(q[3]) * cos(q[4])) - R;
  //cout << "h1=" << y[0] << endl;
}

extern "C" void G1y_(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* G, unsigned int sizeZ, double* param)
{
  double epsilon1 = -1;

  G[3]  = -1.;
  G[9]  =  epsilon1 * R * cos(q[3]) * cos(q[4]);
  G[12] = -epsilon1 * R * sin(q[3]) * sin(q[4]);

  G[1]  =  G[8] = 1.;
  G[11] =  R;
  G[16] = -R;
}

extern "C" void h2y_(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* y, unsigned int sizeZ, double* param)
{
  double epsilon2 = 1;
  y[0] = Wall - (q[1] - epsilon2 * R * sin(q[3]) * cos(q[4])) - R;
  //cout << "h1=" << y[0] << endl;
}

extern "C" void G2y_(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* G, unsigned int sizeZ, double* param)
{
  double epsilon2 = 1;

  G[3]  = -1.;
  G[9]  =  epsilon2 * R * cos(q[3]) * cos(q[4]);
  G[12] = -epsilon2 * R * sin(q[3]) * sin(q[4]);

  G[1]  =  G[8] = 1.;
  G[11] =  R;
  G[16] = -R;
}

/* Interactions with walls with respect to x axis
   We have two interactions for each couple concerning each bead and for each couple you have also interactions with right an left wall
*/

extern "C" void h1x(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* y, unsigned int sizeZ, double* param)
{
  double epsilon1 = -1;
  y[0] = q[0] + epsilon1 * R * sin(q[3]) * sin(q[4]) - R + Wall;
  //cout << "h1=" << y[0] << endl;
}

extern "C" void G1x(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* G, unsigned int sizeZ, double* param)
{
  double epsilon1 = -1;

  G[0]  = 1.;
  G[9]  =  epsilon1 * R * cos(q[3]) * sin(q[4]);
  G[12] =  epsilon1 * R * sin(q[3]) * cos(q[4]);

  G[4]  =  G[8] = 1.;
  G[14] =  R;
  G[16] = -R;
}

extern "C" void h2x(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* y, unsigned int sizeZ, double* param)
{
  double epsilon2 = 1;
  y[0] = q[0] + epsilon2 * R * sin(q[3]) * sin(q[4]) - R + Wall;
  //cout << "h1=" << y[0] << endl;
}

extern "C" void G2x(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* G, unsigned int sizeZ, double* param)
{
  double epsilon2 = 1;

  G[0]  = 1.;
  G[9]  =  epsilon2 * R * cos(q[3]) * sin(q[4]);
  G[12] =  epsilon2 * R * sin(q[3]) * cos(q[4]);

  G[4]  =  G[8] = 1.;
  G[14] =  R;
  G[16] = -R;
}

extern "C" void h1x_(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* y, unsigned int sizeZ, double* param)
{
  double epsilon1 = -1;
  y[0] = Wall - (q[0] + epsilon1 * R * sin(q[3]) * sin(q[4])) - R;
  //cout << "h1=" << y[0] << endl;
}

extern "C" void G1x_(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* G, unsigned int sizeZ, double* param)
{
  double epsilon1 = -1;

  G[0]  = -1.;
  G[9]  = -epsilon1 * R * cos(q[3]) * sin(q[4]);
  G[12] = -epsilon1 * R * sin(q[3]) * cos(q[4]);

  G[4]  =  G[8] = 1.;
  G[14] =  R;
  G[16] = -R;
}

extern "C" void h2x_(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* y, unsigned int sizeZ, double* param)
{
  double epsilon2 = 1;
  y[0] = Wall - (q[0] + epsilon2 * R * sin(q[3]) * sin(q[4])) - R;
  //cout << "h1=" << y[0] << endl;
}

extern "C" void G2x_(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* G, unsigned int sizeZ, double* param)
{
  double epsilon2 = 1;

  G[0]  = -1.;
  G[9]  = -epsilon2 * R * cos(q[3]) * sin(q[4]);
  G[12] = -epsilon2 * R * sin(q[3]) * cos(q[4]);

  G[4]  =  G[8] = 1.;
  G[14] =  R;
  G[16] = -R;
}


/* Interaction between couple of beads */

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
}

/* For contact between couples of beads*/

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

  G[18] = -G[0];
  G[21] = -G[3];
  G[24] = -G[6];


  G[11] =  epsilon1 * R;
  G[13] = -epsilon1 * R;

  G[29] = -epsilon2 * R;
  G[31] =  epsilon2 * R;

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


  G[18] = -G[0];
  G[21] = -G[3];
  G[24] = -G[6];


  G[11] =  epsilon1 * R;
  G[13] = -epsilon1 * R;

  G[29] = -epsilon2 * R;
  G[31] =  epsilon2 * R;
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


  G[18] = -G[0];
  G[21] = -G[3];
  G[24] = -G[6];


  G[11] =  epsilon1 * R;
  G[13] = -epsilon1 * R;

  G[29] = -epsilon2 * R;
  G[31] =  epsilon2 * R;
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


  G[18] = -G[0];
  G[21] = -G[3];
  G[24] = -G[6];


  G[11] =  epsilon1 * R;
  G[13] = -epsilon1 * R;

  G[29] = -epsilon2 * R;
  G[31] =  epsilon2 * R;
}
