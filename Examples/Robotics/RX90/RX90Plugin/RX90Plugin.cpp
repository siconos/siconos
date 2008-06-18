/* Siconos-sample version 3.0.0, Copyright INRIA 2005-2008.
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

#include <iostream>
#include <math.h>

#define PI 3.14159265

extern "C" {
#include "LagrangianModel.h"
#include "ActuationModel.h"
}

extern "C" void mass(unsigned int sizeOfq, const double *q, double *mass, unsigned int sizeZ, double* z)
{
  double q_local[sizeOfq];
  for (unsigned int i = 0; i < sizeOfq; i++)
    q_local[i] = q[i];
  Inertia(mass, q_local);
}

extern "C" void NNL(unsigned int sizeOfq, const double *q, const double *velocity, double *NNL, unsigned int sizeZ, double* z)
{
  double q_local[sizeOfq];
  double v_local[sizeOfq];
  for (unsigned int i = 0; i < sizeOfq; i++)
  {
    q_local[i] = q[i];
    v_local[i] = velocity[i];
  }
  NLEffects(NNL, q_local, v_local);
}

extern "C" void jacobianQNNL(unsigned int sizeOfq, const double *q, const double *velocity, double *jacobq, unsigned int sizeOfZ, double* z)
{
  //dummy function
}

extern "C" void jacobianVNNL(unsigned int sizeOfq, const double *q, const  double *velocity, double *jacobv, unsigned int sizeOfZ, double* z)
{
  //dummy function
}

extern "C" void U(double time, unsigned int sizeOfq, const double *q, const  double *velocity, double *U, unsigned int sizeZ, double* z)
{
  int ndof = sizeOfq;
  int ncont = 2;
  int z_size = 0;
  double *state = 0;
  double *myz = 0;
  double *zdot = 0;
  double q_local[sizeOfq];
  double v_local[sizeOfq];
  for (unsigned int i = 0; i < sizeOfq; i++)
  {
    q_local[i] = q[i];
    v_local[i] = velocity[i];
  }
  actuationDynamics(&time, q_local, v_local, myz, state, &ndof, &ncont, &z_size, zdot, U);
  SpringForce(v_local, q_local);
  for (unsigned int i = 0; i < sizeOfq; i++)
    U[i] = -U[i] + v_local[i];
}

extern "C" void jacobFintQ(double time, unsigned int sizeOfq, const double *q, const  double *velocity, double *jacobFintQ, unsigned int sizeOfZ, double* z)
{
  //dummy function
}
extern "C" void jacobFintV(double time, unsigned int sizeOfq, const double *q, const double *velocity, double *jacobFintV, unsigned int sizeOfZ, double* z)
{
  //dummy function
}
