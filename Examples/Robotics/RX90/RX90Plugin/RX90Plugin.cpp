/* Siconos-sample , Copyright INRIA 2005-2011.
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

#ifdef _WIN32 
#define SICONOS_EXPORT extern "C" __declspec(dllexport) 
#else 
#define SICONOS_EXPORT extern "C" 
#endif  
#include <iostream>
#include <math.h>

#define PI 3.14159265

extern "C" {
#include "LagrangianModel.h"
#include "ActuationModel.h"
}

SICONOS_EXPORT void mass(unsigned int sizeOfq, const double *q, double *mass, unsigned int sizeZ, double* z)
{
  double q_local[sizeOfq];
  for (unsigned int i = 0; i < sizeOfq; i++)
    q_local[i] = q[i];
  Inertia(mass, q_local);
}

SICONOS_EXPORT void FGyr(unsigned int sizeOfq, const double *q, const double *velocity, double *FGyr, unsigned int sizeZ, double* z)
{
  double q_local[sizeOfq];
  double v_local[sizeOfq];
  for (unsigned int i = 0; i < sizeOfq; i++)
  {
    q_local[i] = q[i];
    v_local[i] = velocity[i];
  }
  NLEffects(FGyr, q_local, v_local);
}

SICONOS_EXPORT void jacobianFGyrq(unsigned int sizeOfq, const double *q, const double *velocity, double *jacobq, unsigned int sizeOfZ, double* z)
{
  //dummy function
}

SICONOS_EXPORT void jacobianVFGyr(unsigned int sizeOfq, const double *q, const  double *velocity, double *jacobv, unsigned int sizeOfZ, double* z)
{
  //dummy function
}

SICONOS_EXPORT void U(double time, unsigned int sizeOfq, const double *q, const  double *velocity, double *U, unsigned int sizeZ, double* z)
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

SICONOS_EXPORT void jacobFintQ(double time, unsigned int sizeOfq, const double *q, const  double *velocity, double *jacobFintQ, unsigned int sizeOfZ, double* z)
{
  //dummy function
}
SICONOS_EXPORT void jacobFintV(double time, unsigned int sizeOfq, const double *q, const double *velocity, double *jacobFintV, unsigned int sizeOfZ, double* z)
{
  //dummy function
}
