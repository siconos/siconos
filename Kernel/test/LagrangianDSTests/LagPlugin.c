/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
#include<iostream>

// ===== Lagrangian DS  =====

// Plugins for Fext, FInt, QNLInertia (vectors), Mass, JacobianQNLInertiaq, JacobianQNLInertiaVelocity,
// JacobianFIntq and JacobianFintVelocity (matrices)

extern "C" void computeFInt(unsigned int *sizeOfq, const double *time, double *q, double *velocity, double *fInt)
{
  unsigned int n = *sizeOfq;
  for (unsigned int i = 0; i < n; i++)
    fInt[i] = *time * q[i];
}

extern "C" void computeFExt(unsigned int *sizeOfq, const double *time, double *fExt)
{
  unsigned int n = *sizeOfq;
  for (unsigned int i = 0; i < n; i++)
    fExt[i] = 2 * *time;
}

extern "C" void computeQNLInertia(unsigned int *sizeOfq, double *q, double *velocity, double *Q)
{
  unsigned int n = *sizeOfq;
  for (unsigned int i = 0; i < n; i++)
    Q[i] = q[i];

}


extern "C" void computeMass(unsigned int *sizeOfq, const double *time, double *q, double *mass)
{
  /* input parameter : sizeOfq (size of the vector q); time ; q (pointer to q vector);
   * output parameter : mass (pointer to mass matrix)
   */

  unsigned int n = *sizeOfq;
  for (unsigned int i = 0; i < n; i++)
  {
    for (unsigned int j = 0; j < n; j++)
      mass[i + j * n] = q[i];
  }

}


extern "C" void computeJacobianFIntq(unsigned int *sizeOfq, const double *time, double *q, double *velocity, double *jacob)
{

  unsigned int n = *sizeOfq;
  for (unsigned int i = 0; i < n; i++)
  {
    for (unsigned int j = 0; j < n; j++)
      jacob[i + j * n] = *time * velocity[i];
  }

}

extern "C" void computeJacobianFintVelocity(unsigned int *sizeOfq, const double *time, double *q, double *velocity, double *jacob)
{
  unsigned int n = *sizeOfq;
  for (unsigned int i = 0; i < n; i++)
  {
    for (unsigned int j = 0; j < n; j++)
      jacob[i + j * n] = *time * velocity[i];
  }

}

extern "C" void computeJacobianQNLInertiaq(unsigned int *sizeOfq, double *q, double *velocity, double *jacob)
{
  unsigned int n = *sizeOfq;
  for (unsigned int i = 0; i < n; i++)
  {
    for (unsigned int j = 0; j < n; j++)
      jacob[i + j * n] = 2 * velocity[i];
  }
}

extern "C" void computeJacobianQNLInertiaVelocity(unsigned int *sizeOfq, double *q, double *velocity, double *jacob)
{
  unsigned int n = *sizeOfq;
  for (unsigned int i = 0; i < n; i++)
  {
    for (unsigned int j = 0; j < n; j++)
      jacob[i + j * n] = velocity[i];
  }

}

