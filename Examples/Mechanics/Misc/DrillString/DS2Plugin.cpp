/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2007.
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

#include <stdio.h>
#include <math.h>
const double km = 4.3228;
const double u = 2.5;
const double DTsu = -0.00575;
const double Dbu  = -0.0084;
const double Tsu = 0.37975;
const double Tsl = 0.2781;
const double Tcl = 0.0473;
const double omegasl = 1.4302;
const double deltasl = 2.0575;
const double bu = 2.4245;
const double kth = 0.0775;
const double bl = 0.0105;
const double Ju = 0.4765;
const double Jl = 0.0414;



// ===== Lagrangian Dynamical System =====

// Plugins for Fext, Fint, NNL (vectors), Mass, JacobianQNNL, JacobianVelocityNNL,
// JacobianQFint and JacobianVelocityFint (matrices)
extern "C" void FExt(double time, unsigned int sizeOfq, double *fExt, unsigned int sizeZ, double* z)
{
  fExt[0] = km * u;
  fExt[1] = 0;
  fExt[2] = 0;
  fExt[3] = 0;
}

extern "C" void FInt(double time, unsigned int sizeOfq, const double *q, const double *vel, double *fInt, unsigned int sizeZ, double * z)
{
  fInt[0] = DTsu + Dbu * fabs(vel[0]);
  fInt[1] = 0;
  fInt[2] = Tsu;
  fInt[3] = ((Tsl - Tcl) * exp(-pow(fabs(vel[1] / omegasl), deltasl)) + Tcl);
}

extern "C" void NNL(unsigned int sizeOfq, const double *q, const double *vel, double *NNL, unsigned int sizeZ, double *z)
{
  NNL[0] = bu * vel[0] + kth * q[0] - kth * q[1];
  NNL[1] = bl * vel[1] - kth * q[0] + kth * q[1];
  NNL[2] = 0.0;
  NNL[3] = 0.0;
}

