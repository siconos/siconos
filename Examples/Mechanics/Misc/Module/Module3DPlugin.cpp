/* Siconos version 1.0, Copyright INRIA 2005.
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


const double m = 0.2; // bar2D mass
const double g = 9.8; // gravity


extern "C" double FextFunction(double time)
{
  double res = -0.0;
  return res;
}


extern "C" void module3DFExt(unsigned int sizeOfq, double time, double *fExt, double* param)
{
  /* input parameter : sizeOfq (size of the vector q); time; q (pointer to q vector);
   * output parameter : fExt (pointer to Fext vector)
   */

  int i;
  int n = sizeOfq;
  double T = 0.3;


  for (i = 0; i < n; i++)
  {

    fExt[i] =  FextFunction(time) + param[i] * time / T;

  }

  //        fExt[0] = -20 + fExt[0];

}

extern "C" void groundFExt(unsigned int sizeOfq, double time, double *fExt, double* param)
{
  /* input parameter : sizeOfq (size of the vector q); time; q (pointer to q vector);
   * output patarmeter : fExt (pointer to Fext vector)
   */


  int i;
  int n = sizeOfq;

  for (i = 0; i < n; i++)
    fExt[i] = 0.0;
}

