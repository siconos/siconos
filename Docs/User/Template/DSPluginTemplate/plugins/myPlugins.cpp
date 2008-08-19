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
#include <stdio.h>

extern "C" void myF(double time, unsigned int sizeX, const double* x, double* f, unsigned int sizeZ, double* z)
{
  // Do your stuff ...

  printf("call myF ... \n");

}

extern "C" void myJacobianF(double time, unsigned int sizeX, const double* x, double* jacobF, unsigned int sizeZ, double* z)
{
  // Do your stuff ...

  printf("call myJacobianF ... \n ");

}

extern "C" void myH(unsigned int sizeX, const double* x, unsigned int sizeY, double* Y, unsigned int sizeZ, double* z)
{
  // Do your stuff
  printf("call myH ... \n ");
}

extern "C" void myG(unsigned int sizeY, const double* lambda, unsigned int sizeX, double* g, unsigned int sizeZ, double* z)
{
  // Do your stuff
  printf("call myG ... \n ");
}

extern "C" void myJacobianH(unsigned int sizeX, const double* x, unsigned int sizeY, double* jacob, unsigned int sizeZ, double* z)
{
  // Do your stuff
  printf("call myJacobianH ... \n ");
}

extern "C" void myJacobianG(unsigned int sizeY, const double* lambda, unsigned int sizeX, double* jacob, unsigned int sizeZ, double* z)
{
  // Do your stuff
  printf("call myJacobianG ... \n ");
}

