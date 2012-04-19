/* Siconos-Numerics, Copyright INRIA 2005-2011.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
#ifndef NONSMOOTHNEWTONNEIGH_H
#define NONSMOOTHNEWTONNEIGH_H

#include "NonSmoothNewton.h"

/*!\file NonSmoothNewton.h
  Typedef and functions declarations related to non-smooth Newton solver

  Solve \f$ \phi(z) = 0 \f$ using a Newton method.

  The algorithm is alg 4.1 of the paper of Kanzow and Kleinmichel, "A new class of semismooth Newton-type methods
  for nonlinear complementarity problems", in Computational Optimization and Applications, 11, 227-251 (1998).

  \author Olivier Bonnefon

 */



#ifdef __cplusplus
extern "C"
{
#endif

  double* nonSmoothNewtonNeighInitMemory(int n, double * dWork, int * iWork);

  int nonSmoothNewtonNeigh(int, double*, NewtonFunctionPtr*, NewtonFunctionPtr*, int*, double*);

  int nonSmoothNewtonNeigh_getNbIWork(int n, int m);
  int nonSmoothNewtonNeigh_getNbDWork(int n, int m);

  /*only for debug*/
  void NSNN_thisIsTheSolution(int n, double * z);
  void  NSNN_reset();

#ifdef __cplusplus
}
#endif

#endif
