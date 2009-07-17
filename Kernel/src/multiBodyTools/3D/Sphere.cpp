/* Siconos-Example version 3.0.0, Copyright INRIA 2005-2008.
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
 * Foundation, Inc., 51 Franklin St, Fifth FLOOR, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 *
 */

#include "Sphere.h"

void Sphere::MassSetup()
{
  mass.reset(new PMMass(ndof, ndof));
  mass->resize(ndof, ndof);
  mass->zero();
  (*mass)(0, 0) = (*mass)(1, 1) = (*mass)(2, 2) = massValue;    ;
  (*mass)(3, 3) = (*mass)(4, 4) = (*mass)(5, 5) = 3. / 5 * massValue * radius * radius;
}

Sphere::Sphere(double r, double m,
               const SiconosVector& qinit,
               const SiconosVector& vinit)
  : LagrangianDS(qinit, vinit), radius(r), massValue(m)
{
  ndof = 6;
  MassSetup();
}

Sphere::~Sphere()
{}
