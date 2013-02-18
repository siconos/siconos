/* Siconos-Kernel, Copyright INRIA 2005-2012.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 *
 */

#include "Disk.hpp"

void Disk::MassSetup()
{
  _mass.reset(new SimpleMatrix(_ndof, _ndof));
  //  mass->resize(ndof,ndof);
  _mass->zero();
  (*_mass)(0, 0) = (*_mass)(1, 1) = massValue;
  (*_mass)(2, 2) = massValue * radius * radius / 2.;
}

Disk::Disk(double r, double m,
           SP::SiconosVector qinit,
           SP::SiconosVector vinit)
  : CircularDS(r, m, qinit, vinit)
{
  MassSetup();
}

Disk::~Disk()
{}
