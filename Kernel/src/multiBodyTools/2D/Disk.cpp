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

#include "Disk.h"

void Disk::MassSetup()
{
  MassDisk.reset(new SimpleMatrix(ndof, ndof));
  MassDisk->resize(ndof, ndof);
  MassDisk->zero();
  (*MassDisk)(0, 0) = (*MassDisk)(1, 1) = massDisk;
  (*MassDisk)(2, 2) = massDisk * radiusDisk * radiusDisk / 2.;
  setMassPtr(MassDisk.get());
}

Disk::Disk(int number,
           double r, double m,
           double x, double y)
  : LagrangianDS(), radiusDisk(r), massDisk(m), ndofDisk(3)
{

  QDisk.reset(new SimpleVector(ndofDisk));
  VDisk.reset(new SimpleVector(ndofDisk));
  ADisk.reset(new SimpleVector(ndofDisk));

  setNdof(ndofDisk);
  setNumber(number);

  q.resize(3, NULL);

  QDisk->zero();
  VDisk->zero();
  ADisk->zero();

  QDisk->setValue(0, x);
  QDisk->setValue(1, y);
  setQPtr(QDisk.get());
  setVelocity(*VDisk);

  setQ0(*QDisk);

  setVelocity0(*VDisk);

  // Missing setAccelerationPtr()
  q[2] = new SimpleVector(ndofDisk);
  isAllocatedIn["acceleration"] = true;
  *q[2] = *ADisk;

  MassSetup();
}


Disk::Disk(int number, double r, double m,
           const SiconosVector& qinit,
           const SiconosVector& vinit)
  : LagrangianDS(number, qinit, vinit),
    radiusDisk(r), massDisk(m), ndofDisk(3)
{
  MassSetup();
}

Disk::~Disk()
{}
