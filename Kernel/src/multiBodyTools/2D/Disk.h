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

/*! \file Disk.h
  \brief Definition of a 2D disk - Inherits from LagrangianDS
*/

#ifndef Disk_H
#define Disk_H

#include "LagrangianDS.h"
#include "SiconosPointers.h"

class Disk : public LagrangianDS
{
private:
  double radiusDisk;
  double massDisk;
  unsigned int ndofDisk;

  boost::shared_ptr<SiconosVector> QDisk;
  boost::shared_ptr<SiconosVector> VDisk;
  boost::shared_ptr<SiconosVector> ADisk;

  void MassSetup();

  Disk();

public:

  /** Constructor
      \param radius
      \param mass
      \param x postion
      \param y position
  */

  Disk(double, double, double, double);

  /** Constructor
      \param number
      \param radius
      \param mass
      \param postion vector
      \param velocity vector
  */

  Disk(double, double, const SiconosVector&, const SiconosVector&);

  /** destructor
   */
  ~Disk();

  inline double getQ(unsigned int pos)
  {
    return (*q[0])(pos);
  };
  inline double getVelocity(unsigned int pos)
  {
    return (*q[1])(pos);
  };

  inline double getRadius()
  {
    return radiusDisk;
  };

  inline double getMassValue()
  {
    return massDisk;
  };

};

#endif /* Disk_H */

