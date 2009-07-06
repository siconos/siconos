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

#include "CircularDS.h"

class Disk : public CircularDS, public boost::enable_shared_from_this<Disk>
{
private:

  void MassSetup();

protected:
  Disk();

public:

  /** Constructor
      \param radius
      \param mass
      \param postion vector
      \param velocity vector
  */

  Disk(double, double, const SiconosVector&, const SiconosVector&);

  /** destructor
   */
  ~Disk();


  /** visitors hook
   */
  virtual void accept(SiconosVisitor& tourist)
  {
    tourist.visit(*this);
  }
  virtual void accept(SP::SiconosVisitor tourist)
  {
    tourist->visit(shared_from_this());
  }

};

TYPEDEF_SPTR(Disk);

#endif /* Disk_H */

