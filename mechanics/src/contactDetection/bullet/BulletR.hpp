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
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
*/

#ifndef BulletR_hpp
#define BulletR_hpp

#include "BulletSiconosFwd.hpp"
#include "NewtonEulerFrom3DLocalFrameR.hpp"

class BulletR : public NewtonEulerFrom3DLocalFrameR
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(BulletR);

  const SP::btManifoldPoint _contactPoints;

  const double _y_correction_A;
  const double _y_correction_B;
  const double _scaling;

public:
  BulletR(SP::btManifoldPoint,
          double y_correction_A=0,
          double y_correction_B=0,
          double scaling=1);

  SP::btManifoldPoint contactPoint() const
  {
    return _contactPoints;
  };

  double y_correction_A() { return _y_correction_A; }
  double y_correction_B() { return _y_correction_A; }
  double y_correction() { return _y_correction_A + _y_correction_B; }

  virtual void computeh(double time, BlockVector& q0, SiconosVector& y);

  ACCEPT_STD_VISITORS();
};

#endif
