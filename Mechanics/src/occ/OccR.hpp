/* Siconos-Kernel  Copyright INRIA 2005-2014.
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
 * Contact: Vincent ACARY siconos-team@lists.gforge.inria.fr
 */
/** \file OccR.hpp
    \brief A Siconos Newton Euler 3D friction relation between
    two BRep contact points.
 */

#ifndef OccR_hpp
#define OccR_hpp

#include "MechanicsFwd.hpp"
#include "NewtonEulerFrom3DLocalFrameR.hpp"

class OccR : public NewtonEulerFrom3DLocalFrameR
{
public:
  /** Constructor from contact points.
   *  \param contact1 : the first contact.
   *  \param contact2 : the second contact.
   */
  OccR(const ContactPoint& contact1, const ContactPoint& contact2);

  /** compute h.
   *  \param time : the time.
   *  \param q0 : the state vector.
   *  \param y : output vector.
   */
  virtual void computeh(double time, BlockVector& q0, SiconosVector& y);

protected:
  const ContactPoint& _contact1;
  const ContactPoint& _contact2;
  bool _normalFromFace1;
  bool _offsetp1;
  double _offset;
};

#endif
