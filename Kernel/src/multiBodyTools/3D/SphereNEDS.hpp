/* Siconos-Kernel, Copyright INRIA 2005-2010.
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

/*! \file SphereNEDS

  \brief Definition of a 3D Sphere as a NewtonEulerDS (with
  quaternions).

*/
#ifndef SphereNEDS_h
#define SphereNEDS_h

#include "NewtonEulerDS.hpp"


class SphereNEDS : public NewtonEulerDS, public boost::enable_shared_from_this<SphereNEDS>
{
protected:
  double radius;

public:

  SphereNEDS(double, double, SP::SiconosMatrix, SP::SiconosVector, SP::SiconosVector);

  ~SphereNEDS();

  inline double getQ(unsigned int pos)
  {
    assert(pos < 7);
    return (_q->getValue(pos));
  };

  inline double getVelocity(unsigned int pos)
  {
    assert(pos < 6);
    return (_v->getValue(pos));
  };

  inline double getMassValue()
  {
    return _mass;
  };

  inline double getRadius()
  {
    return radius;
  };

  /** visitors hook
   */
  ACCEPT_SP_VISITORS();

};

TYPEDEF_SPTR(SphereNEDS);

#endif /* SphereNEDS_h */
