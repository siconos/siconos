/* Siconos-Kernel, Copyright INRIA 2005-2011.
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

/*! \file CircularDS
  \brief Definition of a 2D circular shape - Inherits from LagrangianDS
*/


#ifndef CircularDS_h
#define CircularDS_h

#include "LagrangianDS.hpp"

class CircularDS : public LagrangianDS
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(CircularDS);

  double radius;
  double massValue;

  CircularDS(): LagrangianDS() {};

public:

  CircularDS(double, double, SP::SiconosVector, SP::SiconosVector);

  ~CircularDS();

  inline double getQ(unsigned int pos)
  {
    assert(pos < _ndof);
    return (*_q[0])(pos);
  };
  inline double getVelocity(unsigned int pos)
  {
    assert(pos < _ndof);
    return (*_q[1])(pos);
  };

  inline double getMassValue() const
  {
    return massValue;
  };

  inline double getRadius() const
  {
    return radius;
  };

};

TYPEDEF_SPTR(CircularDS);

#endif /* CircularDS_h */
