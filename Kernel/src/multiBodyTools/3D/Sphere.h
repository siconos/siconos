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

/*! \file Sphere
  \brief Definition of a 2D circular shape - Inherits from LagrangianDS
*/


#ifndef Sphere_h
#define Sphere_h

#include "LagrangianDS.h"

class Sphere : public LagrangianDS, public boost::enable_shared_from_this<Sphere>
{
protected:
  double radius;
  double massValue;
  double I;

public:

  Sphere(double, double, const SiconosVector&, const SiconosVector&);

  ~Sphere();

  inline double getQ(unsigned int pos)
  {
    assert(pos < ndof);
    return (*q[0])(pos);
  };
  inline double getVelocity(unsigned int pos)
  {
    assert(pos < ndof);
    return (*q[1])(pos);
  };

  inline double getMassValue()
  {
    return massValue;
  };

  inline double getRadius()
  {
    return radius;
  };

  void computeMass();

  void computeNNL(SP::SiconosVector, SP::SiconosVector);

  void computeNNL();

  void computeJacobianNNL(unsigned int);

  void computeJacobianNNL(unsigned int, SP::SiconosVector, SP::SiconosVector);

  /** visitors hook
   */
  ACCEPT_VISITORS();

};

TYPEDEF_SPTR(Sphere);

#endif /* Sphere_h */
