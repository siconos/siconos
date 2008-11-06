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

/*! \file Sphere.h
  \brief Definition of a Sphere - Inherits from LagrangianDS
*/

#ifndef Sphere_H
#define Sphere_H

#include "LagrangianDS.h"

class LagrangianDS;

/** Sphere

   \author F. Perignon
   \version 3.0.0.
   \date (Creation) May 2008



 */
class Sphere : public LagrangianDS
{
private:

  // radius
  double Radius;

  // mass
  double mass;

  // Number of degrees of freedom
  unsigned int nDof;

  Sphere();

public:

  /** Constructor
      \param radius
      \param mass
      \param vector of initial positions
      \param vector of initial velocities
   */
  Sphere(double, double, const SiconosVector&, const SiconosVector&);

  /** destructor
   */
  ~Sphere();

  /* returns q(pos) */
  inline double getQ(unsigned int pos)
  {
    return (*q[0])(pos);
  };

  /* returns velocity(pos) */
  inline double getVelocity(unsigned int pos)
  {
    return (*q[1])(pos);
  };

  /* returns the radius of the sphere */
  inline double getRadius() const
  {
    return Radius;
  };

  /* draw the sphere: not yet implemented - See BeadsColumn example for a proper way of drawing spheres. */
  void draw();


};

TYPEDEF_SPTR(Sphere);

#endif
