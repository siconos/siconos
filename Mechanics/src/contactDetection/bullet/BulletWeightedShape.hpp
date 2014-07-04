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

#ifndef BULLET_WEIGHTED_SHAPE_HPP
#define BULLET_WEIGHTED_SHAPE_HPP

#include "BulletSiconosFwd.hpp"
#include <SiconosFwd.hpp>

class BulletWeightedShape
{
private:

  double _mass;
  SP::btCollisionShape _shape;
  SP::SimpleMatrix _inertia;

public:
  /* Constructor for an association of a bullet shape and an inertia
   * \param shape the collision shape
   * \param mass the mass
   */
  BulletWeightedShape(SP::btCollisionShape shape, const double& mass);

  /* Get the collision shape
   * \return a btCollisionShape
   */
  const SP::btCollisionShape& collisionShape() const
  {
    return _shape;
  };

  /* Get the mass
   * \return a double
   */
  const double& mass() const
  {
    return _mass;
  };

  /* Get the inertia matrix
   * \return a SP::SimpleMatrix
   */
  SP::SimpleMatrix inertia() const
  {
    return _inertia;
  };

  /* Modify the inertia matrix.
     \param newInertia the new inertia matrix
  */
  void setInertia(SP::SimpleMatrix newInertia)
  {
    _inertia = newInertia;
  }

  /* Modify the inertia matrix.
     \param ix x component
     \param iy y component
     \param iz z component
  */
  void setInertia(double ix, double iy, double iz);

};
#endif
