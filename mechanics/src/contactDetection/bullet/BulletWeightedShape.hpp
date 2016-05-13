/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
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
