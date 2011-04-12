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

#include "BulletDS.hpp"

#include <BulletCollision/CollisionDispatch/btCollisionObject.h>
#include <BulletCollision/CollisionShapes/btBoxShape.h>
#include <BulletCollision/CollisionShapes/btCylinderShape.h>
#include <BulletCollision/CollisionShapes/btCapsuleShape.h>
#include <LinearMath/btVector3.h>

BulletDS::BulletDS(const SP::btCollisionShape& shape,
                   SP::SiconosVector position,
                   SP::SiconosVector velocity,
                   const double& mass) : NewtonEulerDS(), _collisionShape(shape)
{
  _mass = mass;

  btVector3 inertia;

  SP::SimpleMatrix inertiaMatrix(new SimpleMatrix(3, 3));
  inertiaMatrix->eye();

  internalInit(position, velocity, mass, inertiaMatrix);

  _collisionObject.reset(new btCollisionObject());

  btMatrix3x3 basis;
  basis.setIdentity();
  _collisionObject->getWorldTransform().setBasis(basis);

  _collisionObject->setUserPointer(this);
  _collisionObject->setCollisionFlags(_collisionObject->getCollisionFlags() | btCollisionObject::CF_KINEMATIC_OBJECT);

  _collisionShape->calculateLocalInertia(mass, inertia);
  (*_I)(0, 0) = inertia[0];
  (*_I)(1, 1) = inertia[1];
  (*_I)(2, 2) = inertia[2];

  _collisionObject->setCollisionShape(&*_collisionShape);
}


BulletDS::BulletDS(const BroadphaseNativeTypes& shape_type,
                   const SP::SimpleVector& shapeParams,
                   SP::SiconosVector position,
                   SP::SiconosVector velocity,
                   const double& mass) : NewtonEulerDS()
{
  _mass = mass;

  btVector3 inertia;

  SP::SimpleMatrix inertiaMatrix(new SimpleMatrix(3, 3));
  inertiaMatrix->eye();

  internalInit(position, velocity, mass, inertiaMatrix);

  _collisionObject.reset(new btCollisionObject());
  btMatrix3x3 basis;
  basis.setIdentity();
  _collisionObject->getWorldTransform().setBasis(basis);

  _collisionObject->setUserPointer(this);
  _collisionObject->setCollisionFlags(_collisionObject->getCollisionFlags() | btCollisionObject::CF_KINEMATIC_OBJECT);


  switch (shape_type)
  {
  case BOX_SHAPE_PROXYTYPE :
  {
    btVector3 params;
    params[0] = (*shapeParams)(0);
    params[1] = (*shapeParams)(1);
    params[2] = (*shapeParams)(2);
    _collisionShape.reset(new btBoxShape(params));
    _collisionShape->calculateLocalInertia(mass, inertia);
    (*_I)(0, 0) = inertia[0];
    (*_I)(1, 1) = inertia[1];
    (*_I)(2, 2) = inertia[2];

    _collisionObject->setCollisionShape(&*_collisionShape);

    break;
  };
  case CYLINDER_SHAPE_PROXYTYPE :
  {
    btVector3 params;
    params[0] = (*shapeParams)(0);
    params[1] = (*shapeParams)(1);
    params[2] = (*shapeParams)(2);
    _collisionShape.reset(new btCylinderShape(params));
    _collisionShape->calculateLocalInertia(mass, inertia);
    (*_I)(0, 0) = inertia[0];
    (*_I)(1, 1) = inertia[1];
    (*_I)(2, 2) = inertia[2];

    _collisionObject->setCollisionShape(&*_collisionShape);

    break;
  };

  case CAPSULE_SHAPE_PROXYTYPE :
  {
    btScalar radius;
    btScalar height;
    radius = (*shapeParams)(0);
    height = (*shapeParams)(1);
    _collisionShape.reset(new btCapsuleShape(radius, height));
    _collisionShape->calculateLocalInertia(mass, inertia);
    (*_I)(0, 0) = inertia[0];
    (*_I)(1, 1) = inertia[1];
    (*_I)(2, 2) = inertia[2];

    _collisionObject->setCollisionShape(&*_collisionShape);

    break;
  };
  default :
  {
    RuntimeException::selfThrow("BulletDS: unknown shape");
  };
  }
}
