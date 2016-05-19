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


#ifndef BulletSiconosFwd_hpp
#define BulletSiconosFwd_hpp

#include "SiconosPointers.hpp"

DEFINE_SPTR(btCollisionShape);
DEFINE_SPTR(btBoxShape);
DEFINE_SPTR(btCylinderShape);
DEFINE_SPTR(btManifoldPoint);
DEFINE_SPTR(btStaticPlaneShape);
DEFINE_SPTR(btSphereShape);
DEFINE_SPTR(btConvexHullShape);

DEFINE_SPTR(btTriangleIndexVertexMaterialArray);

DEFINE_SPTR(btCollisionObject);
DEFINE_SAPTR(btCollisionObject);

DEFINE_SPTR(btCollisionWorld);
DEFINE_SPTR(btDefaultCollisionConfiguration);
DEFINE_SPTR(btCollisionDispatcher);
DEFINE_SPTR(btVector3);
DEFINE_SPTR(btQuaternion);

DEFINE_SPTR(btPersistentManifold);

DEFINE_SPTR(btBroadphaseInterface);

//DEFINE_SPTR(BulletDS);
DEFINE_SPTR(BulletWeightedShape);
//DEFINE_SPTR(BulletR);
//DEFINE_SPTR(BulletFrom1DLocalFrameR);
//DEFINE_SPTR(BulletSpaceFilter);
//DEFINE_SPTR(BulletTimeStepping);
DEFINE_SPTR(CollisionObjects);
DEFINE_SPTR(StaticObjects);

#include "MechanicsFwd.hpp"

#endif
