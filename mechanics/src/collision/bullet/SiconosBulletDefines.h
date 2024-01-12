/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

/*! \file BodyBulletShapeRecord.hpp
  \brief Body shape record
  The objective of this class is to keep associate ds and static body with the
shape in a contactor.
*/

#ifndef SICOBULLETDEFINES_H

#ifdef BULLET_TIMER
#define BT_ENABLE_PROFILE 1
#include <LinearMath/btQuickprof.h>
#endif

// Comment this to try un-queued static contactor behaviour
#define QUEUE_STATIC_CONTACTORS 1

// We can replace the primitives by alternative implementations.  To date,
// everything works (under test conditions, very tentative) except
// btStaticPlaneShape, so we replace it with a large box.

#define USE_CONVEXHULL_FOR_BOX 1
// #define USE_CONVEXHULL_FOR_SPHERE 1
// #define USE_BOX_FOR_PLANE 1
#define USE_CONVEXHULL_FOR_PLANE 1

// Bullet types (from USE flags above)
#ifdef USE_CONVEXHULL_FOR_BOX
#define BTBOXSHAPE btConvexHullShape
#else
#define BTBOXSHAPE btBoxShape
#endif

#ifdef USE_CONVEXHULL_FOR_SPHERE
#define BTSPHERESHAPE btConvexHullShape
#else
#define BTSPHERESHAPE btSphereShape
#endif

#ifdef USE_BOX_FOR_PLANE
#define BTPLANESHAPE btBoxShape
#else
#ifdef USE_CONVEXHULL_FOR_PLANE
#define BTPLANESHAPE btConvexHullShape
#else
#define BTPLANESHAPE btStaticPlaneShape
#endif
#endif


#endif
