/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef BulletUtils_hpp
#define BulletUtils_hpp

#include "BulletDeclarations.h"

namespace siconos::collision::bullet {
void display_info_collision_object(const btCollisionObject* collisionObject);
void display_info_manifold(const btPersistentManifold& manifold);
void display_info_contact_point(const btManifoldPoint& cp);
void display_info_btConvexHullShape(const btConvexHullShape& cp);
}  // namespace siconos::collision::bullet
#endif
