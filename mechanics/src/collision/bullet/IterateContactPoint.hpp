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

#ifndef ITERCONTACTPOINT_HPP
#define ITERCONTACTPOINT_HPP

/*! \file IterateContactpoints.hpp tools to iterate over contact points in a collision world
  Internal stuff for Bullet collision manager.
 */

#include <BulletCollision/CollisionDispatch/btCollisionWorld.h>

#include <memory>

namespace siconos::collision::bullet::internal {

/** This class allows to iterate over all the contact points in a
 *  btCollisionWorld, returning a tuple containing the two btCollisionObjects
 *  and the btManifoldPoint.  To be called after performDiscreteCollisionDetection().
 */
class IterateContactPoints {
 public:
  std::shared_ptr<btCollisionWorld> world{nullptr};

  IterateContactPoints(std::shared_ptr<btCollisionWorld> wld) : world(wld) {}

  struct ContactPointTuple {
    const btCollisionObject *objectA;
    const btCollisionObject *objectB;
    btManifoldPoint *point;
    btPersistentManifold *manifold;
  };

  class iterator {
   protected:
    std::shared_ptr<btCollisionWorld> world{nullptr};
    ContactPointTuple data;
    int numManifolds{0};
    int numContacts{0};
    int manifold_index{-1};
    int contact_index{-1};
    
   public:
    iterator(std::shared_ptr<btCollisionWorld> wld)
        : world{wld}, numManifolds{wld->getDispatcher()->getNumManifolds()} {
      ++(*this);
    }

    iterator() = default;

    ~iterator() noexcept = default;
    
    const ContactPointTuple &operator*() { return data; };
    const ContactPointTuple *operator->() { return &data; };

    iterator &operator++() {
      if (numManifolds == 0) return *this;
      contact_index++;
      while (contact_index >= numContacts) {
        manifold_index++;
        if (manifold_index < numManifolds) {
          data.manifold = world->getDispatcher()->getManifoldByIndexInternal(manifold_index);
          // display_info_manifold(*data.manifold);
          data.objectA = data.manifold->getBody0();
          data.objectB = data.manifold->getBody1();
          numContacts = data.manifold->getNumContacts();
          contact_index = 0;
        } else {
          numManifolds = 0;
          return *this;
        }
      }
      data.point = &(data.manifold->getContactPoint(contact_index));
      return *this;
    };

    bool operator!=(const iterator &it) {
      if (it.numManifolds == 0) return numManifolds != 0;
      return data.objectA != it.data.objectA || data.objectB != it.data.objectB ||
             data.point != it.data.point;
    };
  };

  iterator begin() { return iterator(world); };

  iterator end() { return iterator(); };
};
}  // namespace siconos::collision::bullet::internal

#endif
