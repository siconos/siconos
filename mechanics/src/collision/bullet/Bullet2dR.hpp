/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

#ifndef Bullet2dR_hpp
#define Bullet2dR_hpp

#include "BulletSiconosFwd.hpp"
#include "Contact2dR.hpp"

class Bullet2dR : public Contact2dR
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(Bullet2dR);

public:
  Bullet2dR();

  virtual ~Bullet2dR() {}

  /* For users that may require extra information about contacts. */
  SP::btCollisionObject btObject[2];
  SP::btCollisionShape btShape[2];

  virtual
  void updateContactPointsFromManifoldPoint(const btPersistentManifold& manifold,
                                            const btManifoldPoint& point,
                                            bool flip, double scaling,
                                            SP::RigidBody2dDS ds1,
                                            SP::RigidBody2dDS ds2);

  ACCEPT_STD_VISITORS();
};

#endif
