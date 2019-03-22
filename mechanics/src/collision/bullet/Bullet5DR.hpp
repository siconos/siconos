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

#ifndef Bullet5DR_hpp
#define Bullet5DR_hpp

#include "BulletSiconosFwd.hpp"
#include "Contact5DR.hpp"

class Bullet5DR : public Contact5DR
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(Bullet5DR);

public:
  Bullet5DR();

  virtual ~Bullet5DR() {}

  /* For users that may require extra information about contacts. */
  SP::btCollisionObject btObject[2];
  SP::btCollisionShape btShape[2];

  virtual
  void updateContactPointsFromManifoldPoint(const btPersistentManifold& manifold,
                                            const btManifoldPoint& point,
                                            bool flip, double scaling,
                                            SP::NewtonEulerDS ds1,
                                            SP::NewtonEulerDS ds2);

  ACCEPT_STD_VISITORS();
};

#endif
