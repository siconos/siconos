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

#ifndef BulletR_hpp
#define BulletR_hpp

#include "BulletSiconosFwd.hpp"
#include "ContactR.hpp"

class BulletR : public ContactR
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(BulletR);

  SP::btManifoldPoint _contactPoints;

public:
  BulletR(const btManifoldPoint &,
          SP::SiconosVector q1, SP::SiconosVector q2,
          bool flip=false,
          double y_correction_A=0,
          double y_correction_B=0,
          double scaling=1);

  /* For users that may require extra information about contacts. */
  SP::btCollisionObject btObject[2];
  SP::btCollisionShape btShape[2];

  // TODO used by BulletSpaceFilter
  SP::btManifoldPoint contactPoint() const
  {
    return _contactPoints;
  };

  // TODO used by BulletSpaceFilter
  void setContactPoint(SP::btManifoldPoint p)
  {
    _contactPoints = p;
  };

  virtual void updateContactPoints(const btManifoldPoint& point,
                                   SP::NewtonEulerDS ds1, SP::NewtonEulerDS ds2);

  ACCEPT_STD_VISITORS();
};

#endif
