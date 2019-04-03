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

#ifndef Bullet1DR_hpp
#define Bullet1DR_hpp

#include "BulletSiconosFwd.hpp"
#include <NewtonEuler1DR.hpp>

class Bullet1DR : public NewtonEuler1DR
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(Bullet1DR);

  SP::btManifoldPoint _contactPoints;

public:
  Bullet1DR(SP::btManifoldPoint);

  SP::btManifoldPoint contactPoint() const
  {
    return _contactPoints;
  };

  virtual void computeh(double time, BlockVector&, SiconosVector&);

  ACCEPT_STD_VISITORS();
};
#endif
