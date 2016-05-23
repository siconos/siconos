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

/*! \file SphereNEDS.hpp

  \brief Definition of a 3D Sphere as a NewtonEulerDS (with
  quaternions).

*/
#ifndef SphereNEDS_h
#define SphereNEDS_h

#include <MechanicsFwd.hpp>
#include "NewtonEulerDS.hpp"


class SphereNEDS : public NewtonEulerDS, public std11::enable_shared_from_this<SphereNEDS>
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(SphereNEDS);

  double radius;

  SphereNEDS() {};

public:

  SphereNEDS(double, double, SP::SiconosMatrix, SP::SiconosVector, SP::SiconosVector);

  ~SphereNEDS();

  inline double getQ(unsigned int pos)
  {
    assert(pos < 7);
    return (_q->getValue(pos));
  };

  inline double getVelocity(unsigned int pos)
  {
    assert(pos < 6);
    return (_v->getValue(pos));
  };

  inline double getRadius() const
  {
    return radius;
  };

  /** visitors hook
   */
  ACCEPT_BASE_VISITORS(NewtonEulerDS);

};
#endif /* SphereNEDS_h */
