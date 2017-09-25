/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2017 INRIA.
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

#ifndef ContactR_hpp
#define ContactR_hpp

#include "MechanicsFwd.hpp"
#include "SiconosVector.hpp"
#include "NewtonEulerFrom3DLocalFrameR.hpp"

class ContactR : public NewtonEulerFrom3DLocalFrameR
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(ContactR);

  double _y_correction_A;
  double _y_correction_B;
  double _scaling;

  ContactR()
    : _y_correction_A(0), _y_correction_B(0), _scaling(1),
      _flip(false), _contactDistance(0)
    {}

public:
  ContactR(SP::SiconosVector q1, SP::SiconosVector q2,
           bool flip=false,
           double y_correction_A=0,
           double y_correction_B=0,
           double scaling=1);

  bool _flip; /* if true, points A and B are swapped from the point
               * view of the Relation. */

  double _contactDistance;

  /* For users that may require extra information about contacts. */
  SP::SiconosVector base[2];
  SP::SiconosShape shape[2];
  SP::SiconosContactor contactor[2];
  SP::BodyDS ds[2];

  double distance() const
  {
    return _contactDistance;
  };

  double y_correction_A() { return _y_correction_A; }
  double y_correction_B() { return _y_correction_A; }
  double y_correction() { return _y_correction_A + _y_correction_B; }

  virtual void computeh(double time, BlockVector& q0, SiconosVector& y);

  virtual void updateContactPoints(SiconosVector& pos1, SiconosVector& pos2,
                                   double distance, SiconosVector& normal,
                                   SP::NewtonEulerDS ds1, SP::NewtonEulerDS ds2);

  virtual void preDelete() {}

  ACCEPT_STD_VISITORS();
};

#endif
