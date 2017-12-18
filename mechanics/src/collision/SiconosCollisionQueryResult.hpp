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

/*! \file SiconosCollisionQueryResult.hpp
\brief Holds one result of a query against the graph of body
contactors maintained by a SiconosCollisionManager.
*/

#ifndef SiconosCollisionQueryResult_h
#define SiconosCollisionQueryResult_h

#include <MechanicsFwd.hpp>

/** \brief Holds one result of a line segment intersection query
    against the graph of body contactors maintained by a
    SiconosCollisionManager. */
class SiconosCollisionQueryResult
{
protected:
  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(SiconosCollisionQueryResult);

public:

  /** Distance from reference point (start of line segment or query center) */
  double distance;

  /** Body owning the contactor that was intersected, may be null for
   * static contactors. */
  SP::BodyDS body;

  /** The shape that was intersected. */
  SP::SiconosShape shape;

  /** The contactor that was intersected. */
  SP::SiconosContactor contactor;

  /** Closest point on contactor in world coordinates. */
  SiconosVector point;
};

#endif /* SiconosCollisionQueryResult.hpp */
