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

/*! \file SiconosContactor.hpp
  \brief Definition of an abstract contactor
*/


#ifndef SiconosContactor_h
#define SiconosContactor_h

#include <vector>
#include <utility>

#include "MechanicsFwd.hpp"

#include <SiconosSerialization.hpp>

#include "SiconosShape.hpp"

/** Class to hold the shape assigned to a body, and to associate each
 *  shape with an offset and collision group. */

class SiconosContactor
{
private:
  SiconosContactor() {};

protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(SiconosContactor);

public:
  SiconosContactor(SP::SiconosShape _shape,
                   SP::SiconosVector _offset = SP::SiconosVector(),
                   int _collision_group = 0);

  SP::SiconosShape shape;
  SP::SiconosVector offset;
  int collision_group;
};

class SiconosContactorSet : public std::vector< SP::SiconosContactor >
{
protected:
  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(SiconosContactorSet);

public:
  typedef std::vector< SP::SiconosContactor >::iterator iterator;

  void append(SP::SiconosContactor b) { push_back(b); }
  void append(std::vector<SP::SiconosContactor> b) { insert(end(), b.begin(), b.end()); }
  void append(const SiconosContactorSet& b) { insert(end(), b.begin(), b.end()); }
  void append(const SP::SiconosContactorSet& b) { insert(end(), b->begin(), b->end()); }
};

#endif /* SiconosContactor_h */
