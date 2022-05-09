/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

/*! \file SiconosSerialization.hpp
  serialization for Siconos
*/

#ifndef SiconosSerialization_hpp
#define SiconosSerialization_hpp

// tells include-what-you-use to keep this file
// and not to suggest boost or alike.
// IWYU pragma: begin_exports 

namespace boost
{
namespace serialization
{
class access;
}
}

/** install serialization hooks. Must be used inside a protected zone
    of class definition
    \parameter a class name
*/
#define ACCEPT_SERIALIZATION(CLASS)                             \
  typedef void serializable;                                    \
  template<typename Archive>                                    \
  friend void siconos_io(Archive&, CLASS&, const unsigned int); \
  friend class boost::serialization::access

// IWYU pragma: end_exports

#endif
