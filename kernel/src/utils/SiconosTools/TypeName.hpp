/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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

/*! \file TypeName.hpp
  \brief get a string name from visitable classes 
*/

#ifndef TypeName_hpp
#define TypeName_hpp

#include "SiconosVisitor.hpp"
#include <string>
namespace Type
{
#undef REGISTER
#define REGISTER(X) case Type:: X : r.reset(new std::string(#X)); break;

#undef REGISTER_STRUCT
#define REGISTER_STRUCT(X) REGISTER(X)
#undef REGISTER_BASE
#undef REGISTER_BASE_EXTERN
#define REGISTER_BASE(X,Y) REGISTER(X)

#define REGISTER_BASE_EXTERN(X,Y) REGISTER_BASE(X,Y)

inline std::shared_ptr<std::string> str(const Siconos& X)
{
  std::shared_ptr<std::string> r;

  switch (X)
  {
    SICONOS_VISITABLES()
  default:
    assert(0);
  }

  return(r);
}


template <class C>
std::string name(const C& c)
{
  return *(Type::str(Type::value(c)));
}

}
#endif
